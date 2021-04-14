#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "GEM_cosmic_tracks.C"
#include "TString.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TSystem.h"
#include <iostream>
#include <fstream> 
#include "TROOT.h"
#include "TClonesArray.h"
#include "TVectorD.h"
#include "TMatrixD.h"
//#include "Math/Functor.h"
#include "TMinuit.h"

vector<double> AsymALL, dAsymALL;
vector<int> APVX_asymALL, APVY_asymALL;
vector<vector<double> > weightX, weightY;

vector<int> countX,countY;

int nAPVmaxX = 12;
int nAPVmaxY = 10;

void chi2_FCN( int &npar, double *gin, double &f, double *par, int flag ){

  //cout << "calculating chi2..." << endl;
  
  double chi2 = 0.0;
  for( int i=0; i<AsymALL.size(); i++ ){
    int iy = APVY_asymALL[i];
    int ix = APVX_asymALL[i];

    bool useasym = false;
    
    double Gx=1.0, Gy=1.0;
    
    if( iy >= 0 && iy < nAPVmaxY && ix >= 0 && ix < nAPVmaxX ){
      Gx = par[ix + nAPVmaxY];
      Gy = par[iy];

      useasym = true;
    } else if( iy >= 0 && iy < nAPVmaxY ){ //using weighted average over all X APVs for Gx:
      double sum_Gx = 0.0;
      double sum_Wx = 0.0;
      for( int j = 0; j<nAPVmaxX; j++ ){
	sum_Gx += weightX[iy][j] * par[j+nAPVmaxY];
	sum_Wx += weightX[iy][j];
      }

      Gx = sum_Gx/sum_Wx;
      Gy = par[iy];

      if( countX[iy] == 0 ) useasym = true;
      
    } else if( ix >= 0 && ix < nAPVmaxX ){ //using weighted average over all Y APVs for Gy:
      double sum_Gy = 0.0;
      double sum_Wy = 0.0;
      for( int j=0; j<nAPVmaxY; j++ ){
	sum_Gy += weightY[ix][j] * par[j];
	sum_Wy += weightY[ix][j];
      }

      Gx = par[ix+nAPVmaxY];
      Gy = sum_Gy/sum_Wy;

      if( countY[ix] == 0 ) useasym = true;
	
    }
      
    double Atheory = (Gx - Gy)/(Gx + Gy);

    if( useasym ) chi2 += pow( (AsymALL[i] - Atheory)/dAsymALL[i], 2 );
    
  }

  //cout << "chi2 = " << chi2 << endl;
  
  f = chi2;
}


void GainRatios( const char *infilename, int nmodules, double chi2cut=100.0, double ADCcut = 1500.0, int nstripxmax=1536, int nstripymax=1280 ){
  gROOT->ProcessLine(".x ~/rootlogon.C");

  TString fname(infilename);

  TChain *C = new TChain("Tout");

  C->Add( infilename );

  TString outfilename = fname;

  outfilename.Prepend("GainRatios_");

  TFile *fout = new TFile( outfilename.Data(), "RECREATE" );
  
  outfilename.ReplaceAll(".root",".txt");

  ofstream outfile(outfilename.Data());

  nAPVmaxX = nstripxmax/128;
  if( nstripxmax % 128 > 0 ) nAPVmaxX++;

  nAPVmaxY = nstripymax/128;
  if( nstripymax % 128 > 0 ) nAPVmaxY++;
  
  TClonesArray *hADCasym_vs_APVXY = new TClonesArray( "TH1D", nAPVmaxX*nAPVmaxY*nmodules );
  TClonesArray *hADCasym_vs_APVX = new TClonesArray( "TH1D", nAPVmaxX*nmodules );
  TClonesArray *hADCasym_vs_APVY = new TClonesArray( "TH1D", nAPVmaxY*nmodules );
  
  for( int imodule=0; imodule<nmodules; imodule++ ){
    TString hname;
    for( int iAPVx = 0; iAPVx < nAPVmaxX; iAPVx++ ){
      hname.Form( "hADCasym_vs_APVX%d_mod%d", iAPVx, imodule );
      new( (*hADCasym_vs_APVX)[iAPVx+nAPVmaxX*imodule] ) TH1D( hname.Data(), "", 250, -1.01, 1.01 );
      
      
      for( int iAPVy = 0; iAPVy < nAPVmaxY; iAPVy++ ){
	if( iAPVx==0 ){
	  hname.Form( "hADCasym_vs_APVY%d_mod%d", iAPVy, imodule );
	  new( (*hADCasym_vs_APVY)[iAPVy+nAPVmaxY*imodule] ) TH1D( hname.Data(), "", 250, -1.01, 1.01 );
	}
	
	hname.Form("hADCasym_vs_APV_mod%d_x%d_y%d",imodule,iAPVx,iAPVy);
	
	new( (*hADCasym_vs_APVXY)[iAPVy+nAPVmaxY*iAPVx+nAPVmaxX*nAPVmaxY*imodule] ) TH1D( hname.Data(), "", 250,-1.01,1.01);
      }
    }
  }
  
  GEM_cosmic_tracks *T = new GEM_cosmic_tracks(C);

  long nevent=0;

  TH2D *hADCasym_module = new TH2D("hADCasym_module","",nmodules,-0.5,nmodules-0.5,500,-1.01,1.01);
  TH2D *hNstripX_module = new TH2D("hNstripX_module","",nmodules,-0.5,nmodules-0.5,12,0.5,12.5);
  TH2D *hNstripY_module = new TH2D("hNstripY_module","",nmodules,-0.5,nmodules-0.5,12,0.5,12.5);

  //int nAPVmax = 
  
  while( C->GetEntry( nevent++ ) ){
    //loop over hits:
    if( T->TrackChi2NDF < chi2cut ){
      for( int ihit=0; ihit<T->TrackNhits; ihit++ ){
	if( sqrt(T->HitADCX[ihit]*T->HitADCY[ihit]) >= ADCcut ){
	  hADCasym_module->Fill( T->HitModule[ihit], T->HitADCasym[ihit] );
	  hNstripX_module->Fill( T->HitModule[ihit], T->HitNstripX[ihit] );
	  hNstripY_module->Fill( T->HitModule[ihit], T->HitNstripY[ihit] );

	  int ixlo = T->HitXstripLo[ihit];
	  int ixhi = T->HitXstripHi[ihit];
	  int ixmax = T->HitXstripMax[ihit];

	  int iylo = T->HitYstripLo[ihit];
	  int iyhi = T->HitYstripHi[ihit];
	  int iymax = T->HitYstripMax[ihit];

	  int xAPVmax = ixmax/128;
	  int yAPVmax = iymax/128;
	  int xAPVlo = ixlo/128;
	  int yAPVlo = iylo/128;

	  int xAPVhi = ixhi/128;
	  int yAPVhi = iyhi/128;

	  if( xAPVlo == xAPVmax && xAPVhi == xAPVmax &&
	      yAPVlo == yAPVmax && yAPVhi == yAPVmax &&
	      T->HitNstripX[ihit] >= 2 && T->HitNstripY[ihit] >= 2 ){

	    ( (TH1D*) (*hADCasym_vs_APVXY)[yAPVmax + nAPVmaxY*xAPVmax+nAPVmaxX*nAPVmaxY*T->HitModule[ihit]] )->Fill( T->HitADCasym[ihit] );

	    ( (TH1D*) (*hADCasym_vs_APVX)[xAPVmax + nAPVmaxX*T->HitModule[ihit]] )->Fill( T->HitADCasym[ihit] );
	    ( (TH1D*) (*hADCasym_vs_APVY)[yAPVmax + nAPVmaxY*T->HitModule[ihit]] )->Fill( T->HitADCasym[ihit] );
	    
	  }
	}
	////
      }
    }
  }

  // 

  double asympeak[nmodules];
  double R[nmodules];
  
  TH1D *htemp;

  //outfile << "mod_RYX     ";

  TCanvas *c1 = new TCanvas("c1","c1",2000,1000);
  c1->Divide(2,1);
  c1->cd(1);

  hADCasym_module->Draw("colz");

  gPad->Modified();
  c1->Update();
  
  c1->cd(2);

  
  
  for( int i=0; i<nmodules; i++ ){
    TString hnametemp;
    hnametemp.Form("ADCasym_module%d",i);
    htemp = hADCasym_module->ProjectionY(hnametemp.Data(), i+1, i+1 );

    double max = htemp->GetMaximum();
    int binmax,binlo,binhi;

    binmax = htemp->GetMaximumBin();
    binlo = binmax;
    binhi = binmax;

    while( htemp->GetBinContent(binlo) >= 0.5*max && binlo > 1 ){binlo--;}
    while( htemp->GetBinContent(binhi) >= 0.5*max && binhi < 500 ){binhi++; }

    if( htemp->GetEntries() >= 100 ){
    
      htemp->Fit("gaus","S","",htemp->GetBinCenter(binlo),htemp->GetBinCenter(binhi) );
      
      TF1 *fitfunc = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");
      
      asympeak[i] = fitfunc->GetParameter(1);
      R[i] = (1.0-asympeak[i])/(1.0+asympeak[i]);
      
      gPad->Modified();
      c1->Update();
      
      //      gSystem->Sleep(250);
      TString outstring;
    //outfile << outstring.Format(" %12.6g ", R[i] );

    //Within each module, we measure an asymmetry for all possible combinations of X APV (or "U") and Y APV (or "V"):

    //Let us find the gain coefficients for each individual APV card that minimize the chi^2 defined as the
    // sum of squared differences between the modified ADC asymmetries and zero:
    // If APV card X has relative gain Gx and card Y has relative gain Gy, then the observed ADC asymmetry would be:
    // ASYM = ( Gx - Gy )/(Gx + Gy ) = (1 - Ryx_i)/(1 + Ryx_i)
    // Let's choose some reference APVX as having gain = 1. So all gains are measured relative to this common
    // reference. 
    // Ryx = (1-A)/(1+A),
    // For sufficiently small A, Ryx - 1 = (1-A)/(1+A) - 1 = -2A/(1+A) ~= -2A
    //We have asymmetries A_ij 
    }
  }

  for( int i=0; i<nmodules; i++ ){
    int xAPV_ref = nAPVmaxX/2;
      
      //First, determine all the Y gains relative to the reference APV,
    // A_{ij} = (1-R_{ij})/(1+R_{ij})
    // Ryx_{ij] = (1-A_{ij})/(1+A_{ij}) = Gyi/Gxj
    // 
    //chi^2 = sum_i,j (A_{ij}^(measured) - A_{ij}^{calc})^2/(Delta A_{ij}^{measured})^2
    // Is there a way to make this problem linear?
    // Maybe:
    //For sufficiently small values of R_{yx} - 1, we have:
    
    //A_{ij}^{calc} = (1 - Ryx_{ij})/(1+ Ryx_{ij}) 

    vector<double> Asym, dAsym;
    vector<double> Ryx, dRyx;

    vector<double> Xgain, dXgain, Ygain, dYgain; //by APV within this module:

    double minAsymerr = 1000.0;

    int count=0;

    AsymALL.clear();
    dAsymALL.clear();

    TMinuit *gainfit = new TMinuit( nAPVmaxX+nAPVmaxY );

    TString fitname;
    gainfit->SetName( fitname.Format("gainfit_module%d",i) ); 
    
    gainfit->SetFCN(chi2_FCN);
    
    double arglist[10];
    
    int ierflg=0;
    
    for( int ipar=0; ipar<nAPVmaxY; ipar++ ){
      TString parname;
      parname.Form("Gain_APVY%d",ipar);
      
      gainfit->mnparm( ipar, parname.Data(), 1.0, 0.1,0,0,ierflg ); 
      
    }
    
    for( int ipar=0; ipar<nAPVmaxX; ipar++ ){
      TString parname;
      parname.Form("Gain_APVX%d",ipar);
      
      gainfit->mnparm( ipar+nAPVmaxY, parname.Data(), 1.0, 0.1,0,0,ierflg ); 
      
    }
    AsymALL.clear();
    dAsymALL.clear();

    APVX_asymALL.clear();
    APVY_asymALL.clear();
    
    weightX.clear();
    weightY.clear();

    weightX.resize(nAPVmaxY);
    weightY.resize(nAPVmaxX);
    //Shall we just create a separate TMinuit object for each module? Why TF not?

    countX.resize(nAPVmaxY);
    countY.resize(nAPVmaxX); 
    
    //Card-average asymmetries:
    for( int ix=0; ix<nAPVmaxX; ix++ ){
      weightY[ix].resize(nAPVmaxY);

      countY[ix] = 0;
      
      TH1D *htemp = ( (TH1D*) (*hADCasym_vs_APVX)[ix+nAPVmaxX*i] );
      if( htemp->GetEntries() >= 100 ){
	double Amean = htemp->GetMean();
	double dAmean = htemp->GetMeanError();

	double Arms = htemp->GetRMS();
	double dArms = htemp->GetRMSError();

	htemp->Fit("gaus","SQ","",Amean-Arms, Amean+Arms);
	 
	double Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	double dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);

	double Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	double dArms_fit = htemp->GetFunction("gaus")->GetParError(2);

	//re-run the fit with a tighter range set by the result of the first fit:

	htemp->Fit( "gaus", "SQ", "", Amean_fit - 2.0*Arms_fit, Amean_fit + 2.0*Arms_fit );

	Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);
	  
	Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	dArms_fit = htemp->GetFunction("gaus")->GetParError(2);

	double Ryxtemp = (1.0-Amean_fit)/(1.0+Amean_fit);

	double Aup = Amean_fit + dAmean_fit;
	double Adown = Amean_fit - dAmean_fit;
	  
	double Ryx_Aup = (1.0 - Aup)/(1.0 + Aup );
	double Ryx_Adown = (1.0 - Adown)/(1.0 + Adown);

	AsymALL.push_back( Amean_fit );
	dAsymALL.push_back( dAmean_fit );

	APVX_asymALL.push_back( ix );
	APVY_asymALL.push_back( -1 );
      } else { //If less than 100 entries for this entire APV, fix the gain to 1:
	gainfit->FixParameter( ix + nAPVmaxY );
      }
    }

    for( int iy=0; iy<nAPVmaxY; iy++ ){
      weightX[iy].resize( nAPVmaxX );

      countX[iy] = 0;
      
      TH1D *htemp = ( (TH1D*) (*hADCasym_vs_APVY)[iy+nAPVmaxY*i] );
      if( htemp->GetEntries() >= 100 ){
	double Amean = htemp->GetMean();
	double dAmean = htemp->GetMeanError();

	double Arms = htemp->GetRMS();
	double dArms = htemp->GetRMSError();

	htemp->Fit("gaus","SQ","",Amean-Arms, Amean+Arms);
	 
	double Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	double dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);

	double Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	double dArms_fit = htemp->GetFunction("gaus")->GetParError(2);

	//re-run the fit with a tighter range set by the result of the first fit:

	htemp->Fit( "gaus", "SQ", "", Amean_fit - 2.0*Arms_fit, Amean_fit + 2.0*Arms_fit );

	Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);
	  
	Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	dArms_fit = htemp->GetFunction("gaus")->GetParError(2);

	double Ryxtemp = (1.0-Amean_fit)/(1.0+Amean_fit);

	double Aup = Amean_fit + dAmean_fit;
	double Adown = Amean_fit - dAmean_fit;
	  
	double Ryx_Aup = (1.0 - Aup)/(1.0 + Aup );
	double Ryx_Adown = (1.0 - Adown)/(1.0 + Adown);

	AsymALL.push_back( Amean_fit );
	dAsymALL.push_back( dAmean_fit );

	APVX_asymALL.push_back( -1 );
	APVY_asymALL.push_back( iy );
      } else {
	gainfit->FixParameter( iy );
      }
    }
    
    for( int ix = 0; ix<nAPVmaxX; ix++ ){
      
      for( int iy = 0; iy<nAPVmaxY; iy++ ){
	TH1D *htemp = ( (TH1D*) (*hADCasym_vs_APVXY)[iy + nAPVmaxY*ix + nAPVmaxY*nAPVmaxX*i] );

	
	
	weightX[iy][ix] = htemp->GetEntries();
	weightY[ix][iy] = htemp->GetEntries();
	
	if( htemp->GetEntries() >= 100 ){

	  countY[ix]++;
	  countX[iy]++;
	  
	  double Amean = htemp->GetMean();
	  double dAmean = htemp->GetMeanError();

	  double Arms = htemp->GetRMS();
	  double dArms = htemp->GetRMSError();

	  htemp->Fit("gaus","SQ","",Amean-Arms, Amean+Arms);
	 
	  double Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	  double dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);

	  double Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	  double dArms_fit = htemp->GetFunction("gaus")->GetParError(2);

	  //re-run the fit with a tighter range set by the result of the first fit:

	  htemp->Fit( "gaus", "SQ", "", Amean_fit - 2.0*Arms_fit, Amean_fit + 2.0*Arms_fit );

	  Amean_fit = htemp->GetFunction("gaus")->GetParameter("Mean");
	  dAmean_fit = htemp->GetFunction("gaus")->GetParError(1);
	  
	  Arms_fit = htemp->GetFunction("gaus")->GetParameter("Sigma");
	  dArms_fit = htemp->GetFunction("gaus")->GetParError(2);

	  double Ryxtemp = (1.0-Amean_fit)/(1.0+Amean_fit);

	  double Aup = Amean_fit + dAmean_fit;
	  double Adown = Amean_fit - dAmean_fit;
	  
	  double Ryx_Aup = (1.0 - Aup)/(1.0 + Aup );
	  double Ryx_Adown = (1.0 - Adown)/(1.0 + Adown);

	  
	  
	  Asym.push_back( Amean_fit );
	  dAsym.push_back( dAmean_fit );
	  Ryx.push_back( Ryxtemp );
	  dRyx.push_back( 0.5*fabs(Ryx_Aup - Ryx_Adown) );

	  //ONLY populate the list of asymmetries IF we had a successful fit:
	  AsymALL.push_back( Amean_fit );
	  dAsymALL.push_back( dAmean_fit );

	  APVX_asymALL.push_back( ix );
	  APVY_asymALL.push_back( iy );
	  
	  if( count == 0 || dAmean_fit < minAsymerr ){
	    minAsymerr = dAmean_fit;
	    xAPV_ref = ix;
	    count++;
	  }
	} else {
	  Asym.push_back( -10.0 );
	  dAsym.push_back( 1.0 );
	  Ryx.push_back( 1.0 );
	  dRyx.push_back( 1.0 );
	}
      }
    }

    arglist[0] = 5000;
    arglist[1] = 0.1;
    gainfit->mnexcm("MIGRAD",arglist,2,ierflg);
    
    //First calculate all the Y gains relative to the reference x APV

    outfile << "mod_Ygain  " << i << "   " << nAPVmaxY << "     ";
    
    for( int iy=0; iy<nAPVmaxY; iy++ ){
      int yidx= iy + nAPVmaxY * xAPV_ref;
      
      // Ygain.push_back( Ryx[yidx] );
      // dYgain.push_back( dRyx[yidx] );

      double Gy, dGy;
      gainfit->GetParameter( iy, Gy, dGy );
      Ygain.push_back( Gy );
      dYgain.push_back( dGy );
      
      cout << "module " << i << ", Y APV " << iy << ", Relative gain = " << Ygain.back() << " +/- " << dYgain.back() << endl;

      outfile << Ygain.back() << "  ";
    }
    outfile << endl;

    outfile << "mod_Xgain  " << i << "   " << nAPVmaxX << "     ";
    
    for( int ix=0; ix<nAPVmaxX; ix++ ){

      double Gx, dGx;
      gainfit->GetParameter( ix + nAPVmaxY, Gx, dGx );
      Xgain.push_back( Gx );
      dXgain.push_back( dGx );
      // if( ix != xAPV_ref ){
      // 	double sum_xgain = 0.0;
      // 	double sum_weights = 0.0;
      // 	//to get relative X gains, compute a weighted average of the relative X gains from
      // 	//the relative Y gains: 
      // 	for( int iy=0; iy<nAPVmaxY; iy++ ){
      // 	  int idx = iy + nAPVmaxY*ix;
      // 	  if( fabs( Asym[idx] ) <= 1.0 ){ //then we have successfully fit the asymmetry for this combination:
      // 	    double ygaintemp = Ygain[iy];

      // 	    double dygaintemp = dYgain[iy];

      // 	    double xgaintemp = Ygain[iy]/Ryx[idx];
      // 	    double dxgaintemp = fabs(xgaintemp)*sqrt(pow( dygaintemp/ygaintemp, 2 ) + pow( dRyx[idx]/Ryx[idx],2));

      // 	    sum_xgain += xgaintemp * pow( dxgaintemp, -2 );
      // 	    sum_weights += pow( dxgaintemp, -2 );
      // 	  }
      // 	}

      // 	if( sum_xgain > 0.0 ){
      // 	  Xgain.push_back( sum_xgain/sum_weights );
      // 	  dXgain.push_back( 1.0/sqrt(sum_weights) );
      // 	} else {
      // 	  Xgain.push_back( 1.0 );
      // 	  dXgain.push_back( 1.0 );
      // 	}
      // } else {
      // 	Xgain.push_back( 1.0 );
      // 	dXgain.push_back( 0.0 );
      // }
      cout << "module " << i << ", X APV " << ix << ", Relative gain = " << Xgain.back() << " +/- " << dXgain.back() << endl;

      outfile << Xgain.back() << "  ";
      
    }
    outfile << endl;
    
  }

  outfile << endl;

  fout->Write();

}
