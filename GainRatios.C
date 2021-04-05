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

  int nAPVmaxX = nstripxmax/128;
  if( nstripxmax % 128 > 0 ) nAPVmaxX++;

  int nAPVmaxY = nstripymax/128;
  if( nstripymax % 128 > 0 ) nAPVmaxY++;
  
  TClonesArray *hADCasym_vs_APVXY = new TClonesArray( "TH1D", nAPVmaxX*nAPVmaxY*nmodules );

  for( int imodule=0; imodule<nmodules; imodule++ ){
    for( int iAPVx = 0; iAPVx < nAPVmaxX; iAPVx++ ){
      for( int iAPVy = 0; iAPVy < nAPVmaxY; iAPVy++ ){
	
	TString hname;
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
	}

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
	  
	}
	
      }
    }
  }

  

  double asympeak[nmodules];
  double R[nmodules];
  
  TH1D *htemp;

  outfile << "mod_RYX     ";

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

    htemp->Fit("gaus","S","",htemp->GetBinCenter(binlo),htemp->GetBinCenter(binhi) );

    TF1 *fitfunc = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");
    
    asympeak[i] = fitfunc->GetParameter(1);
    R[i] = (1.0-asympeak[i])/(1.0+asympeak[i]);

    gPad->Modified();
    c1->Update();

    gSystem->Sleep(250);
    TString outstring;
    outfile << outstring.Format(" %12.6g ", R[i] );
  }

  outfile << endl;

  fout->Write();

}
