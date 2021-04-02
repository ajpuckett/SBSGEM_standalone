#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;

void GEM_GainMatch( const char *rootfilename, int nmodules ){

  gStyle->SetOptFit();
  gStyle->SetNdivisions(505,"XYZ");
  
  TFile *f = new TFile(rootfilename,"READ");

  TString outfilename = rootfilename;

  outfilename.ReplaceAll(".root",".txt");
  outfilename.Prepend("GainMatching_");

  ofstream outfile(outfilename.Data());
  

  TString pdffilename = outfilename;
  pdffilename.ReplaceAll(".txt",".pdf");
  
  TH1D *hclustersumX_all, *hclustersumY_all;

  f->GetObject("hADCsumXclust_max",hclustersumX_all);
  f->GetObject("hADCsumYclust_max",hclustersumY_all);

  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);

  c1->Divide(2,1);

  c1->cd(1);
  hclustersumX_all->Draw();
  hclustersumX_all->GetXaxis()->SetRangeUser(0,20000);
  hclustersumX_all->Fit("landau","","",1000,20000);
  hclustersumX_all->GetFunction("landau")->SetNpx(1000);

  c1->cd(2);
  hclustersumY_all->Draw();
  hclustersumY_all->GetXaxis()->SetRangeUser(0,20000);
  hclustersumY_all->Fit("landau","","",1000,20000);
  hclustersumY_all->GetFunction("landau")->SetNpx(1000);

  TString printname = pdffilename + "(";
  c1->Print(printname.Data(), "pdf");
  
  double peakposX_all = hclustersumX_all->GetFunction("landau")->GetParameter("MPV");
  double peakposY_all = hclustersumY_all->GetFunction("landau")->GetParameter("MPV");

  cout << "Cluster sum X ADC peak position, all modules = " << peakposX_all << endl;
  cout << "Cluster sum Y ADC peak position, all modules = " << peakposY_all << endl;

  int maxAPVs_per_module = 27;

  
  //c1->Divide(4,3);

  TH2D *hADCX_vs_module,*hADCY_vs_module;
  f->GetObject( "hADCsumXclust_module", hADCX_vs_module );
  f->GetObject( "hADCsumYclust_module", hADCY_vs_module );

  bool firstplot = true;

  map<int,vector<double> > mod_relative_gainX_vsAPV;
  map<int,vector<double> > mod_relative_gainY_vsAPV;
  
  for( int imodule=0; imodule<nmodules; imodule++ ){
    
    //if( imodule == 0 ) pdffilename += "(";
    //Do X strips first:
    TH2D *htemp;

    TString hname;

    hname.Form( "hADCX_mod%d" , imodule );

    //Plot module averages:
    c1->Clear();
    c1->Divide(1,2);
    
    TH1D *hADCX_modavg =  hADCX_vs_module->ProjectionX(hname.Data(), imodule+1, imodule+1 );
    
    c1->cd(1);
    hADCX_modavg->GetXaxis()->SetRangeUser(0,20000);
    hADCX_modavg->DrawCopy();
    hADCX_modavg->Fit("landau","","",1000,20000);
    hADCX_modavg->GetFunction("landau")->SetNpx(1000);
    
    hname.Form( "hADCY_mod%d", imodule );
    TH1D *hADCY_modavg = hADCY_vs_module->ProjectionX(hname.Data(), imodule+1, imodule+1 );
    c1->cd(2);
    hADCY_modavg->GetXaxis()->SetRangeUser(0,20000);
    hADCY_modavg->DrawCopy();
    hADCY_modavg->Fit("landau","","",1000,20000);
    hADCY_modavg->GetFunction("landau")->SetNpx(1000);
    
    double Xpeak_modavg = hADCX_modavg->GetFunction("landau")->GetParameter("MPV");
    double Ypeak_modavg = hADCY_modavg->GetFunction("landau")->GetParameter("MPV");

    double TargetADC = 0.5*(Xpeak_modavg+Ypeak_modavg);

    double RelativeGain_Xavg = TargetADC/Xpeak_modavg;
    double RelativeGain_Yavg = TargetADC/Ypeak_modavg;
    
    // if( imodule == 0 ){
    //   TString printname = pdffilename + "(";
    //   c1->Print( printname.Data(), "pdf" );
    // } else {
    c1->Print(pdffilename.Data(), "pdf" );
      //}

    c1->Clear();
    
    
    hname.Form( "hmod_clustersumX_vs_APV_module%d", imodule );

    f->GetObject( hname.Data(), htemp );

    int nAPVs = std::min( maxAPVs_per_module, htemp->GetNbinsX() );

    //The formula for canvas division should be 
    int canv_nrows = int(std::round( sqrt(double(nAPVs)) ) );

    int canv_ncols = nAPVs/canv_nrows;
    if ( nAPVs%canv_nrows != 0 ) canv_ncols++;

    c1->Divide( canv_ncols, canv_nrows );

    for( int bin=1; bin<=nAPVs; bin++ ){
      hname.Form( "hADCX_mod%d_APV%d", imodule, bin );
      TH1D *hADCtemp = htemp->ProjectionY( hname.Data(), bin, bin );

      c1->cd( bin );
      hADCtemp->GetXaxis()->SetRangeUser(0,20000);
      hADCtemp->DrawCopy();

      double APV_peakpos = 0.0;

      if( hADCtemp->GetEntries() >= 1000 ){
      
	hADCtemp->Fit("landau","","",1000,20000);
	hADCtemp->GetFunction("landau")->SetNpx(1000);
	
	APV_peakpos = hADCtemp->GetFunction("landau")->GetParameter("MPV");
      }
	
      if( hADCtemp->GetEntries() < 1000 || APV_peakpos < 1000.0 ){
	mod_relative_gainX_vsAPV[imodule].push_back( RelativeGain_Xavg );
      } else {
	mod_relative_gainX_vsAPV[imodule].push_back( TargetADC/APV_peakpos );
      }
      
    }

    c1->Print( pdffilename.Data(), "pdf");

    hname.Form( "hmod_clustersumY_vs_APV_module%d", imodule );

    f->GetObject( hname.Data(), htemp );

    nAPVs = std::min( maxAPVs_per_module, htemp->GetNbinsX() );

    //The formula for canvas division should be 
    canv_nrows = int(std::round( sqrt(double(nAPVs)) ) );

    canv_ncols = nAPVs/canv_nrows;
    if ( nAPVs%canv_nrows != 0 ) canv_ncols++;
    
    //Now do Y's:
    c1->Clear();
    c1->Divide( canv_ncols, canv_nrows );
    
    for( int bin=1; bin<=nAPVs; bin++ ){
      hname.Form( "hADCY_mod%d_APV%d", imodule, bin );
      TH1D *hADCtemp = htemp->ProjectionY( hname.Data(), bin, bin );

      c1->cd( bin );
      hADCtemp->GetXaxis()->SetRangeUser(0,20000);
      hADCtemp->DrawCopy();

      double APV_peakpos = 0.0;
      
      if( hADCtemp->GetEntries() >= 1000 ){
      
	hADCtemp->Fit("landau","","",1000,20000);
	
	hADCtemp->GetFunction("landau")->SetNpx(1000);
	
	APV_peakpos = hADCtemp->GetFunction("landau")->GetParameter("MPV");
      } 


      if( hADCtemp->GetEntries() < 1000 || APV_peakpos < 1000.0 ){
	mod_relative_gainY_vsAPV[imodule].push_back( RelativeGain_Yavg );
      } else{
	mod_relative_gainY_vsAPV[imodule].push_back( TargetADC/APV_peakpos );
      }      
    }

    if( imodule+1 == nmodules ){
      TString printname = pdffilename +")";
      c1->Print(printname.Data(), "pdf" );
    } else {
      c1->Print( pdffilename.Data(), "pdf" );
    }


    //Now write out
    outfile << "mod_Xgain  " << imodule << "   " << mod_relative_gainX_vsAPV[imodule].size() << "     ";
    for( int iAPV = 0; iAPV<mod_relative_gainX_vsAPV[imodule].size(); iAPV++ ){
      outfile << mod_relative_gainX_vsAPV[imodule][iAPV] << "  ";
    }
    outfile << endl;

    outfile << "mod_Ygain  " << imodule << "   " << mod_relative_gainY_vsAPV[imodule].size() << "     ";
    for( int iAPV = 0; iAPV<mod_relative_gainY_vsAPV[imodule].size(); iAPV++ ){
      outfile << mod_relative_gainY_vsAPV[imodule][iAPV] << "  ";
    }
    outfile << endl;
    
  }
  


}
