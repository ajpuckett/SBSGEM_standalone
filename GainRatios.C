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

void GainRatios( const char *infilename, int nmodules, double chi2cut=100.0, double ADCcut = 1500.0 ){
  gROOT->ProcessLine(".x ~/rootlogon.C");

  TString fname(infilename);

  TChain *C = new TChain("Tout");

  C->Add( infilename );

  TString outfilename = fname;

  outfilename.Prepend("GainRatios_");
  outfilename.ReplaceAll(".root",".txt");

  ofstream outfile(outfilename.Data());

  GEM_cosmic_tracks *T = new GEM_cosmic_tracks(C);

  long nevent=0;

  TH2D *hADCasym_module = new TH2D("hADCasym_module","",nmodules,-0.5,nmodules-0.5,500,-1.01,1.01);
  TH2D *hNstripX_module = new TH2D("hNstripX_module","",nmodules,-0.5,nmodules-0.5,12,0.5,12.5);
  TH2D *hNstripY_module = new TH2D("hNstripY_module","",nmodules,-0.5,nmodules-0.5,12,0.5,12.5);
  
  while( C->GetEntry( nevent++ ) ){
    //loop over hits:
    if( T->TrackChi2NDF < chi2cut ){
      for( int ihit=0; ihit<T->TrackNhits; ihit++ ){
	if( sqrt(T->HitADCX[ihit]*T->HitADCY[ihit]) >= ADCcut ){
	  hADCasym_module->Fill( T->HitModule[ihit], T->HitADCasym[ihit] );
	  hNstripX_module->Fill( T->HitModule[ihit], T->HitNstripX[ihit] );
	  hNstripY_module->Fill( T->HitModule[ihit], T->HitNstripY[ihit] );
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

  

}
