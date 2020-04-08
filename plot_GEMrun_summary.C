void plot_GEMrun_summary(const char *filename, int nlayers=4, int nmodules=12, double chi2cut=50.0, double CALOCUT=0.0){
  gROOT->ProcessLine(".x ~/rootlogon.C");

  gStyle->SetOptStat(0);

  gStyle->SetTitleSize( .05, "XYZ");
  gStyle->SetLabelSize( .05, "XYZ");
  gStyle->SetNdivisions(505,"XYZ");

  gStyle->SetStatW(0.25);
  //gStyle->SetStatH(0.45);
  
  TGaxis::SetMaxDigits(3);

  gROOT->ForceStyle();
  
  TFile *fin = new TFile(filename,"READ");

  TCanvas *c1 = new TCanvas("c1","c1",1400,1050);

  c1->Divide(2,2,.001,.001);

  double lmargin=0.18,rmargin=0.12,tmargin=0.09,bmargin=0.12;

  TH2D *hNstripsX_layer,*hNstripsY_layer, *hNclustX_layer;
  TH1D *hNlayers_2Dclust;
  
  fin->GetObject("hNstripsX_layer",hNstripsX_layer);
  fin->GetObject("hNstripsY_layer",hNstripsY_layer);
  fin->GetObject("hNclustX_layer",hNclustX_layer);
  fin->GetObject("hNlayers_2Dclust",hNlayers_2Dclust);
  
  hNstripsX_layer->SetTitle("");
  hNstripsX_layer->GetXaxis()->SetRangeUser(-0.5,20.5);
  hNstripsX_layer->GetXaxis()->SetTitle("Number of X strips fired per layer");
  hNstripsX_layer->GetYaxis()->SetTitle("Layer");

  hNstripsY_layer->SetTitle("");
  hNstripsY_layer->GetXaxis()->SetRangeUser(-0.5,20.5);
  hNstripsY_layer->GetXaxis()->SetTitle("Number of Y strips fired per layer");
  hNstripsY_layer->GetYaxis()->SetTitle("Layer");

  hNclustX_layer->SetTitle("");
  hNclustX_layer->GetXaxis()->SetRangeUser(-0.5,20.5);
  hNclustX_layer->GetXaxis()->SetTitle("Number of reconstructed 2D clusters per layer");
  hNclustX_layer->GetYaxis()->SetTitle("Layer");

  hNlayers_2Dclust->SetTitle("");
  hNlayers_2Dclust->GetXaxis()->SetTitle("Number of layers with 2D cluster reconstructed");
  
		 
		 
  lmargin=0.12;
  rmargin=0.15;
  bmargin=0.12;
  tmargin=0.09;
  
  c1->cd(1);
  gPad->SetLeftMargin(lmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);

  hNstripsX_layer->Draw("colz");
  
  c1->cd(2);
  gPad->SetLeftMargin(lmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);

  hNstripsY_layer->Draw("colz");

  c1->cd(3);
  gPad->SetLeftMargin(lmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);

  hNclustX_layer->Draw("colz");

  c1->cd(4);
  
  gPad->SetLeftMargin(lmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);

  hNlayers_2Dclust->Draw();

  TString pdffilename = filename;
  pdffilename.ReplaceAll( ".root", ".pdf");

  pdffilename.Prepend("plots/");
  
  TString openfilename = pdffilename+"[";
  TString closefilename = pdffilename+"]";
  
  c1->Print(openfilename);
  c1->Print(pdffilename);

  TTree *Tout;
  fin->GetObject("Tout",Tout);

  TH1D *hclustwidthX = new TH1D("hclustwidthX","",12,-0.5,11.5);
  TH1D *hclustwidthY = new TH1D("hclustwidthY","",12,-0.5,11.5);
  TH1D *hNclustPerLayer = new TH1D("hNclustPerLayer","",21,-0.5,20.5);
  TH2D *hClust2D_NstripX_vs_NstripY;
  fin->GetObject("hClust2D_NstripX_vs_NstripY",hClust2D_NstripX_vs_NstripY);
  
  Tout->Project("hclustwidthX","HitNstripX");
  Tout->Project("hclustwidthY","HitNstripY");
  Tout->Project("hNclustPerLayer","Ncluster[HitLayer]");

  hclustwidthX->GetXaxis()->SetTitle("Width of cluster in X strips");
  hclustwidthY->GetXaxis()->SetTitle("Width of cluster in Y strips");
  hNclustPerLayer->GetXaxis()->SetTitle("Number of clusters per layer (layer on track)");
  hClust2D_NstripX_vs_NstripY->GetXaxis()->SetTitle("Cluster width in Y strips");
  hClust2D_NstripX_vs_NstripY->GetYaxis()->SetTitle("Cluster width in X strips");
  
  c1->cd(1);
  hclustwidthX->Draw();
  c1->cd(2);
  hclustwidthY->Draw();
  c1->cd(3);
  hNclustPerLayer->Draw();
  c1->cd(4);
  hClust2D_NstripX_vs_NstripY->Draw("colz");

  c1->Print(pdffilename);

  TH2D *hClust2D_ADCXvsY,*hClust2D_T0XvsY;
  TH1D *hClust2D_ADCasym,*hClust2D_Tdiff;

  fin->GetObject( "hClust2D_ADCXvsY", hClust2D_ADCXvsY );
  fin->GetObject( "hClust2D_T0XvsY", hClust2D_T0XvsY );
  fin->GetObject( "hClust2D_ADCasym", hClust2D_ADCasym );
  fin->GetObject( "hClust2D_Tdiff", hClust2D_Tdiff );

  TH1D *hClust2D_ADCdiff,*hNtracks_found, *hNhitspertrack, *hTrackChi2NDF, *hClust_corr, *hStrip_maxcor,*hTrackXp,*hTrackYp;
  TH2D *hClust2D_ADCasym_vs_ADCavg, *hTrackXresid_vs_layer, *hTrackYresid_vs_layer, *hTrackXresid_vs_module, *hTrackYresid_vs_module, *hTrackXY;
  TH2D *hTrackXeresid_vs_layer, *hTrackYeresid_vs_layer, *hTrackXeresid_vs_module, *hTrackYeresid_vs_module;

  
  fin->GetObject( "hClust2D_ADCdiff", hClust2D_ADCdiff );
  fin->GetObject( "hNtracks_found", hNtracks_found );
  fin->GetObject( "hNhitspertrack", hNhitspertrack );
  fin->GetObject( "hTrackChi2NDF", hTrackChi2NDF );
  fin->GetObject( "hClust_corr", hClust_corr );
  fin->GetObject( "hStrip_maxcor", hStrip_maxcor );
  fin->GetObject( "hTrackXY", hTrackXY );
  fin->GetObject( "hTrackXp", hTrackXp );
  fin->GetObject( "hTrackYp", hTrackYp );
  fin->GetObject( "hClust2D_ADCasym_vs_ADCavg", hClust2D_ADCasym_vs_ADCavg );
  fin->GetObject( "hTrackXresid_vs_layer", hTrackXresid_vs_layer );
  fin->GetObject( "hTrackYresid_vs_layer", hTrackYresid_vs_layer );
  fin->GetObject( "hTrackXresid_vs_module", hTrackXresid_vs_module );
  fin->GetObject( "hTrackYresid_vs_module", hTrackYresid_vs_module );
  fin->GetObject( "hTrackXeresid_vs_layer", hTrackXeresid_vs_layer );
  fin->GetObject( "hTrackYeresid_vs_layer", hTrackYeresid_vs_layer );
  fin->GetObject( "hTrackXeresid_vs_module", hTrackXeresid_vs_module );
  fin->GetObject( "hTrackYeresid_vs_module", hTrackYeresid_vs_module );

  TClonesArray *heff_layers = new TClonesArray( "TH2D", nlayers );
  TClonesArray *hxyhit_layers = new TClonesArray( "TH2D", nlayers );

  for( int ilayer=0; ilayer<nlayers; ilayer++ ){
    TH2D *htemp;
    TString histname;
    fin->GetObject( histname.Format("heff_layer%d",ilayer), htemp );

    new( (*heff_layers)[ilayer] ) TH2D( *htemp );

    fin->GetObject( histname.Format("hxyhit_layer%d",ilayer), htemp );
    
    new( (*hxyhit_layers)[ilayer] ) TH2D( *htemp );
  }
  
  // TH2D *heff_layer0,*heff_layer1,*heff_layer2,*heff_layer3;
  // TH2D *hxyhit_layer0,*hxyhit_layer1,*hxyhit_layer2,*hxyhit_layer3;

  // fin->GetObject( "heff_layer0", heff_layer0 );
  // fin->GetObject( "heff_layer1", heff_layer1 );
  // fin->GetObject( "heff_layer2", heff_layer2 );
  // fin->GetObject( "heff_layer3", heff_layer3 );
  // fin->GetObject( "hxyhit_layer0", hxyhit_layer0 );
  // fin->GetObject( "hxyhit_layer1", hxyhit_layer1 );
  // fin->GetObject( "hxyhit_layer2", hxyhit_layer2 );
  // fin->GetObject( "hxyhit_layer3", hxyhit_layer3 );
  
  lmargin = 0.15;
  
  c1->cd(1)->SetGrid();
  gPad->SetLeftMargin(lmargin);
  hClust2D_T0XvsY->SetTitle("");
  hClust2D_T0XvsY->GetXaxis()->SetTitle("Cluster t_{mean} (ns), X strips");
  hClust2D_T0XvsY->GetYaxis()->SetTitle("Cluster t_{mean} (ns), Y strips");
  
  hClust2D_T0XvsY->Draw("colz");

  c1->cd(2)->SetGrid();
  gPad->SetLeftMargin(lmargin);

  hClust2D_ADCXvsY->SetTitle("");
  hClust2D_ADCXvsY->GetXaxis()->SetTitle("Cluster ADC sum, X strips");
  hClust2D_ADCXvsY->GetYaxis()->SetTitle("Cluster ADC sum, Y strips");
  hClust2D_ADCXvsY->GetXaxis()->SetRangeUser(750.0,3e4);
  hClust2D_ADCXvsY->GetYaxis()->SetRangeUser(750.0,3e4);
  
  hClust2D_ADCXvsY->Draw("colz");

  c1->cd(3);
  gPad->SetLeftMargin(lmargin);
  hClust2D_Tdiff->Draw();
  hClust2D_Tdiff->SetTitle("");
  hClust2D_Tdiff->GetXaxis()->SetTitle( "#Delta t = t_{x} - t_{y} (ns)");
  hClust2D_Tdiff->GetXaxis()->SetRangeUser(-50,50);
  gStyle->SetOptFit();

  hClust2D_Tdiff->Fit("gaus","","",-6,6);
  
  c1->cd(4);
  gPad->SetLeftMargin(lmargin);
  hClust2D_ADCasym->SetTitle("");
  hClust2D_ADCasym->GetXaxis()->SetTitle("(ADCX-ADCY)/(ADCX+ADCY)");
  hClust2D_ADCasym->Draw();

  c1->Print(pdffilename);

  c1->cd(1);
  hClust2D_ADCasym_vs_ADCavg->SetTitle("");
  hClust2D_ADCasym_vs_ADCavg->GetYaxis()->SetRangeUser(-1,1);
  hClust2D_ADCasym_vs_ADCavg->GetXaxis()->SetRangeUser(0,4e4);
  hClust2D_ADCasym_vs_ADCavg->GetXaxis()->SetTitle("Cluster #sqrt{ADCX*ADCY}");
  hClust2D_ADCasym_vs_ADCavg->GetYaxis()->SetTitle("(ADCX-ADCY)/(ADCX+ADCY)");
  
  hClust2D_ADCasym_vs_ADCavg->Draw("col");

  c1->cd(2);
  hClust2D_ADCdiff->SetTitle("");
  hClust2D_ADCdiff->GetXaxis()->SetTitle("Cluster ADCX - ADCY");
  
  hClust2D_ADCdiff->Draw();

  c1->cd(3);
  hClust_corr->SetTitle("");
  hClust_corr->GetXaxis()->SetTitle("Cluster XY correlation coefficient");
  hClust_corr->Draw();

  c1->cd(4);
  hStrip_maxcor->SetTitle("");
  hStrip_maxcor->GetXaxis()->SetTitle("Best XY strip correlation coefficient in cluster");
  hStrip_maxcor->Draw();
  
  c1->Print(pdffilename);

  c1->Clear();
  c1->Divide(3,2,.001,.001);

  lmargin = 0.15;
  
  c1->cd(1);
  gPad->SetLeftMargin(lmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);
  
  hNtracks_found->SetTitle("");
  hNtracks_found->GetXaxis()->SetRangeUser(-0.5,5.5);

  double effavg = 1.0-hNtracks_found->GetBinContent(1)/hNtracks_found->GetEntries();

  TString htitle;
  hNtracks_found->SetTitle( htitle.Format( "Average track finding efficiency = %6.3g %%", effavg*100.0 ) );

  hNtracks_found->GetXaxis()->SetTitle("Number of tracks found per event");
  
  hNtracks_found->Draw();

  c1->cd(2);
  gPad->SetLeftMargin(lmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);

  hNhitspertrack->SetTitle("");
  hNhitspertrack->GetXaxis()->SetTitle("Number of hits on track");
  
  hNhitspertrack->Draw();

  c1->cd(3);
  gPad->SetLeftMargin(lmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);

  hTrackChi2NDF->SetTitle("");
  hTrackChi2NDF->GetXaxis()->SetRangeUser( 0.0, 150.0 );
  hTrackChi2NDF->GetXaxis()->SetTitle("Track #chi^{2}/dof");
  
  hTrackChi2NDF->Draw();

  c1->cd(4);
  gPad->SetLeftMargin(lmargin+.06);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);

  hTrackXY->GetYaxis()->SetTitle("Track X(Z=0) (mm)");
  hTrackXY->GetXaxis()->SetTitle("Track Y(Z=0) (mm)");
  
  hTrackXY->Draw("col");

  c1->cd(5);
  gPad->SetLeftMargin(lmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);
  hTrackXp->SetTitle("");
  hTrackXp->GetXaxis()->SetTitle("Track dx/dz fit");
  hTrackXp->Draw();

  c1->cd(6);
  gPad->SetLeftMargin(lmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);
  hTrackYp->SetTitle("");
  hTrackYp->GetXaxis()->SetTitle("Track dy/dz fit");
  //hTrackYp->GetXaxis()->SetRangeUser(-0.5,0.5);
  
  hTrackYp->Draw();
  
  c1->Print(pdffilename);

  c1->cd(1);
  hTrackXresid_vs_layer->SetTitle("");
  hTrackXresid_vs_layer->GetYaxis()->SetRangeUser(-3,3);
  hTrackXresid_vs_layer->GetYaxis()->SetTitle("Track X residual (mm)");
  hTrackXresid_vs_layer->GetXaxis()->SetTitle("layer");
  
  hTrackXresid_vs_layer->Draw("colz");

  c1->cd(4);

  hTrackYresid_vs_layer->SetTitle("");
  hTrackYresid_vs_layer->GetYaxis()->SetRangeUser(-3,3);
  hTrackYresid_vs_layer->GetYaxis()->SetTitle("Track Y residual (mm)");
  hTrackYresid_vs_layer->GetXaxis()->SetTitle("layer");
  hTrackYresid_vs_layer->Draw("colz");

  c1->cd(2);

  hTrackXresid_vs_module->SetTitle("");
  hTrackXresid_vs_module->GetYaxis()->SetRangeUser(-3,3);
  hTrackXresid_vs_module->GetYaxis()->SetTitle("Track X residual (mm)");
  hTrackXresid_vs_module->GetXaxis()->SetTitle("module");
  
  hTrackXresid_vs_module->Draw("colz");

  c1->cd(5);

  hTrackYresid_vs_module->SetTitle("");
  hTrackYresid_vs_module->GetYaxis()->SetRangeUser(-3,3);
  hTrackYresid_vs_module->GetYaxis()->SetTitle("Track Y residual (mm)");
  hTrackYresid_vs_module->GetXaxis()->SetTitle("module");
  
  hTrackYresid_vs_module->Draw("colz");

  c1->cd(3);

  TH1D *hXresid = new TH1D("hXresid","",1000,-3,3);
  TH1D *hYresid = new TH1D("hYresid","",1000,-3,3);

  TString schi2cut;
  TCut cut = schi2cut.Format("TrackChi2NDF<%g",chi2cut).Data();
  
  Tout->Project("hXresid","HitXresid",cut);
  Tout->Project("hYresid","HitYresid",cut);

  double sig = 0.5*(hXresid->GetRMS()+hYresid->GetRMS());

  TString stitle;

  int binmax = hXresid->GetMaximumBin();

  int binlow = binmax, binhigh = binmax;
  
  while( hXresid->GetBinContent(binlow--)>0.6*hXresid->GetBinContent(binmax) ){};
  while( hXresid->GetBinContent(binhigh++)>0.6*hXresid->GetBinContent(binmax) ){};
  
  hXresid->GetXaxis()->SetTitle(stitle.Format("Track X residual (mm), #chi^{2}/dof<%5.3g",chi2cut));
  hXresid->Draw();
  
  hXresid->Fit("gaus","","",hXresid->GetBinCenter(binlow),hXresid->GetBinCenter(binhigh));
  
  c1->cd(6);

  binmax = hYresid->GetMaximumBin();

  binlow = binmax; binhigh = binmax;
  while( hYresid->GetBinContent(binlow--)>0.6*hYresid->GetBinContent(binmax) ){};
  while( hYresid->GetBinContent(binhigh++)>0.6*hYresid->GetBinContent(binmax) ){};
  
  hYresid->GetXaxis()->SetTitle(stitle.Format("Track Y residual (mm), #chi^{2}/dof<%5.3g",chi2cut));
  hYresid->Draw();
  hYresid->Fit("gaus","","",hYresid->GetBinCenter(binlow),hYresid->GetBinCenter(binhigh));

  c1->Print(pdffilename);

  hTrackXeresid_vs_layer->SetYTitle("Track X excl. residual (mm)");
  hTrackYeresid_vs_layer->SetYTitle("Track Y excl. residual (mm)");
  hTrackXeresid_vs_layer->SetXTitle("Layer");
  hTrackYeresid_vs_layer->SetXTitle("Layer");

  hTrackXeresid_vs_module->SetYTitle("Track X excl. residual (mm)");
  hTrackYeresid_vs_module->SetYTitle("Track Y excl. residual (mm)");
  hTrackXeresid_vs_module->SetXTitle("Module");
  hTrackYeresid_vs_module->SetXTitle("Module");

  c1->cd(1);
  hTrackXeresid_vs_layer->Draw("colz");

  c1->cd(2);
  hTrackXeresid_vs_module->Draw("colz");

  c1->cd(4);
  hTrackYeresid_vs_layer->Draw("colz");

  c1->cd(5);
  hTrackYeresid_vs_module->Draw("colz");

  TH1D *hXeresid = new TH1D("hXeresid","",1000,-6,6);
  TH1D *hYeresid = new TH1D("hYeresid","",1000,-6,6);
  
  Tout->Project("hXeresid","HitXresidE",cut);
  Tout->Project("hYeresid","HitYresidE",cut);

  binmax = hXeresid->GetMaximumBin();
  binlow = binmax; binhigh = binmax;

  while( hXeresid->GetBinContent(binlow--)>0.6*hXeresid->GetBinContent(binmax) ){};
  while( hXeresid->GetBinContent(binhigh++)>0.6*hXeresid->GetBinContent(binmax) ){};

  
  c1->cd(3);
  hXeresid->Fit("gaus","","",hXeresid->GetBinCenter(binlow),hXeresid->GetBinCenter(binhigh));
  hXeresid->SetXTitle(stitle.Format("Track X excl. resid (mm), #chi^{2}/dof<%5.3g",chi2cut));
  hXeresid->Draw();

  binmax = hYeresid->GetMaximumBin();
  binlow = binmax; binhigh = binmax;

  while( hYeresid->GetBinContent(binlow--)>0.6*hYeresid->GetBinContent(binmax) ){};
  while( hYeresid->GetBinContent(binhigh++)>0.6*hYeresid->GetBinContent(binmax) ){};

  c1->cd(6);
  hYeresid->Fit("gaus","","",hYeresid->GetBinCenter(binlow),hYeresid->GetBinCenter(binhigh));
  hYeresid->SetXTitle(stitle.Format("Track Y excl. resid (mm), #chi^{2}/dof<%5.3g",chi2cut));
  hYeresid->Draw();
  
  
  c1->Print(pdffilename);
  
  TCanvas *c3 = new TCanvas("c3","c3",2400,1200);
  c3->Divide(2,1,.001,.001);

  lmargin=0.1;
  rmargin=0.05;
  tmargin=0.06;
  
  c3->cd(1);

  gPad->SetLeftMargin(lmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);
  
  hXresid->Draw();
  c3->cd(2);

  gPad->SetLeftMargin(lmargin);
  gPad->SetRightMargin(rmargin);
  gPad->SetBottomMargin(bmargin);
  gPad->SetTopMargin(tmargin);
  
  hYresid->Draw();

  
  
  TCanvas *c2 = new TCanvas("c2","c2",2400,1200);

  lmargin = 0.2;
  bmargin = 0.12;
  rmargin = 0.18;
  tmargin = 0.06;

  gStyle->SetNdivisions(510,"XYZ");
  gStyle->SetLabelSize(.06,"XYZ");
  gStyle->SetTitleSize(.06,"XYZ");

  gROOT->ForceStyle();
  
  c2->Divide(nlayers,1,.001,.001);

  for( int ilayer = 0; ilayer<nlayers; ilayer++ ){
    TH2D *htemp = ( (TH2D*) (*heff_layers)[ilayer] );

    c2->cd(ilayer+1);

    gPad->SetLeftMargin(lmargin);
    gPad->SetRightMargin(rmargin);
    gPad->SetBottomMargin(bmargin);
    gPad->SetTopMargin(tmargin);

    htemp->SetTitleSize(.06,"XYZ");
    htemp->SetLabelSize(.06,"XYZ");
    htemp->SetTitleOffset(1.4,"Y");
    htemp->SetNdivisions(510,"YZ");
    htemp->GetXaxis()->SetTitle("Y (mm)");
    htemp->GetYaxis()->SetTitle("X (mm)");
    htemp->Draw("colz");
  }

  c2->Print(pdffilename);

  for( int ilayer = 0; ilayer<nlayers; ilayer++ ){
    TH2D *htemp = ( (TH2D*) (*hxyhit_layers)[ilayer] );

    c2->cd(ilayer+1);

    gPad->SetLeftMargin(lmargin);
    gPad->SetRightMargin(rmargin);
    gPad->SetBottomMargin(bmargin);
    gPad->SetTopMargin(tmargin);

    htemp->SetTitleSize(.06,"XYZ");
    htemp->SetLabelSize(.06,"XYZ");
    htemp->SetTitleOffset(1.4,"Y");
    htemp->SetNdivisions(510,"YZ");
    htemp->GetXaxis()->SetTitle("Y (mm)");
    htemp->GetYaxis()->SetTitle("X (mm)");
    htemp->Draw("colz");
  }

  c2->Print(pdffilename);

  //Make a set of canonical plots module-by-module:
  TCanvas *c4 = new TCanvas("c4","c4",1760,1320);

  //1. Number of strips fired X and Y (same plot, overlay)
  //2. Number of clusters (track passed through)
  //3. Cluster width X and Y (same plot, overlay)
  //4. Max time bin X and Y (same plot, overlay)
  //5. ADC asymmetry
  //6. X-Y time difference
  //7. ADC sum for max X and Y strips (same plot, overlay)
  //8. Max ADC sample for X and Y strips (same plot, overlay)
  //9. Cluster ADC sum X and Y (same plot, overlay)
  //10. seed pixel correlation coeff. and cluster corr. coeff. (same plot, overlay)
  //11. X inclusive residuals with fit 
  //12. Y inclusive residuals with fit 
  //13. X exclusive residuals with fit
  //14. Y exclusive residuals with fit
  //15. Hit map
  //16. Efficiency map

  TH2D *hNstripsX_module, *hNstripsY_module, *hNclust_module, *hclustwidth_x_module,*hclustwidth_y_module;
  TH2D *hmaxtimebin_xstrip_module,*hmaxtimebin_ystrip_module;
  TH2D *hADCasym_vs_module,*hdT_vs_module;
  TH2D *hADCsum_Xstrip_max_module,*hADCsum_Ystrip_max_module;
  TH2D *hADCsampmax_Xstrip_module,*hADCsampmax_Ystrip_module;
  TH2D *hADCprodXYstrip_max_module,*hADCprodXYsamp_max_module,*hADCprodXYclust_module;
  TH2D *hADCsumXclust_module,*hADCsumYclust_module;
  TH2D *hStrip_maxcor_module,*hClust_corr_module;

  fin->GetObject("hNstripsX_module",hNstripsX_module);
  fin->GetObject("hNstripsY_module",hNstripsY_module);
  fin->GetObject("hNclust_module",hNclust_module);
  fin->GetObject("hclustwidth_x_module",hclustwidth_x_module);
  fin->GetObject("hclustwidth_y_module",hclustwidth_y_module);
  fin->GetObject("hmaxtimebin_xstrip_module",hmaxtimebin_xstrip_module);
  fin->GetObject("hmaxtimebin_ystrip_module",hmaxtimebin_ystrip_module);
  fin->GetObject("hADCasym_vs_module",hADCasym_vs_module);
  fin->GetObject("hdT_vs_module",hdT_vs_module);
  fin->GetObject("hADCsum_Xstrip_max_module",hADCsum_Xstrip_max_module);
  fin->GetObject("hADCsum_Ystrip_max_module",hADCsum_Ystrip_max_module);
  fin->GetObject("hADCsampmax_Xstrip_module",hADCsampmax_Xstrip_module);
  fin->GetObject("hADCsampmax_Ystrip_module",hADCsampmax_Ystrip_module);
  fin->GetObject("hADCprodXYstrip_max_module",hADCprodXYstrip_max_module);
  fin->GetObject("hADCprodXYsamp_max_module",hADCprodXYsamp_max_module);
  fin->GetObject("hADCprodXYclust_module",hADCprodXYclust_module);
  fin->GetObject("hADCsumXclust_module",hADCsumXclust_module);
  fin->GetObject("hADCsumYclust_module",hADCsumYclust_module);
  fin->GetObject("hStrip_maxcor_module",hStrip_maxcor_module);
  fin->GetObject("hClust_corr_module",hClust_corr_module);
  
  //TClonesArray *hxyhit_modules = new TClonesArray("TH2D",nmodules);
  //TClonesArray *heff_modules = new TClonesArray("TH2D",nmodules);

  //THStack *hoverlay = new THStack("hoverlay","");
  
  for( int imodule=0; imodule<nmodules; imodule++ ){
    //int module = imodule-1;
    TString hnametemp;
    hnametemp.Form("hxyhit_module%d",imodule);
    
    TH2D *hxytemp;
    fin->GetObject( hnametemp, hxytemp );
    //new( (TH2D*) (*hxyhit_modules)[imodule] ) TH2D( *htemp );

    hnametemp.Form("heff_module%d",imodule);

    TH2D *hefftemp;
    fin->GetObject( hnametemp, hefftemp );
    //new( (TH2D*) (*heff_modules)[imodule] ) TH2D( *htemp );

    c4->Clear();
    c4->Divide(4,4,.001,.001);
    c4->Update();

    TH1D *htemp1,*htemp2,*htemp3;
    
    c4->cd(1);

    htemp1 = hNstripsX_module->ProjectionX( hnametemp.Format("hnstripx_mod%d",imodule), imodule+1,imodule+1 );
    htemp2 = hNstripsY_module->ProjectionX( hnametemp.Format("hnstripy_mod%d",imodule), imodule+1,imodule+1 );

    htemp1->SetLineColor(2);
    //   htemp1->SetLineStyle(1);
    htemp2->SetLineColor(4);
    //htemp2->SetLineStyle(10);

    if( htemp1->GetMaximum() > htemp2->GetMaximum() ){
      htemp1->Draw();
      htemp1->SetTitle("");
      htemp1->GetXaxis()->SetRangeUser(-0.5,20.5);
      htemp2->Draw("SAME");
    } else {
      htemp2->Draw();
      htemp2->SetTitle("");
      htemp2->GetXaxis()->SetRangeUser(-0.5,20.5);
      htemp1->Draw("SAME");
    }

    c4->cd(2);
    htemp1 = hNclust_module->ProjectionX( hnametemp.Format("hnclust_mod%d",imodule), imodule+1,imodule+1);
    htemp1->SetTitle("");
    htemp1->Draw();
    htemp1->GetXaxis()->SetRangeUser(-0.5,20.5);

    c4->cd(3);
    htemp1 = hclustwidth_x_module->ProjectionX( hnametemp.Format("hclwidthx_mod%d",imodule), imodule+1,imodule+1);
    htemp2 = hclustwidth_y_module->ProjectionX( hnametemp.Format("hclwidthy_mod%d",imodule), imodule+1,imodule+1);

    htemp1->SetLineColor(2);
    //   htemp1->SetLineStyle(1);
    htemp2->SetLineColor(4);
    //htemp2->SetLineStyle(10);
    
    if( htemp1->GetMaximum() > htemp2->GetMaximum() ){
      htemp1->Draw();
      htemp2->Draw("SAME");
    } else {
      htemp2->Draw();
      htemp1->Draw("SAME");
    }

    c4->cd(4);

    htemp1 = hmaxtimebin_xstrip_module->ProjectionX(hnametemp.Format("hmaxbinx_mod%d",imodule),imodule+1,imodule+1);
    htemp2 = hmaxtimebin_ystrip_module->ProjectionX(hnametemp.Format("hmaxbiny_mod%d",imodule),imodule+1,imodule+1);

    htemp1->SetLineColor(2);
    //   htemp1->SetLineStyle(1);
    htemp2->SetLineColor(4);
    //htemp2->SetLineStyle(10);
    
    if( htemp1->GetMaximum() > htemp2->GetMaximum() ){
      htemp1->Draw();
      htemp2->Draw("SAME");
    } else {
      htemp2->Draw();
      htemp1->Draw("SAME");
    }

    c4->cd(5);
    htemp1 = hADCasym_vs_module->ProjectionY(hnametemp.Format("hADCasym_mod%d",imodule),imodule+1,imodule+1);
    htemp1->Draw();

    c4->cd(6);
    htemp1 = hdT_vs_module->ProjectionY(hnametemp.Format("hdT_mod%d",imodule),imodule+1,imodule+1);
    htemp1->Draw();
    
    c4->cd(7);

    hxytemp->Draw("colz");

    c4->cd(8);

    hefftemp->SetMinimum(0);
    hefftemp->SetMaximum(1);
    hefftemp->Draw("colz");

    c4->cd(9);

    htemp1 = hADCsampmax_Xstrip_module->ProjectionX(hnametemp.Format("hsxmax_mod%d",imodule),imodule+1,imodule+1);
    htemp2 = hADCsampmax_Ystrip_module->ProjectionX(hnametemp.Format("hsymax_mod%d",imodule),imodule+1,imodule+1);
    htemp3 = hADCprodXYsamp_max_module->ProjectionX(hnametemp.Format("hsxymax_mod%d",imodule),imodule+1,imodule+1);
    htemp1->SetLineColor(2);
    //   htemp1->SetLineStyle(1);
    htemp2->SetLineColor(4);
    //htemp2->SetLineStyle(10);
    htemp3->SetLineColor(kGreen+1);
    
    THStack *hstack_samp = new THStack("hstack_samp","");
    hstack_samp->Add(htemp1);
    hstack_samp->Add(htemp2);
    hstack_samp->Add(htemp3);

    hstack_samp->Draw("nostack");
    hstack_samp->GetXaxis()->SetRangeUser(0,2.75e3);
    
    c4->cd(10);

    htemp1 = hADCsum_Xstrip_max_module->ProjectionX(hnametemp.Format("hstrxmax_mod%d",imodule),imodule+1,imodule+1);
    htemp2 = hADCsum_Ystrip_max_module->ProjectionX(hnametemp.Format("hstrymax_mod%d",imodule),imodule+1,imodule+1);
    htemp3 = hADCprodXYstrip_max_module->ProjectionX(hnametemp.Format("hstrxymax_mod%d",imodule),imodule+1,imodule+1);

    htemp1->SetLineColor(2);
    //   htemp1->SetLineStyle(1);
    htemp2->SetLineColor(4);
    //htemp2->SetLineStyle(10);
    htemp3->SetLineColor(kGreen+1);
    
    THStack *hstack_strip = new THStack("hstack_strip","");
    hstack_strip->Add(htemp1);
    hstack_strip->Add(htemp2);
    hstack_strip->Add(htemp3);

    hstack_strip->Draw("nostack");
    hstack_strip->GetXaxis()->SetRangeUser(0.0,1.25e4);

    c4->cd(11);
    
    htemp1 = hADCsumXclust_module->ProjectionX(hnametemp.Format("hclx_mod%d",imodule),imodule+1,imodule+1);
    htemp2 = hADCsumYclust_module->ProjectionX(hnametemp.Format("hcly_mod%d",imodule),imodule+1,imodule+1);
    htemp3 = hADCprodXYclust_module->ProjectionX(hnametemp.Format("hclxy_mod%d",imodule),imodule+1,imodule+1);

    htemp1->SetLineColor(2);
    //   htemp1->SetLineStyle(1);
    htemp2->SetLineColor(4);
    //htemp2->SetLineStyle(10);
    htemp3->SetLineColor(kGreen+1);
    
    THStack *hstack_clust = new THStack("hstack_clust","");
    hstack_clust->Add(htemp1);
    hstack_clust->Add(htemp2);
    hstack_clust->Add(htemp3);

    hstack_clust->Draw("nostack");
    hstack_clust->GetXaxis()->SetRangeUser(0,2.5e4);

    c4->cd(12)->SetLogy();

    htemp1 = hStrip_maxcor_module->ProjectionX(hnametemp.Format("hstrcorr%d",imodule),imodule+1,imodule+1);
    htemp2 = hClust_corr_module->ProjectionX(hnametemp.Format("hclcorr%d",imodule),imodule+1,imodule+1);

    htemp1->SetLineColor(2);
    htemp2->SetLineColor(4);
    
    THStack *hstack_corr = new THStack("hstack_corr","");

    hstack_corr->Add(htemp1);
    hstack_corr->Add(htemp2);
    
    hstack_corr->Draw("nostack");

    c4->cd(13);
    htemp1 = hTrackXresid_vs_module->ProjectionY(hnametemp.Format("hresidx%d",imodule),imodule+1,imodule+1);

    htemp1->Draw();
    htemp1->GetXaxis()->SetRangeUser(-3.,3.);
    htemp1->Fit("gaus","","",-0.5*htemp1->GetRMS(),0.5*htemp1->GetRMS());
    
    c4->cd(15);

    htemp1 = hTrackXeresid_vs_module->ProjectionY(hnametemp.Format("heresidx%d",imodule),imodule+1,imodule+1);

    htemp1->Draw();
    htemp1->GetXaxis()->SetRangeUser(-3.,3.);
    htemp1->Fit("gaus","","",-0.5*htemp1->GetRMS(),0.5*htemp1->GetRMS());

    c4->cd(14);
    htemp1 = hTrackYresid_vs_module->ProjectionY(hnametemp.Format("hresidy%d",imodule),imodule+1,imodule+1);

    htemp1->Draw();
    htemp1->GetXaxis()->SetRangeUser(-3.,3.);
    htemp1->Fit("gaus","","",-0.5*htemp1->GetRMS(),0.5*htemp1->GetRMS());

    c4->cd(16);
    htemp1 = hTrackYeresid_vs_module->ProjectionY(hnametemp.Format("heresidy%d",imodule),imodule+1,imodule+1);

    htemp1->Draw();
    htemp1->GetXaxis()->SetRangeUser(-3.,3.);
    htemp1->Fit("gaus","","",-0.5*htemp1->GetRMS(),0.5*htemp1->GetRMS());
    
    
    gPad->Modified();
    c4->Update();
    c4->Print(pdffilename);
  }
  
  
  // heff_layer1->SetTitleSize(.06,"XYZ");
  // heff_layer1->SetLabelSize(.06,"XYZ");
  // heff_layer1->SetTitleOffset(1.4,"Y");
  // heff_layer1->SetNdivisions(510,"YZ");
  
  // heff_layer2->SetTitleSize(.06,"XYZ");
  // heff_layer2->SetLabelSize(.06,"XYZ");
  // heff_layer2->SetTitleOffset(1.4,"Y");
  // heff_layer2->SetNdivisions(510,"YZ");

  // heff_layer3->SetTitleSize(.06,"XYZ");
  // heff_layer3->SetLabelSize(.06,"XYZ");
  // heff_layer3->SetTitleOffset(1.4,"Y");
  // heff_layer3->SetNdivisions(510,"YZ");

  // hxyhit_layer0->SetTitleSize(.06,"XYZ");
  // hxyhit_layer0->SetLabelSize(.06,"XYZ");
  // hxyhit_layer0->SetTitleOffset(1.4,"Y");
  // hxyhit_layer0->SetNdivisions(510,"YZ");

  // hxyhit_layer1->SetTitleSize(.06,"XYZ");
  // hxyhit_layer1->SetLabelSize(.06,"XYZ");
  // hxyhit_layer1->SetTitleOffset(1.4,"Y");
  // hxyhit_layer1->SetNdivisions(510,"YZ");
  
  // hxyhit_layer2->SetTitleSize(.06,"XYZ");
  // hxyhit_layer2->SetLabelSize(.06,"XYZ");
  // hxyhit_layer2->SetTitleOffset(1.4,"Y");
  // hxyhit_layer2->SetNdivisions(510,"YZ");

  // hxyhit_layer3->SetTitleSize(.06,"XYZ");
  // hxyhit_layer3->SetLabelSize(.06,"XYZ");
  // hxyhit_layer3->SetTitleOffset(1.4,"Y");
  // hxyhit_layer3->SetNdivisions(510,"YZ");

  // hxyhit_layer0->GetXaxis()->SetTitle("Y (mm)");
  // hxyhit_layer0->GetYaxis()->SetTitle("X (mm)");

  // hxyhit_layer1->GetXaxis()->SetTitle("Y (mm)");
  // hxyhit_layer1->GetYaxis()->SetTitle("X (mm)");

  // hxyhit_layer2->GetXaxis()->SetTitle("Y (mm)");
  // hxyhit_layer2->GetYaxis()->SetTitle("X (mm)");

  // hxyhit_layer3->GetXaxis()->SetTitle("Y (mm)");
  // hxyhit_layer3->GetYaxis()->SetTitle("X (mm)");
  
  // c2->cd(1);

  

  

  // c2->cd(2);

  // gPad->SetLeftMargin(lmargin);
  // gPad->SetRightMargin(rmargin);
  // gPad->SetBottomMargin(bmargin);
  // gPad->SetTopMargin(tmargin);

  // heff_layer1->GetXaxis()->SetTitle("Y (mm)");
  // heff_layer1->GetYaxis()->SetTitle("X (mm)");
  // heff_layer1->Draw("colz");

  // c2->cd(3);

  // gPad->SetLeftMargin(lmargin);
  // gPad->SetRightMargin(rmargin);
  // gPad->SetBottomMargin(bmargin);
  // gPad->SetTopMargin(tmargin);

  // heff_layer2->GetXaxis()->SetTitle("Y (mm)");
  // heff_layer2->GetYaxis()->SetTitle("X (mm)");
  // heff_layer2->Draw("colz");

  // c2->cd(4);

  // gPad->SetLeftMargin(lmargin);
  // gPad->SetRightMargin(rmargin);
  // gPad->SetBottomMargin(bmargin);
  // gPad->SetTopMargin(tmargin);

  // heff_layer3->GetXaxis()->SetTitle("Y (mm)");
  // heff_layer3->GetYaxis()->SetTitle("X (mm)");
  // heff_layer3->Draw("colz");

  // c2->Print(pdffilename);

  // c2->cd(1);

  // hxyhit_layer0->Draw("col");

  // c2->cd(2);

  // hxyhit_layer1->Draw("col");

  // c2->cd(3);

  // hxyhit_layer2->Draw("col");

  // c2->cd(4);

  // hxyhit_layer3->Draw("col");
  //  hxyhit_layer0->Draw("colz");

  //c2->Print(pdffilename);
  
  c1->Print(closefilename);
}
