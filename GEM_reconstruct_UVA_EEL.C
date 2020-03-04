#include "TTree.h"
#include "TChain.h"
#include "GEMHit_tree_UVA_EEL.C"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TF1.h"
#include "TF2.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include "TString.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TPad.h"
#include "TROOT.h"
#include "TMath.h"
#include "TGraph2DErrors.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMarker.h"
#include "TClonesArray.h"
#include "TLine.h"
//#include "TApplication.h"

//TODO: don't hard-code the number of APV25 samples. But for this we need a better decoded hit format anyway

//Note that our definition of "X" and "Y" will differ from simprox.cpp method. What we are calling the "X" direction is the "vertical" (long) axis of the
//layer, with more strips (planeID == 1)
//What we are calling "Y" is the horizontal dimension with fewer strips (planeID == 0).

//We need to modify this so that we no longer assume that every module has the same numbers of strips or strip orientations:
//The only assumption we will retain is the assumption that each module has two non-parallel strip orientations, generically denoted "U" and "V"
//We will also make the strip pitch along each direction a configurable parameter:
int nstripsx = 1280;
int nstripsy = 1024;
int nmodules = 12;
int nlayers  = 4;
int nADCsamples = 6;

double strip_pitch = 0.4; //mm

const double PI = TMath::Pi();

//These default values are guesses:
//int maxnhitspercluster=10;
int maxnstripspercluster=7;
int maxnstripXpercluster=9;
int maxnstripYpercluster=7;
//double clustersigma=0.4; //mm
//double clustertau=50.0; //ns
//double sigmasample=20.0; //individual sample noise width
//double localmaxthreshold_nsigma=3.0; //number of sigmas to create a new local max/hit:
//double sigmaT_2Dmatch=26.0; //guess
//double sigmaADC_2Dmatch=18.0; //guess
double cluster2Dmatch_asymcut=0.5; //dimensionless
double cluster2Dmatch_tcut = 25.0; //ns
double maxstripcorthreshold = 0.0;
//double maxADCXYthreshold = 3000.0; 
double stripcorthreshold = 0.0; //correlation coefficient must be positive:
double clustcorthreshold = 0.0;

double thresh_maxsample = 150.0; //max ADC sample on a strip must exceed this to be included
double thresh_stripsum =  300.0; //threshold on the sum of all samples on a strip:
double thresh_clustersum = 1000.0; //threshold on the summed ADCs of all strips in a cluster:

double sigma_hitpos=0.15; //mm
double TrackMaxSlopeX = 1.0; //
double TrackMaxSlopeY = 0.5;
double TrackChi2Cut=1.e4; //Let's initially assume the chambers are aligned to within 1 cm = 10 mm. Then the max residual at any plane should be about 1 cm. the chi^2 contribution of said hit would be
double TrackFindingMaxRadius=20.0; //max radius from projection of track to find hits for adding to track:
int maxnhitcombinations=100000; //need to make user-configurable; max number of possible hit combinations to attempt tracking; give up on "too noisy" events
//flags to decide whether hit template function cluster width and timing parameters are fixed or variable:
// int varyclustersigma=0;
// int varyclustertau=0;
// int fittemplatefunctions=0;

//Default values for walk correction parameters:
//tstrip_cor = tstrip + pow( ADC0/ADCstrip, exp )
double walkcor_mean_const = 1.87334;
double walkcor_mean_ADC0=2818.0;
double walkcor_mean_exp =1.133;

double walkcor_sigma_const = 1.725;
double walkcor_sigma_ADC0=6448.0;
double walkcor_sigma_exp=0.78896;

//Width of strip timing cut for addition to cluster, in standard deviations:
double tstripcut_nsigma = 5.0;

//Walk correction for strip times:
TF1 *walkcor_mean_func = new TF1("walkcor_mean_func","[0]-pow([1]/x,[2])",250.0,2e4);
TF1 *walkcor_sigma_func = new TF1("walkcor_sigma_func","[0]+pow([1]/x,[2])",250.0,2e4);

//(optional) histograms containing cluster mapping functions for improved hit reconstruction:
bool usehitmaps=false;
TH2D *hXmom_int;
TH2D *hYmom_int; 

struct moduledata_t { //for now, we will keep the "x" and "y" notation for these data structures:
  int modindex,layerindex;
  set<int> xstrips;
  set<int> ystrips;
  map<int,double> ADCsum_xstrips;
  map<int,double> ADCsum_ystrips;
  map<int,vector<double> > ADCsamp_xstrips;
  map<int,vector<double> > ADCsamp_ystrips;

  map<int,int> Bestmatch_xstrips; //Best y match by correlation coefficient for x strips:
  map<int,int> Bestmatch_ystrips; //Best x match by correlation coefficient for y strips:

  map<int,double> BestCor_xstrips; //Best correlation coefficient, x strips
  map<int,double> BestCor_ystrips; //Best correlation coefficient, y strips

  //list of x and y strips, filtered by correlation coefficient match. This array is guaranteed to have equal size in x and y:
  set<int> xstrips_filtered;
  set<int> ystrips_filtered;
  
  //mapping by strip, max. individual ADC sample.
  map<int,double> ADCmax_xstrips;
  map<int,double> ADCmax_ystrips;

  //mapping by strip, index of max individual ADCsample:
  map<int,int> isampmax_xstrips;
  map<int,int> isampmax_ystrips;
  

  map<int,double> Tmean_xstrips;
  map<int,double> Tsigma_xstrips;

  map<int,double> Tmean_ystrips;
  map<int,double> Tsigma_ystrips;

  map<int,double> Tmean_xstrips_walkcor;
  map<int,double> Tsigma_xstrips_walkcor;
  map<int,double> Tmean_ystrips_walkcor;
  map<int,double> Tsigma_ystrips_walkcor;

  
  
};

struct clusterdata_t { //1D and 2D clustering results by module:

  int modindex,layerindex;
  
  int nclust2D;
  vector<int> itrack_clust2D;
  vector<int> ixclust2D;
  vector<int> iyclust2D;
  vector<int> nstripx2D;
  vector<int> nstripy2D;
  vector<double> xclust2D;
  vector<double> yclust2D;
  vector<double> xclust2Dcorr;
  vector<double> yclust2Dcorr;
  vector<double> Eclust2D;
  vector<double> tclust2D;
  vector<double> tclust2Dwalkcorr; //walk-corrected cluster time
  vector<double> dEclust2D;
  vector<double> dtclust2D;
  vector<double> dtclust2Dwalkcorr; //"uncertainty" of weighted average cluster time:
  vector<double> CorrCoeff2D; //Correlation coefficient
  vector<double> CorrCoeffMaxStrips; //Correlation coefficient between X and Y strips used to "seed" the cluster:
  vector<double> xglobal2D;
  vector<double> yglobal2D;
  vector<double> zglobal2D; //calculate these ONCE:
  vector<double> xmom2D; // xmean - xstrip max (use to refine hit position reconstruction)
  vector<double> ymom2D; // ymean - ystrip max (use to refine hit position reconstruction)
  vector<bool> keepclust2D; //flag to include or not include cluster in tracking analysis: default to true; set by "prune" filtering algorithm:
  vector<vector<double> > ADCsamp_xclust; // Cluster-summed ADC time samples, X strips:
  vector<vector<double> > ADCsamp_yclust; // Cluster-summed ADC time samples, Y strips;
  
  
  int nclustx;
  vector<int> nstripx;
  vector<int> ixstriplo;
  vector<int> ixstriphi;
  vector<int> ixstripmax;
  vector<double> xmean;
  vector<double> xsigma;
  vector<double> totalchargex;
  vector<double> txmean;
  vector<double> txsigma;

  int nclusty;
  vector<int> nstripy;
  vector<int> iystriplo;
  vector<int> iystriphi;
  vector<int> iystripmax;
  vector<double> ymean;
  vector<double> ysigma;
  vector<double> totalchargey;
  vector<double> tymean;
  vector<double> tysigma;

  // int nhitx1D;
  // vector<int>    clustidx_xhit;
  // vector<double> x0hit1D;
  // vector<double> tx0hit1D;
  // vector<double> Axhit1D;
  // vector<double> chi2ndf_hitx1D;
  // vector<double> sigx0hit1D; //for fits with variable cluster width
  // vector<double> taux0hit1D; //for fits with variable cluster time constant
  
  // int nhity1D;
  // vector<int>    clustidx_yhit;
  // vector<double> y0hit1D;
  // vector<double> ty0hit1D;
  // vector<double> Ayhit1D;
  // vector<double> chi2ndf_hity1D;
  // vector<double> sigy0hit1D; //for fits with variable cluster width
  // vector<double> tauy0hit1D; //for fits with variable cluster time constant

  // int nhit2D;
  // vector<double> xhit2D;
  // vector<double> yhit2D;
  // vector<double> thit2D;
  // vector<double> dthit2D;
  // vector<double> Ahit2D;
  // vector<double> dAhit2D;
  
  // vector<double> ixhit2D;
  // vector<double> iyhit2D;
  // vector<double> chi2match2D;
  
  
};

struct trackdata_t {
  int ntracks;
  vector<int> nhitsontrack;
  vector<vector<int> > modlist_track; //module index of hits on track
  vector<vector<int> > hitlist_track; //index in 2D cluster array within module of hits on track
  vector<vector<double> > residx_hits; //X residuals of track wrt hits:
  vector<vector<double> > residy_hits; //Y residuals of track wrt hits
  vector<vector<double> > eresidx_hits; //X residuals excluding the hit in question from the track
  vector<vector<double> > eresidy_hits; //Y residuals excluding the hit in question from the track
  vector<double> Xtrack,Ytrack,Xptrack,Yptrack;
  vector<double> Chi2NDFtrack;
  //  vector<double> NDFtrack; 
  
};

//Detector maps:
map<int,int> mod_layer; //key = module index, value = layer (0..3) logical layer index for modules.
map<int,double> mod_x0; //key = module index, value = x center position
map<int,double> mod_y0; //key = module index, value = y center position
map<int,double> mod_z0; //key = module index, value = z position
map<int,double> mod_ax; //key = module index, value = X axis rotation (yaw)
map<int,double> mod_ay; //key = module index, value = Y axis rotation (pitch)
map<int,double> mod_az; //key = module index, value = Z axis rotation (roll)
map<int,int> mod_nstripsu; //key = module index, value = total number of strips along "U" direction 
map<int,int> mod_nstripsv; //key = module index, value = total number of strips along "V" direction
map<int,double> mod_uangle; //key = module index, value = angle between module x axis and "U" direction
map<int,double> mod_vangle; //key = module index, value = angle between modyle x axis and "V" direction
map<int,int>    mod_uplaneID; //might as well make this configurable as well:
map<int,int>    mod_vplaneID;
map<int,double> mod_Pxu;    //cos(uangle);
map<int,double> mod_Pyu;    //sin(uangle);
map<int,double> mod_Pxv;    //cos(vangle);
map<int,double> mod_Pyv;    //cos(vangle);
map<int,double> mod_ustrip_pitch; //strip pitch in "U" direction
map<int,double> mod_vstrip_pitch; //strip pitch in "V" direction
map<int,double> mod_Lx;     //module physical extent in "X" direction
map<int,double> mod_Ly;     //module physical extent in "Y" direction
map<int,TRotation> mod_Rot; //module rotation matrix: Compute ONCE!
map<int,TRotation> mod_Rotinv; //inverse of module rotation matrix:

// Tracking layer combinatorics: populate these arrays once so we don't do it every event:
//For each possible number of layers from 3 up to the total number of layers, we list all possible combinations of n layers
map<int,vector<vector<int> > > layercombos;


vector<double> zavg_layer;

// double HitTemplateFunc( double *x, double *par ){
//   double X = x[0];
//   double T = x[1];

//   int Nhit = int(par[0]);

//   double ADCsum = 0.0;

//   int ipar=1;
//   for( int ihit=0; ihit<Nhit; ihit++ ){  
//     double A = par[5*ihit+1];
//     double x0 = par[5*ihit+2];
//     double tmax = par[5*ihit+3];
//     double sigma = par[5*ihit+4]; //ordinarily fixed
//     double tau   = par[5*ihit+5]; //ordinarily fixed
    
//     double t0 = tmax - tau;
    
//     //sum the contributions of individual hits:
//     ADCsum += A/(1.0+pow((X-x0)/sigma,2))*TMath::Max(0.0,(T-t0)/tau*exp(-(T-t0)/tau+1.0));
//   }
  
//   return ADCsum;
// }

//TF2 for hit template function fitting:
//TF2 *HitFunc = new TF2("HitFunc",HitTemplateFunc, -0.55*nstripsx,0.55*nstripsx, 0.0,150.0,5*maxnhitspercluster);

void filter_strips_by_module( moduledata_t &mod_data, double cor_threshold=0.0 ){

  mod_data.Bestmatch_xstrips.clear();
  mod_data.Bestmatch_ystrips.clear();
  mod_data.xstrips_filtered.clear();
  mod_data.ystrips_filtered.clear();
  
  for( set<int>::iterator ix=mod_data.xstrips.begin(); ix!=mod_data.xstrips.end(); ++ix ){ //Compute a correlation coefficient for all possible combinations of X and Y strips:
    int ixstrip=*ix;
    
    mod_data.Bestmatch_xstrips[ixstrip] = -1;
    mod_data.BestCor_xstrips[ixstrip] = -1.1;
    
    for( set<int>::iterator iy=mod_data.ystrips.begin(); iy!=mod_data.ystrips.end(); ++iy ){
      int iystrip=*iy;

      if( mod_data.Bestmatch_ystrips.find(iystrip) == mod_data.Bestmatch_ystrips.end() ){ //first time considering this y strip
	mod_data.Bestmatch_ystrips[iystrip] = -1;
	mod_data.BestCor_ystrips[iystrip] = -1.1;
      }
      
      double sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;

      for( int isamp=0; isamp<6; isamp++ ){
	double ADCx = mod_data.ADCsamp_xstrips[ixstrip][isamp];
	double ADCy = mod_data.ADCsamp_ystrips[iystrip][isamp];
	
	sumx += ADCx;
	sumy += ADCy;
	sumx2 += pow(ADCx,2);
	sumy2 += pow(ADCy,2);
	sumxy += ADCx*ADCy;
      }

      double mux = sumx/6.0;
      double muy = sumy/6.0;
      double varx = sumx2/6.0-pow(mux,2);
      double vary = sumy2/6.0-pow(muy,2);
      double sigx = sqrt(varx);
      double sigy = sqrt(vary);

      double ccor = ( sumxy - 6.0*mux*muy )/( 6.0*sigx*sigy );

      if( mod_data.Bestmatch_xstrips[ixstrip] < 0 || ccor > mod_data.BestCor_xstrips[ixstrip] ){
	mod_data.Bestmatch_xstrips[ixstrip] = iystrip;
	mod_data.BestCor_xstrips[ixstrip] = ccor;
      }

      if( mod_data.Bestmatch_ystrips[iystrip] < 0 || ccor > mod_data.BestCor_ystrips[iystrip] ){
	mod_data.Bestmatch_ystrips[iystrip] = ixstrip;
	mod_data.BestCor_ystrips[iystrip] = ccor;
      }
    }
  }

  // Every x strip and y strip will have at least one "best match" assigned to it. We loop over whichever array is SMALLER.
  // Each strip in the SMALLER array could be associated with (in principle) more than one strip in the LARGER array
  // if( mod_data.xstrips.size() > mod_data.ystrips.size() ){ //more X strips: loop over Y strips:
  //   for( set<int>::iterator iy=mod_data.ystrips.begin(); iy!=mod_data.ystrips.end(); ++iy ){
  //     int iystrip = *iy;

  //     int bestx = mod_data.Bestmatch_ystrips[iystrip];

  //     if( mod_data.BestCor_ystrips[iystrip] > stripcorthreshold ){ //match with correlation coefficient above threshold:
  // 	//We always insert one x strip with one y strip:
  // 	mod_data.xstrips_filtered.insert( bestx );
  // 	mod_data.ystrips_filtered.insert( iystrip );
  //     }
  //   }
  // } else { //loop over x hits:
  //   for( set<int>::iterator ix=mod_data.xstrips.begin(); ix!=mod_data.xstrips.end(); ++ix ){
  //     int ixstrip = *ix;
  //     int besty = mod_data.Bestmatch_xstrips[ixstrip];

  //     if( mod_data.BestCor_xstrips[ixstrip] > stripcorthreshold ){ //match with correlation coefficient above threshold:
  // 	mod_data.xstrips_filtered.insert( ixstrip );
  // 	mod_data.ystrips_filtered.insert( besty );
  //     }
  //   }
  // }

  if( mod_data.xstrips.size() > 0 && mod_data.ystrips.size() > 0 ){
    for( set<int>::iterator iy=mod_data.ystrips.begin(); iy!=mod_data.ystrips.end(); ++iy ){
      int iystrip = *iy;
      int bestx = mod_data.Bestmatch_ystrips[iystrip];
      if( mod_data.BestCor_ystrips[iystrip] > stripcorthreshold ){ //match with correlation coefficient above threshold:
	mod_data.xstrips_filtered.insert( bestx );
	mod_data.ystrips_filtered.insert( iystrip );
      }
    }

    for( set<int>::iterator ix=mod_data.xstrips.begin(); ix!=mod_data.xstrips.end(); ++ix ){
      int ixstrip = *ix;
      int besty = mod_data.Bestmatch_xstrips[ixstrip];
      if( mod_data.BestCor_xstrips[ixstrip] > stripcorthreshold ){ //match with correlation coefficient above threshold:
	mod_data.xstrips_filtered.insert( ixstrip );
	mod_data.ystrips_filtered.insert( besty );
      }
    }
  }
}


int find_clusters_by_module( moduledata_t mod_data, clusterdata_t &clust_data ){
  //Don't do 1D clustering any more; instead, now that we have filtered strips by XY matching correlation coefficient, the plan is to
  // group together contiguous 2D pixels: start with max correlation coefficient; build out from there:

  //cout << "Starting strip filtration by correlation coefficient:" << endl;
  filter_strips_by_module( mod_data, stripcorthreshold );
  //  cout << "Filtered strips successfully" << endl;
  
  int module = mod_data.modindex;
  int layer = mod_data.layerindex;
  
  int nclust=0;
  double cormax=-1.1;

  bool foundclust=true;

  //
  map<int,bool> pixelused;
  map<int,bool> stripxused,stripyused;
  //for( set<int>::iterator ix=mod_data.xstrips_filtered.begin(); ix!=mod_data.xstrips_filtered.end(); ++ix ){
  //Let's initialize this flag for ALL fired strips, regardless of whether they passed a correlation threshold:
  for( set<int>::iterator ix=mod_data.xstrips.begin(); ix != mod_data.xstrips.end(); ++ix ){
    int ixstrip = *ix;
    stripxused[ixstrip] = false;
  }
  for( set<int>::iterator iy=mod_data.ystrips.begin(); iy != mod_data.ystrips.end(); ++iy ){
    int iystrip = *iy;
    stripyused[iystrip] = false;
  }

  int returnval = 0;
  
  while( foundclust ){

    foundclust = false;
    
    int nhitclust=0;

    //choose starting pixel by largest ADC value, largest ADC sum, or largest correlation coefficient?
    //initially, let's go with largest correlation coefficient:

    double maxADCXY=-1.0;
    int ixmax=-1,iymax=-1,pixelmax=-1;

    double maxcor=-1.1;
    
    //Find the xy pair with largest correlation coefficient, try to build a cluster around it:

    int ncombos = mod_data.xstrips_filtered.size()*mod_data.ystrips_filtered.size();

    // We probably shouldn't hard-code this, but make it user-configurable: in high-rate conditions we may want to allow much larger numbers of possible strip combinations:
    if( ncombos > 100000 ){
      //     cout << "too many strip combos this module, giving up, module..." << module << endl;
      returnval = -1;
      break;
    }

    //The "filtered" lists of x (y) strips are those whose "best" correlation coefficient with any y (x) strips exceeds "stripcorthreshold"
    
    for( set<int>::iterator ix=mod_data.xstrips_filtered.begin(); ix != mod_data.xstrips_filtered.end(); ++ix ){
      //for( set<int>::iterator iy=mod_data.ystrips_filtered.begin(); iy != mod_data.ystrips_filtered.end(); ++iy ){
      int ixstrip=*ix;
      for( set<int>::iterator iy=mod_data.ystrips_filtered.begin(); iy != mod_data.ystrips_filtered.end(); ++iy ){
	int iystrip = *iy;
	//	int pixel = ixstrip+iystrip*nstripsy;
	
	//	double ADCXY = mod_data.ADCmax_xstrips[ixstrip] * mod_data.ADCmax_ystrips[iystrip]; 

	//Question: do we really need to compute these again? I think maybe not, but to be on the safe side let's do it anyway:
	
	double sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;

	for( int isamp=0; isamp<6; isamp++ ){
	  double ADCx = mod_data.ADCsamp_xstrips[ixstrip][isamp];
	  double ADCy = mod_data.ADCsamp_ystrips[iystrip][isamp];
	  
	  sumx += ADCx;
	  sumy += ADCy;
	  sumx2 += pow(ADCx,2);
	  sumy2 += pow(ADCy,2);
	  sumxy += ADCx*ADCy;
	}
	
	double mux = sumx/6.0;
	double muy = sumy/6.0;
	double varx = sumx2/6.0-pow(mux,2);
	double vary = sumy2/6.0-pow(muy,2);
	double sigx = sqrt(varx);
	double sigy = sqrt(vary);
	
	double ccor = ( sumxy - 6.0*mux*muy )/( 6.0*sigx*sigy );

	double ADCXY = sqrt( mod_data.ADCsum_xstrips[ixstrip]*mod_data.ADCsum_ystrips[iystrip] );

	// double txcorr = mod_data.Tmean_xstrips[ixstrip]+walkcor_mean_func->Eval( mod_data.ADCsum_xstrips[ixstrip] );
      // 	double tycorr = mod_data.Tmean_ystrips[iystrip]+walkcor_mean_func->Eval( mod_data.ADCsum_ystrips[iystrip] );

	
      // //      double tavgcorr = 0.5*(txcorr+tycorr);  
      // 	double tsigmaxcorr = walkcor_sigma_func->Eval( mod_data.ADCsum_xstrips[ixstrip] );
      // 	double tsigmaycorr = walkcor_sigma_func->Eval( mod_data.ADCsum_ystrips[iystrip] );

	double dtcorr = mod_data.Tmean_xstrips_walkcor[ixstrip]-mod_data.Tmean_ystrips_walkcor[iystrip];
	double sigdt = sqrt(pow(mod_data.Tsigma_xstrips_walkcor[ixstrip],2)+pow(mod_data.Tsigma_ystrips_walkcor[iystrip],2));
	
	//only seed a cluster from a local maximum if the correlation coefficient is reasonable:
	if( !stripxused[ixstrip] && !stripyused[iystrip] && ccor >= stripcorthreshold &&
	    fabs( dtcorr ) <= tstripcut_nsigma*sigdt ){
	  if( ixmax < 0 || ccor > maxcor ){
	  //if( ixmax < 0 || ADCXY > maxADCXY ){
	    maxcor = ccor;
	    maxADCXY = ADCXY;
	    ixmax = ixstrip;
	    iymax = iystrip;
	  }
	}
      }
    }

    //    if( ixmax >= 0 && maxADCXY >= thresh_stripsum ){ //local maximum established:
    if( ixmax >= 0 ){
      //      cout << "found local maximum, (ixmax, iymax, maxcor)=(" << ixmax << ", " << iymax << ", " << maxcor << ")" << endl;
      //     cout << "found local maximum, (ixmax,iymax,pixelmax)=(" << ixmax << ", " << iymax << ", " << pixelmax << ")" << endl;
      //    if( ixmax >= 0 && maxcor >= maxstripcorthreshold ){ //start a new cluster:
      
      nhitclust = 1;

      //nclust++;

      double sum_tcorr = 0.0;
      double sumw_tcorr = 0.0;
      
      double txcorr = mod_data.Tmean_xstrips_walkcor[ixmax];
      double tycorr = mod_data.Tmean_ystrips_walkcor[iymax];


      //      double tavgcorr = 0.5*(txcorr+tycorr);  
      double tsigmaxcorr = mod_data.Tsigma_xstrips_walkcor[ixmax];
      double tsigmaycorr = mod_data.Tsigma_ystrips_walkcor[iymax];
      
      sum_tcorr += txcorr*pow(tsigmaxcorr,-2) + tycorr*pow(tsigmaycorr,-2);
      sumw_tcorr += pow(tsigmaxcorr,-2)+pow(tsigmaycorr,-2);

      double tcorr_mean = sum_tcorr/sumw_tcorr;
      
      //We'll use these later to compute cluster moments:
      set<int> xstriplistclust; 
      set<int> ystriplistclust;

      //x strips by "pixel index"
      vector<int> ixhit;
      //y strips by "pixel index"
      vector<int> iyhit;
      
      //int iymax = mod_data.Bestmatch_xstrips[ixmax];
      //pixelused[pixelmax] = true;
      stripxused[ixmax] = true;
      stripyused[iymax] = true;
      
      xstriplistclust.insert( ixmax );
      ystriplistclust.insert( iymax );
      
      int ixpixel=ixmax, iypixel=iymax, ihit=0;

      ixhit.push_back(ixpixel);
      iyhit.push_back(iypixel);
      
      while( ihit < nhitclust ){
	ixpixel = ixhit[ihit];
	iypixel = iyhit[ihit];
	
	int ixleft = ixpixel-1;
	int ixright = ixpixel+1;
	int iyleft = iypixel-1;
	int iyright = iypixel+1;

	//check x left:
	//if( mod_data.xstrips_filtered.find( ixleft ) != mod_data.xstrips_filtered.end() ){ //x strip to the left is in the list of filtered strips:
	if( mod_data.xstrips.find( ixleft ) != mod_data.xstrips.end() && xstriplistclust.size() < maxnstripXpercluster ){ //the strip to the immediate left fired:
	  int bestymatch = mod_data.Bestmatch_xstrips[ixleft];

	  double txcorr = mod_data.Tmean_xstrips_walkcor[ixleft];
	  double sigtxcorr = mod_data.Tsigma_xstrips_walkcor[ixleft];
	  
	  
	  if( !stripxused[ixleft] && abs( bestymatch - iypixel ) < 3 &&
	      fabs( txcorr - tcorr_mean ) <= tstripcut_nsigma*sigtxcorr ){
	    xstriplistclust.insert( ixleft );
	    ixhit.push_back( ixleft );
	    iyhit.push_back( iypixel );
	    //pixelused[ipixel_xleft] = true;
	    stripxused[ixleft] = true;

	    sum_tcorr += txcorr * pow( sigtxcorr, -2 );
	    sumw_tcorr += pow( sigtxcorr, -2 );

	    tcorr_mean = sum_tcorr/sumw_tcorr;
	    
	    nhitclust++;
	  }
	}

	//check x right:
	//if( mod_data.xstrips_filtered.find( ixright ) != mod_data.xstrips_filtered.end() ){ //x strip to the right is in the list of filtered strips:
	if( mod_data.xstrips.find( ixright ) != mod_data.xstrips.end() && xstriplistclust.size() < maxnstripXpercluster ){
	  
	  int bestymatch = mod_data.Bestmatch_xstrips[ixright];

	  double txcorr = mod_data.Tmean_xstrips_walkcor[ixright];   
	  double sigtxcorr = mod_data.Tsigma_xstrips_walkcor[ixright];
	  
	  if( !stripxused[ixright] && abs( bestymatch - iypixel ) < 3 &&
	      fabs( txcorr - tcorr_mean ) <= tstripcut_nsigma*sigtxcorr ){
	    xstriplistclust.insert( ixright );
	    ixhit.push_back( ixright );
	    iyhit.push_back( iypixel );
	    stripxused[ixright] = true;

	    sum_tcorr += txcorr * pow( sigtxcorr, -2 );
	    sumw_tcorr += pow( sigtxcorr, -2 );

	    tcorr_mean = sum_tcorr/sumw_tcorr;
	    
	    nhitclust++;
	  }
	}
      
	//check y left:
	//if( mod_data.ystrips_filtered.find( iyleft ) != mod_data.ystrips_filtered.end() ){
	if( mod_data.ystrips.find( iyleft ) != mod_data.ystrips.end() && ystriplistclust.size() < maxnstripYpercluster ){

	  int bestxmatch = mod_data.Bestmatch_ystrips[iyleft];

	  double tycorr = mod_data.Tmean_ystrips_walkcor[iyleft];
	  double sigtycorr = mod_data.Tsigma_ystrips_walkcor[iyleft];
	  
	  if( !stripyused[iyleft] && abs( bestxmatch - ixpixel ) < 3 &&
	      fabs( tycorr - tcorr_mean ) <= tstripcut_nsigma*sigtycorr ){
	    ystriplistclust.insert( iyleft );
	    ixhit.push_back( ixpixel );
	    iyhit.push_back( iyleft );
	    stripyused[iyleft] = true;

	    sum_tcorr += tycorr * pow( sigtycorr, -2 );
	    sumw_tcorr += pow( sigtycorr, -2 );

	    tcorr_mean = sum_tcorr/sumw_tcorr;
	    
	    nhitclust++;
	  }
	}

	//check y right:
	//	if( mod_data.ystrips_filtered.find( iyright ) != mod_data.ystrips_filtered.end() ){
	if( mod_data.ystrips.find( iyright ) != mod_data.ystrips.end() && ystriplistclust.size() < maxnstripYpercluster ){
	  int bestxmatch = mod_data.Bestmatch_ystrips[iyright];

	  double tycorr = mod_data.Tmean_ystrips_walkcor[iyright];
	  double sigtycorr = mod_data.Tsigma_ystrips_walkcor[iyright];
	  
	  if( !stripyused[iyright] && abs( bestxmatch - ixpixel ) < 3 &&
	      fabs( tycorr - tcorr_mean ) <= tstripcut_nsigma*sigtycorr ){
	    ystriplistclust.insert( iyright );
	    ixhit.push_back( ixpixel );
	    iyhit.push_back( iyright );
	    stripyused[iyright] = true;

	    sum_tcorr += tycorr * pow(sigtycorr,-2);
	    sumw_tcorr += pow( sigtycorr, -2 );

	    tcorr_mean = sum_tcorr/sumw_tcorr;
	    
	    nhitclust++;
	  }
	}
      
	ihit++;
      }

      //      cout << "finished adding strips to cluster, (nclust, nstripx,nstripy)=(" << nclust << ", " << xstriplistclust.size() << ", "
      //	   << ystriplistclust.size() << ")" << endl;
      
      //here we have a cluster: compute cluster properties and add it to the arrays:

      //we also want to compute a 2D cluster six-sample correlation coefficient :
      
      double ADCsumy=0.0;
      double ADCsumx=0.0, xsum=0.0, ysum=0.0, xsum2=0.0, ysum2=0.0, txsum=0.0, txsum2=0.0,tysum=0.0,tysum2=0.0;      

      int ixlo=nstripsx+1,ixhi=-1,iylo=nstripsy+1,iyhi=-1;

      //cluster ADC x and ADC y sums for individual time samples:
      double sumxsamp[6],sumysamp[6];

      for( int isamp=0; isamp<6; isamp++ ){
	sumxsamp[isamp] = 0.0;
	sumysamp[isamp] = 0.0;
      }

      double maxADCx=0.0, maxADCy=0.0;
      int ixmaxadc=-1,iymaxadc=-1;
      double xstripmaxadc=-1e9, ystripmaxadc=-1e9;
      
      //cout << "nclust, nhitclust, ixmax = " << nclust << ", " << nhitclust << ", " <<  ixmax << endl;
      
      for( set<int>::iterator ix=xstriplistclust.begin(); ix != xstriplistclust.end(); ++ix ){
	int ixstrip = *ix;

	//	cout << "x strip " << ixstrip << endl;

	//	cout << "found in list of strips in module = " << (mod_data.xstrips.find(ixstrip) != mod_data.xstrips.end()) << endl;
	
	ixlo = (ixstrip < ixlo ) ? ixstrip : ixlo;
	ixhi = (ixstrip > ixhi ) ? ixstrip : ixhi;


	double ADCtemp = mod_data.ADCsum_xstrips[ixstrip];
	ADCsumx += ADCtemp;

	//Don't hard-code the strip pitch any more:
	double xlocal = (ixstrip + 0.5 - 0.5*mod_nstripsu[module])*mod_ustrip_pitch[module];
	double tstrip = mod_data.Tmean_xstrips[ixstrip];

	if( ixmaxadc < 0 || ADCtemp > maxADCx ){
	  maxADCx = ADCtemp;
	  ixmaxadc = ixstrip;
	  xstripmaxadc = xlocal;
	}
	
	xsum += ADCtemp * xlocal;
	xsum2 += ADCtemp * pow(xlocal,2);

	txsum += ADCtemp * tstrip;
	txsum2 += ADCtemp * pow(tstrip,2);

	for( int isamp=0; isamp<6; isamp++ ){
	  sumxsamp[isamp] += mod_data.ADCsamp_xstrips[ixstrip][isamp];
	}
      }

      for( set<int>::iterator iy=ystriplistclust.begin(); iy != ystriplistclust.end(); ++iy ){
	int iystrip = *iy;

	//	cout << "y strip " << iystrip << endl;
	
	iylo = (iystrip < iylo ) ? iystrip : iylo;
	iyhi = (iystrip > iyhi ) ? iystrip : iyhi;
	
	double ADCtemp = mod_data.ADCsum_ystrips[iystrip];
	ADCsumy += ADCtemp;
	//don't hard-code the strip pitch any more:
	double ylocal = (iystrip + 0.5 - 0.5*mod_nstripsv[module])*mod_vstrip_pitch[module];
	double tstrip = mod_data.Tmean_ystrips[iystrip];

	if( iymaxadc < 0 || ADCtemp > maxADCy ){
	  maxADCy = ADCtemp;
	  iymaxadc = iystrip;
	  ystripmaxadc = ylocal;
	}
	
	ysum += ADCtemp * ylocal;
	ysum2 += ADCtemp * pow(ylocal,2);

	tysum += ADCtemp * tstrip;
	tysum2 += ADCtemp * pow(tstrip,2);
	for( int isamp=0; isamp<6; isamp++ ){
	  sumysamp[isamp] += mod_data.ADCsamp_ystrips[iystrip][isamp];
	}
      }

      vector<double> xADCsamples,yADCsamples;
      
      double csumxsamp =0.0, csumysamp=0.0, csumx2samp=0.0, csumy2samp=0.0, csumxysamp=0.0;
      for( int isamp=0; isamp<6; isamp++ ){
	csumxsamp += sumxsamp[isamp];
	csumysamp += sumysamp[isamp];
	csumx2samp += pow(sumxsamp[isamp],2);
	csumy2samp += pow(sumysamp[isamp],2);
	csumxysamp += sumxsamp[isamp]*sumysamp[isamp];
	xADCsamples.push_back( sumxsamp[isamp] );
	yADCsamples.push_back( sumysamp[isamp] );
      }

      double mux = csumxsamp/6.0;
      double muy = csumysamp/6.0;
      double varx = csumx2samp/6.0-pow(mux,2);
      double vary = csumy2samp/6.0-pow(muy,2);
      double sigx = sqrt(varx);
      double sigy = sqrt(vary);

      double ccor = (csumxysamp - 6.0*mux*muy)/(6.0*sigx*sigy);

      //      if( ccor >= clustcorthreshold && ADCsumx+ADCsumy > thresh_clustersum*2.0 ){
      
      clust_data.nstripx.push_back( xstriplistclust.size() );
      clust_data.ixstriplo.push_back( ixlo );
      clust_data.ixstriphi.push_back( ixhi );
      clust_data.ixstripmax.push_back( ixmaxadc );
      clust_data.xmean.push_back( xsum/ADCsumx );
      clust_data.xsigma.push_back( xsum2/ADCsumx - pow( xsum/ADCsumx,2 ) );
      clust_data.totalchargex.push_back( ADCsumx );
      clust_data.txmean.push_back( txsum/ADCsumx );
      clust_data.txsigma.push_back( txsum2/ADCsumx - pow( txsum/ADCsumx,2 ) );

      clust_data.nstripy.push_back( ystriplistclust.size() );
      clust_data.iystriplo.push_back( iylo );
      clust_data.iystriphi.push_back( iyhi );
      clust_data.iystripmax.push_back( iymaxadc );
      clust_data.ymean.push_back( ysum/ADCsumy );
      clust_data.ysigma.push_back( ysum2/ADCsumy - pow( ysum/ADCsumy,2 ) );
      clust_data.totalchargey.push_back( ADCsumy );
      clust_data.tymean.push_back( tysum/ADCsumy );
      clust_data.tysigma.push_back( tysum2/ADCsumy - pow( tysum/ADCsumy,2 ) );

      clust_data.itrack_clust2D.push_back( -1 );
      clust_data.ixclust2D.push_back( nclust );
      clust_data.iyclust2D.push_back( nclust );
      clust_data.nstripx2D.push_back( xstriplistclust.size() );
      clust_data.nstripy2D.push_back( ystriplistclust.size() );
      clust_data.xclust2D.push_back(clust_data.xmean[nclust]);
      clust_data.yclust2D.push_back(clust_data.ymean[nclust]);
      clust_data.Eclust2D.push_back(0.5*(clust_data.totalchargex[nclust]+clust_data.totalchargey[nclust]));
      clust_data.dEclust2D.push_back( clust_data.totalchargex[nclust]-clust_data.totalchargey[nclust] );
      clust_data.tclust2D.push_back( 0.5*(clust_data.txmean[nclust]+clust_data.tymean[nclust]));
      clust_data.tclust2Dwalkcorr.push_back( tcorr_mean );
      clust_data.dtclust2D.push_back( clust_data.txmean[nclust]-clust_data.tymean[nclust] );
      clust_data.dtclust2Dwalkcorr.push_back( pow( sumw_tcorr, -0.5 ) );
      //clust_data.CorrCoeff2D.push_back( mod_data.BestCor_xstrips[ixmax] );
      clust_data.CorrCoeff2D.push_back( ccor );
      clust_data.CorrCoeffMaxStrips.push_back( maxcor );
      clust_data.keepclust2D.push_back( true );
      clust_data.ADCsamp_xclust.push_back( xADCsamples );
      clust_data.ADCsamp_yclust.push_back( yADCsamples );
      

      clust_data.xmom2D.push_back( (clust_data.xmean[nclust] - xstripmaxadc)/mod_ustrip_pitch[module] );
      clust_data.ymom2D.push_back( (clust_data.ymean[nclust] - ystripmaxadc)/mod_vstrip_pitch[module] );

      double xclust_corr = clust_data.xmean[nclust];
      double yclust_corr = clust_data.ymean[nclust];
      if( usehitmaps ){
	double fracx,fracy;

	//cout << "correcting cluster X coordinates based on hit maps" << endl;

	//hXmom_int->Print();
	
	if( clust_data.nstripx2D[nclust] >= 2 && clust_data.nstripx2D[nclust] <= hXmom_int->GetNbinsX() ){

	  int bin = hXmom_int->FindBin( clust_data.nstripx2D[nclust], clust_data.xmom2D[nclust] );
	  int binx,biny,binz;
	  //cout << "bin = " << bin << endl;
	  
	  hXmom_int->GetBinXYZ( bin, binx, biny, binz );
	  //Each bin content is the integral from 1 up to that bin, inclusive. Therefore, to interpolate within the bin, we should linearly interpolate between the bin itself and the previous bin:
	  //if( biny > 1 ){
	  int binlo = hXmom_int->GetBin( binx, biny-1 );

	  //cout << "bin, binlo " << bin << ", " << binlo << endl;
	  
	  double binfrac = (clust_data.xmom2D[nclust] - hXmom_int->GetYaxis()->GetBinLowEdge( biny ) )/hXmom_int->GetYaxis()->GetBinWidth(biny);
	  //linearly interpolate within bin: maybe later we'd like to use some kind of spline interpolation but this seems simplest for now:

	  //cout << "binfrac = " << binfrac << endl;
	  fracx = hXmom_int->GetBinContent( binlo ) * binfrac + hXmom_int->GetBinContent( bin ) * (1.0 - binfrac);
	  xclust_corr = xstripmaxadc + mod_ustrip_pitch[module] * ( fracx - 0.5 );
	}
	  
	//	cout << "correcting cluster Y coordinates based on hit maps" << endl;

	if( clust_data.nstripy2D[nclust] >= 2 && clust_data.nstripy2D[nclust] <= hYmom_int->GetNbinsY() ){
	  int bin = hYmom_int->FindBin( clust_data.nstripy2D[nclust], clust_data.ymom2D[nclust] );
	  int binx,biny,binz;
	  hYmom_int->GetBinXYZ( bin, binx, biny, binz );
	  int binlo = hYmom_int->GetBin( binx, biny-1 );
	  double binfrac = (clust_data.ymom2D[nclust] - hYmom_int->GetYaxis()->GetBinLowEdge( biny ) )/hYmom_int->GetYaxis()->GetBinWidth(biny);
	  //linearly interpolate within bin: maybe later we'd like to use some kind of spline interpolation but this seems simplest for now:
	  fracy = hYmom_int->GetBinContent( binlo ) * binfrac + hYmom_int->GetBinContent( bin ) * (1.0 - binfrac );
	  yclust_corr = ystripmaxadc + mod_vstrip_pitch[module] * ( fracy - 0.5 );
	}

	// if( clust_data.nstripx2D[nclust] >= 2 && clust_data.nstripy2D[nclust] >= 2 ){
	//   cout << "(xfrac,xclustcorr-xclust)=(" << fracx << ", " << xclust_corr - clust_data.xmean[nclust] << ")" << endl;
	//   cout << "(yfrac,yclustcorr-yclust)=(" << fracy << ", " << yclust_corr - clust_data.ymean[nclust] << ")" << endl;
	// }
      }

     
      
      clust_data.xclust2Dcorr.push_back( xclust_corr );
      clust_data.yclust2Dcorr.push_back( yclust_corr );
      
      //Now, global coordinate computation is going to be modified to handle (almost) arbitrary strip orientations:
      // U strip angle is assumed to be measured FROM the X axis toward the Y axis:
      // therefore, the "X" strips are actually measuring "U = X cos alphau + Y sin alphau = X*Pxu + Y*Pyu"
      //  AND       the "Y" strips are actually measuring "V = X cos alphav + Y sin alphav = X*Pxv + Y*Pyv"
      // To get X and Y from U and V we need the inverse transformation:
      // double Utemp = clust_data.xmean[nclust];
      // double Vtemp = clust_data.ymean[nclust];

      double Utemp = clust_data.xclust2Dcorr[nclust];
      double Vtemp = clust_data.yclust2Dcorr[nclust];
      
      double det = mod_Pxu[module]*mod_Pyv[module] - mod_Pyu[module]*mod_Pxv[module]; //cos( alphau) * sin(alphav) - sin(alphau)*cos(alphav) = 1 for alphau = 0, alphav = 90
      double Xtemp = (mod_Pyv[module]*Utemp - mod_Pyu[module]*Vtemp)/det; //(sin(alphav)*U - sin(alphau)*V)/det = U = X for alphau = 0, alphav = 90
      double Ytemp = (mod_Pxu[module]*Vtemp - mod_Pxv[module]*Utemp)/det; //(cos(alphau)*V - cos(alphav)*U)/det = V = Y for alphau = 0, alphav = 90

      //      cout << "(module, U, V, X, Y)=(" << module << ", " << Utemp << ", " << Vtemp << ", " << Xtemp << ", " << Ytemp << ")"
      //   << endl;
      
      //compute global hit coordinates ONCE:
      //TVector3 hitpos_local( clust_data.xmean[nclust], clust_data.ymean[nclust], 0.0 );
      TVector3 hitpos_local( Xtemp, Ytemp, 0.0 );
     

      TVector3 modcenter_global(mod_x0[module], mod_y0[module], mod_z0[module] );
      TVector3 hitpos_global = mod_Rot[module] * hitpos_local + modcenter_global;

      clust_data.xglobal2D.push_back( hitpos_global.X() );
      clust_data.yglobal2D.push_back( hitpos_global.Y() );
      clust_data.zglobal2D.push_back( hitpos_global.Z() );

      
      
      nclust++;
				    
      clust_data.nclustx = nclust;
      clust_data.nclusty = nclust;
      clust_data.nclust2D = nclust;

      foundclust = true;
      
    }
  }

  clust_data.modindex=module;
  clust_data.layerindex=layer;

  // clust_data.nhitx1D = 0;
  //clust_data.nhity1D = 0;
  //clust_data.nhit2D = 0;

  clust_data.nclustx = nclust;
  clust_data.nclusty = nclust;
  clust_data.nclust2D = nclust;
  
  //cout << "finished cluster finding, nclust = " << nclust << endl;

  return returnval;
}


void prune_clusters( clusterdata_t &clusttemp ){
  //Here we carry out a series of loops over all the clusters in this module. Within each loop, we check one (or perhaps more)
  //criteria, including cluster sum, max. strip correlation, total cluster correlation, ADC asymmetry, time difference, number of strips, etc:
  //within each loop, if at least one cluster passes the criterion, we reject all clusters that don't pass the criterion, guaranteeing that we will always keep at least one
  //good cluster per module:
  //The "keep" flag for each cluster is always initialized to true:
  //order of criteria evaluation will have some effect on selection.
  //At each stage of the pruning algorithm, 
  
  int ngood = 0;

  //sqrt(ADCX*ADCY)>=threshold_clustersum
  for( int pass=0; pass<2; pass++ ){
    for( int iclust=0; iclust<clusttemp.nclust2D; iclust++ ){
      double ADCXY = sqrt(clusttemp.totalchargex[clusttemp.ixclust2D[iclust]]*clusttemp.totalchargey[clusttemp.iyclust2D[iclust]] );
      if( ADCXY >= thresh_clustersum && clusttemp.keepclust2D[iclust] && pass == 0 ){
	ngood++;
      }

      if( pass == 1 && ngood > 0 && ADCXY < thresh_clustersum ){
	clusttemp.keepclust2D[iclust] = false;
      }
    }
  }

  ngood = 0;

  //Time correlation:
  for( int pass=0; pass<2; pass++ ){
    for( int iclust=0; iclust<clusttemp.nclust2D; iclust++ ){
      if( clusttemp.keepclust2D[iclust] && fabs( clusttemp.dtclust2D[iclust] ) <= cluster2Dmatch_tcut && pass == 0 ){ //cluster passed all previous cuts and passes current cut; 
	ngood++;
	//keepclust2D[iclust] = true;
      }

      if( pass == 1 && ngood > 0 && fabs( clusttemp.dtclust2D[iclust] ) > cluster2Dmatch_tcut ){
	clusttemp.keepclust2D[iclust] = false;
      }
    }
  }

  ngood = 0;
  //ADC asymmetry:
  for( int pass=0; pass<2; pass++ ){
    for( int iclust=0; iclust<clusttemp.nclust2D; iclust++ ){
      double ADCasym = clusttemp.dEclust2D[iclust]/(2.0*clusttemp.Eclust2D[iclust]);
      if( clusttemp.keepclust2D[iclust] && fabs( ADCasym ) <= cluster2Dmatch_asymcut && pass == 0 ){
	ngood++;
	//clusttemp.keepclust2D[iclust] = true;
      }
      
      if( pass == 1 && ngood > 0 && fabs( ADCasym ) > cluster2Dmatch_asymcut ){
	clusttemp.keepclust2D[iclust] = false;
      }
    }
  }

  ngood = 0;
  //Max. strip correlation coefficient:
  for( int pass=0; pass<2; pass++ ){
    for( int iclust=0; iclust<clusttemp.nclust2D; iclust++ ){
      if( clusttemp.keepclust2D[iclust] && clusttemp.CorrCoeffMaxStrips[iclust] >= maxstripcorthreshold && pass == 0 ){
	ngood++;
      }

      if( pass == 1 && ngood > 0 && clusttemp.CorrCoeffMaxStrips[iclust] < maxstripcorthreshold ){
	clusttemp.keepclust2D[iclust] = false;
      }
    }
  }
  
  ngood = 0;
  //Cluster correlation coefficient:
  for( int pass=0; pass<2; pass++ ){
    for( int iclust=0; iclust<clusttemp.nclust2D; iclust++ ){
      if( clusttemp.keepclust2D[iclust] && clusttemp.CorrCoeff2D[iclust] >= clustcorthreshold && pass == 0 ){
	ngood++;
      }

      if( pass == 1 && ngood > 0 && clusttemp.CorrCoeff2D[iclust] < clustcorthreshold ){
	clusttemp.keepclust2D[iclust] = false;
      }
    }
  }
}

void find_tracks( map<int,clusterdata_t> mod_clusters, trackdata_t &trackdata ){
  //only attempt tracking if we have at least three layers with at least one 2D matched hit passing the XY and T correlation cuts:

  //  int nhitsrequired = nlayers;
  
  set<int> layers_2Dmatch;
  map<int,int> N2Dhits_layer; //number of 2D matched hits per layer
  map<int,vector<int> > modindexhit2D; //module index of hits mapped by layer
  map<int,vector<int> > clustindexhit2D; //index position of hits within 2D cluster array
  map<int,vector<bool> > hitused2D; //flag to indicate hit is used in a track:

  //Here we are populating the hit arrays, mapped by layer, that won't change throughout the track-finding iterations:
  for( map<int,clusterdata_t>::iterator imod=mod_clusters.begin(); imod != mod_clusters.end(); ++imod ){
    int module = imod->first;
    //clusterdata_t clusttemp = mod_clusters[module];
    int layer = mod_layer[module];

    prune_clusters( mod_clusters[module] );
    
    if( mod_clusters[module].nclust2D > 0 ){
      for( int iclust=0; iclust<mod_clusters[module].nclust2D; iclust++ ){
	// double ADCasym = mod_clusters[module].dEclust2D[iclust]/(2.*mod_clusters[module].Eclust2D[iclust]);
	// double Tdiff   = mod_clusters[module].dtclust2D[iclust];

	// double corrcoeff = mod_clusters[module].CorrCoeff2D[iclust];

	
	// if( fabs( ADCasym ) < cluster2Dmatch_asymcut && fabs( Tdiff ) < cluster2Dmatch_tcut &&
	//     fabs(mod_clusters[module].dEclust2D[iclust])<1000.0 && mod_clusters[module].nstripx2D[iclust]<=5 &&
	//     mod_clusters[module].nstripy2D[iclust]<=5 ){ // GOOD match:
	// if( fabs( ADCasym ) < cluster2Dmatch_asymcut && fabs( Tdiff ) < cluster2Dmatch_tcut &&
	//     corrcoeff > clustcorthreshold && mod_clusters[module].Eclust2D[iclust] > thresh_clustersum ){

	if( mod_clusters[module].keepclust2D[iclust] ){
	
	  layers_2Dmatch.insert( layer );
	  modindexhit2D[layer].push_back( module );
	  clustindexhit2D[layer].push_back( iclust );
	  hitused2D[layer].push_back( false );
	  
	  N2Dhits_layer[layer] = modindexhit2D[layer].size();
	}  
	
      }
    }
  }

  if( layers_2Dmatch.size() >= 3 ){
    //trackdata_t trackdatatemp;

    trackdata.ntracks = 0;
    
    bool foundtrack=true;
    while( foundtrack ){ //consider all possible combinations of one (2D matched) hit per layer:
      // This "brute force" approach will work reasonably well for cosmic data, but we'll need to be smarter when dealing with real data
      //under high-rate conditions
      //set<int> layers_2Dmatch; //list of layers with at least one 2D-matched hit:
      foundtrack = false;

      int nhitsrequired=layers_2Dmatch.size(); //on first iteration, require nhits = total number of layers with unused hits:

      while( nhitsrequired >= 3 ){ //first iteration: determine number of layers with (unused) hits, populate lists of unused hits, count number of combinations,
	//etc:
	
	foundtrack = false;

	//here we populate the lists of free hits by layer: if we found a track on the previous iteration, some hits will have been marked as used
	//if we didn't find a track, then no hits will have been marked as used, but the number of hits required to
	//make a track will have been decremented:
	long ncombosremaining=1;

	map<int,int> Nfreehits_layer; //free hit count mapped by layer
	set<int> layerswithfreehits;  //list of layers with free hits
	map<int,vector<int> > freehitlist_layer; //list of free hits mapped by layer: index in the unchanging arrays defined above:
	map<int,int> freehitcounter; //counter for looping over combos:
	
	for(set<int>::iterator ilay=layers_2Dmatch.begin(); ilay!=layers_2Dmatch.end(); ++ilay ){
	  int layer = *ilay;
	  Nfreehits_layer[layer] = 0;
	  
	  for( int ihit=0; ihit<N2Dhits_layer[layer]; ihit++ ){
	    if( !hitused2D[layer][ihit] ) {
	      Nfreehits_layer[layer]++;
	      freehitlist_layer[layer].push_back( ihit );
	    }
	  }
 
	  if( Nfreehits_layer[layer] > 0 ){
	    ncombosremaining *= Nfreehits_layer[layer];
	    layerswithfreehits.insert(layer);
	    freehitcounter[layer] = 0;
	  }
	}

	//this will get us stuck in an infinite loop if we don't find a track: comment out
	//nhitsrequired = layerswithfreehits.size();

	if( ncombosremaining > maxnhitcombinations || ncombosremaining < 0 ){
	  cout << "too many hit combos, skipping tracking for this event..." << endl;
	  return;
	}

	//	cout << "Starting tracking, ncombosremaining = " << ncombosremaining << endl;
	
	if( layerswithfreehits.size() >= nhitsrequired ){ //If the number of layers with free hits exceeds the number of hits required to make a track, proceed:
	  

	  //number of layers to in total:
	  int nlayerstot = layerswithfreehits.size();

	  //int nlayercombos=pow(2,nlayerstot); //ALL possible ON/OFF combos of the layers with available hits:

	  //vector<vector<int> > layercombos; //all possible combinations of nhitsrequired layers:
	  //This we don't need anymore.
	  
	  // vector<int> listoflayerswithfreehits; //list of all the layers with available hits:
	  // for( set<int>::iterator ilay=layerswithfreehits.begin(); ilay!=layerswithfreehits.end(); ++ilay ){
	  //   listoflayerswithfreehits.push_back( *ilay );
	  // }


	  //Now we only populate the layer combinations once for each nhitsrequired from 3 up to nlayers, after the configuration file is parsed
	  // for( int icombo=0; icombo<pow(2,nlayerstot); icombo++ ){

	  //   vector<bool> onoff(nlayerstot);

	  //   vector<int> layercombo;
	    
	  //   int nlayersoncombo=0;
	  //   for( int ilayer=0; ilayer<nlayerstot; ilayer++ ){ //loop over all layers with available hits:
	  //     int testbit = pow(2,ilayer);
	      
	  //     onoff[ilayer] = ( (testbit & icombo) != 0 ); //bitwise AND of icombo and 2^ilayer non-zero  
	  //     if( onoff[ilayer] ) { //layer ON:
	  // 	nlayersoncombo++; //increment number of layers in this combo:
	  // 	layercombo.push_back( ilayer ); //add layer to the combo:
	  //     }
	  //   }

	  //   if( nlayersoncombo == nhitsrequired ){
	  //     layercombos.push_back( layercombo );
	  //   }
	  // }

	  bool nextcomboexists=true;

	  //vector<int> hitcombo(layerswithfreehits.size());
	  map<int,int> hitcombo; //hits mapped by layer:
	  map<int,bool> ontrack; //flag to indicate whether this layer actually ended up on the track; to handle cases when the number of layers with free hits exceeds the number of layers that end up on the track:
	  map<int,double> xresid_layer; //"U" residual (generalized "X") for test combos:
	  map<int,double> yresid_layer; //"V" residual (generalized "Y") for test combos:
	  bool first = true; //flag to indicate first hit combo to be considered

	  map<int,int> besthitcombo; //map to store "BEST" combo according to Track chi2
	  map<int,bool> onbesttrack; //flag to indicate whether this layer actually ended up on the track; to handle cases when the number of layers with free hits exceeds the number of layers that end up on the track; i.e., when nhitsrequired < n layers with free hits
	  map<int,double> bestxresid; //"U" residual (generalized "X") for best hit combo, mapped by LAYER?
	  map<int,double> bestyresid; //"V" residual (generalized "Y") for best hit combo, mapped by LAYER?
	  
	  //	  map<int,double> xhit,yhit,zhit,xresid,yresid; no longer used

	  double BestTrack[4]; //X,Xp,Y,Yp of best track

	  double bestchi2; //chi2/NDF of best track
	  
	  //cout << "looping over hit combinations, ncombos = " << ncombosremaining << endl;
	  // while( nextcomboexists ){
	    // nextcomboexists=false;

	  int nhits = 0; //number of hits on candidate track
	  int nhitbestcombo = 0;
	  
	  int ngoodcombos=0; //number of plausible track candidates found:
	  
	  for( int icombo=0; icombo<layercombos[nhitsrequired].size(); icombo++ ){ //loop over all possible combinations of nhitsrequired LAYERS with available hits:

	    nextcomboexists = true;
	    
	    set<int> layerstotest; //populate the list of layers to be tested:
	    for( int ihit=0; ihit<nhitsrequired; ihit++ ){

	      int layeri = layercombos[nhitsrequired][icombo][ihit];

	      if( layerswithfreehits.find(layeri) != layerswithfreehits.end() ){
	      
		layerstotest.insert( layeri );
	      //at the beginning of each layer combination, reset the free hit counter to zero for each layer being tested:
		freehitcounter[layeri] = 0;
	      }
	    }

	    if( layerstotest.size() < nhitsrequired ) nextcomboexists = false;

	    first = true;
	    
	    while( nextcomboexists ){ //this will loop over all possible combos:
	      
	      //sums for computation of tracks "on the fly":
	      double sumxhits=0.0,sumyhits=0.0,sumzhits=0.0,sumxzhits=0.0,sumyzhits=0.0,sumz2hits=0.0;
	      
	      double xtrtemp,ytrtemp,ztrtemp,xptrtemp,yptrtemp; //temporary storage for track parameters: 

	      double varx,vary,varxp,varyp,covxxp,covyyp;
	      
	      //start at first layer:
	      //get first possible combination:
	      //int nhits=0;
	      //	      set<int>::iterator nextlayercounter = layerswithfreehits.begin(); //"Next" layer counter for "odometer algorithm:
	      set<int>::iterator nextlayercounter = layerstotest.begin(); //"Next" layer counter for "odometer algorithm:
	      //Reset number of hits found:
	      nhits = 0;

	      for( set<int>::iterator layercounter = layerstotest.begin(); layercounter != layerstotest.end(); ++layercounter ){
		int layer = *layercounter;
		int nextlayer = *nextlayercounter;
		
		ontrack[layer] = false;
	      
		if( layer == nextlayer && !first ){
		  if( freehitcounter[layer]+1 < Nfreehits_layer[layer] ){
		    //increment free hit counter:
		    freehitcounter[layer]++;
		  } else { //reached last hit in current layer; roll back to first hit in this layer and increment hit counter in next layer:
		    // Note that this means on the next iteration of layercounter,
		    // layer == nextlayer will evaluate to true, and we will attempt to increment freehitcounter
		    // for that layer as long as another free hit is available. Meanwhile, since the
		    // free hit counter in the current layer has rolled back to 0, the next time we try to
		    // populate a unique hit combination, starting from the first layer, we will loop over the hits
		    // in the first layer again, having incremented the free hit counter for the next layer.
		    // This process will repeat itself until we reach the last hit
		    // in the last layer, at which point nextlayercounter will evaluate to layerswithfreehits.end,
		    // and the iteration will stop.
		    freehitcounter[layer]=0;
		    ++nextlayercounter; 
		  }
		}
		
		if( nextlayercounter == layerstotest.end() ) nextcomboexists = false;
		
		//set hit combo:
		hitcombo[layer] = freehitlist_layer[layer][freehitcounter[layer]]; //hit combo = index in unchanging hit arrays defined at the beginning:
		int module = modindexhit2D[layer][hitcombo[layer]]; 
		int iclust = clustindexhit2D[layer][hitcombo[layer]];
		
		//Get global hit coordinates:
		double xhittemp = mod_clusters[module].xglobal2D[iclust];
		double yhittemp = mod_clusters[module].yglobal2D[iclust];
		double zhittemp = mod_clusters[module].zglobal2D[iclust]; 
		
		if( nhits == 0 ){ //first hit: use as seed to filter subsequent hits based on fit to straight lines:
		  xtrtemp = xhittemp;
		  ytrtemp = yhittemp;
		  ztrtemp = zhittemp;
		  
		  sumxhits += xhittemp;
		  sumyhits += yhittemp;
		  sumzhits += zhittemp;
		  sumxzhits += xhittemp*zhittemp;
		  sumyzhits += yhittemp*zhittemp;
		  sumz2hits += pow(zhittemp,2);
		  
		  //nhits == 1 after this step:
		  nhits++;
		  ontrack[layer] = true;
		} else if( nhits == 1 ){ //second hit: check within range and update slope:
		  //only if the second hit falls within the plausible range for tracks of the first do we keep:
		  if( fabs( xhittemp - xtrtemp ) < TrackMaxSlopeX * (zhittemp - ztrtemp ) &&
		      fabs( yhittemp - ytrtemp ) < TrackMaxSlopeY * (zhittemp - ztrtemp ) ){
		    //Now compute candidate track slope from line through two hits:
		    xptrtemp = (xhittemp - xtrtemp)/(zhittemp - ztrtemp);
		    yptrtemp = (yhittemp - ytrtemp)/(zhittemp - ztrtemp);
		    xtrtemp = xhittemp - xptrtemp*zhittemp;
		    ytrtemp = yhittemp - yptrtemp*zhittemp;
		    
		    //increment sums:
		    sumxhits += xhittemp;
		    sumyhits += yhittemp;
		    sumzhits += zhittemp;
		    sumxzhits += xhittemp*zhittemp;
		    sumyzhits += yhittemp*zhittemp;
		    sumz2hits += pow(zhittemp,2);
		    
		    //nhits == 2 after this step:
		    nhits++;
		    ontrack[layer] = true;
		  }
		} else { //third and subsequent hits: 
		  //check within range and update track fit:
		  //compute squared distance between this hit and the projection of track from previous hits:
		  double r2track =
		    pow( xhittemp-(xtrtemp + xptrtemp*zhittemp), 2 ) +
		    pow( yhittemp-(ytrtemp + yptrtemp*zhittemp), 2 );

		  varx = pow(sigma_hitpos,2)*sumz2hits/(sumz2hits*nhits - pow(sumzhits,2));
		  varxp = pow(sigma_hitpos,2)*nhits/(sumz2hits*nhits - pow(sumzhits,2));
		  vary = varx;
		  varyp = varxp;
		  covxxp = -pow(sigma_hitpos,2)*sumzhits/(sumz2hits*nhits-pow(sumzhits,2));
		  covyyp = covxxp;

		  //error in x projection = 
		  double dxproj2 = varx + varxp*pow(zhittemp,2) + 2.0*covxxp*zhittemp;
		  double dyproj2 = vary + varyp*pow(zhittemp,2) + 2.0*covyyp*zhittemp;

		  double drproj2 = dxproj2 + dyproj2;
		  
		  //cout << "nhits = " << nhits << ", x projection uncertainty = " << sqrt(dxproj2) << ", y projection uncertainty = " << sqrt(dyproj2)
		  //     << endl;
		  //max. residual corresponding to max. chi2 cut = chi^2 = maxresid^2/sigma^2 --> maxresid = sigma * sqrt(TrackChi2Cut)
		  
		  if( r2track < pow(TrackFindingMaxRadius,2) ){ //This hit falls within plausible range for track (user-configurable radius):
		    
		    //increment sums and number of hits; update track best-fit params:
		    sumxhits += xhittemp;
		    sumyhits += yhittemp;
		    sumzhits += zhittemp;
		    sumxzhits += xhittemp*zhittemp;
		    sumyzhits += yhittemp*zhittemp;
		    sumz2hits += pow(zhittemp,2);
		    
		    nhits++;
		    ontrack[layer] = true;
		    
		    //Update track best-fit params:
		    double denom = (sumz2hits*nhits - pow(sumzhits,2));
		    
		    xptrtemp = (nhits*sumxzhits - sumxhits*sumzhits)/denom;
		    yptrtemp = (nhits*sumyzhits - sumyhits*sumzhits)/denom;
		    xtrtemp = (sumz2hits*sumxhits - sumzhits*sumxzhits)/denom;
		    ytrtemp = (sumz2hits*sumyhits - sumzhits*sumyzhits)/denom;
		  }
		}
	      }
	      
	      if( first ) first = false;
	      
	      //Proceed only if the number of found hits in this combo exceeds the currently required number:
	      if( nhits >= nhitsrequired ){ //compute chi^2 for this hit combo:
		
		//Arrays of global hit coordinates for tracking:
		//no longer useD
		// map<int,double> xtemp,ytemp,ztemp;
		// map<int,double> xresidtemp,yresidtemp;
		
		double chi2 = 0.0;

		//now, we need to modify the chi2 calculation to account for possible different strip orientations:
		//Should we look at U and V residuals instead of "X" and "Y" residuals? YES:
		
		int ndf = 2*nhits-4;
		
		for(set<int>::iterator ilayer=layerswithfreehits.begin(); ilayer!=layerswithfreehits.end(); ++ilayer){
		  int layer = *ilayer;
		  
		  if( ontrack[layer] ){
		    int ihit=hitcombo[layer];
		    
		    int module = modindexhit2D[layer][ihit];
		    int iclust= clustindexhit2D[layer][ihit];
		    
		    //global hit coordinates now stored during cluster finding, no need to recompute:
		    double xhittemp = mod_clusters[module].xglobal2D[iclust];
		    double yhittemp = mod_clusters[module].yglobal2D[iclust];
		    double zhittemp = mod_clusters[module].zglobal2D[iclust]; 

		    //These are the local coordinates of the hit, measured in "Strip" coordinates:
		    double uhittemp = mod_clusters[module].xclust2Dcorr[iclust];
		    double vhittemp = mod_clusters[module].yclust2Dcorr[iclust];

		    //double utracktemp = (xtrtemp+xptrtemp*zhittemp)*mod_Pxu[module] + (ytrtemp+yptrtemp*zhittemp)*mod_Pyu[module];
		    //double vtracktemp = (xtrtemp+xptrtemp*zhittemp)*mod_Pxv[module] + (ytrtemp+yptrtemp*zhittemp)*mod_Pyv[module];

		    TVector3 trackpos_global( xtrtemp+xptrtemp*zhittemp, ytrtemp+yptrtemp*zhittemp, zhittemp );
		    TVector3 modcenter_global( mod_x0[module], mod_y0[module], mod_z0[module] );
		   
		    
		    TVector3 trackpos_local = mod_Rotinv[module]*(trackpos_global - modcenter_global);

		    //Compute the local "u" and "v" coordinates of the track within the module, to allow for the
		    //possibility of different strip orientations than X/Y:
		    double utracktemp = trackpos_local.X()*mod_Pxu[module] + trackpos_local.Y()*mod_Pyu[module];
		    double vtracktemp = trackpos_local.X()*mod_Pxv[module] + trackpos_local.Y()*mod_Pyv[module];
		    
		    // chi2 += pow( (xhittemp - (xtrtemp + xptrtemp*zhittemp))/sigma_hitpos, 2 ) +
		    //   pow( (yhittemp - (ytrtemp + yptrtemp*zhittemp))/sigma_hitpos, 2 );
		    
		    chi2 += pow( (uhittemp-utracktemp)/sigma_hitpos, 2 ) + pow( (vhittemp - vtracktemp)/sigma_hitpos, 2 );

		    xresid_layer[layer] = uhittemp - utracktemp;
		    yresid_layer[layer] = vhittemp - vtracktemp;
		    
		  }
		}
		
		if( ngoodcombos == 0 || chi2/double(ndf) < bestchi2 ){
		  bestchi2 = chi2/double(ndf);
		  besthitcombo = hitcombo;
		  onbesttrack = ontrack;
		  BestTrack[0] = xtrtemp;
		  BestTrack[1] = xptrtemp;
		  BestTrack[2] = ytrtemp;
		  BestTrack[3] = yptrtemp;

		  bestxresid = xresid_layer;
		  bestyresid = yresid_layer;
		  
		  nhitbestcombo = nhits;
		}
		
		ngoodcombos++;
		
	      } //if (nhits > nhitsrequired)
	    } //while (nextcomboexists): end loop over all possible hit combos for current value of nhitsrequired:
	  } //end loop over all possible combinations of nhitsrequired layers:  
	  if( ngoodcombos > 0 && bestchi2 < TrackChi2Cut ){ //add track to global track arrays and mark hits as used:
	    foundtrack = true;
	    trackdata.nhitsontrack.push_back( nhitbestcombo );
	    
	    vector<int> modlisttemp,hitlisttemp;
	    vector<double> residxtemp,residytemp;
	    
	    for( map<int,int>::iterator ihit=besthitcombo.begin(); ihit != besthitcombo.end(); ++ihit ){
	      int layer = ihit->first;
	      
	      if( onbesttrack[layer] ){
		int hit = besthitcombo[layer];
		int module = modindexhit2D[layer][hit];
		int iclust = clustindexhit2D[layer][hit];
		
		hitused2D[layer][hit] = true;
		
		// double xhittemp = mod_clusters[module].xglobal2D[iclust];
		// double yhittemp = mod_clusters[module].yglobal2D[iclust];
		// double zhittemp = mod_clusters[module].zglobal2D[iclust]; 

		// double uhittemp = mod_clusters[module].xclust2D[iclust];
		// double vhittemp = mod_clusters[module].yclust2D[iclust];

		
		// TVector3 trackpos_global (BestTrack[0] + BestTrack[1]*zhittemp, BestTrack[2] + BestTrack[3]*zhittemp, zhittemp );
		// TVector3 modcenter_global( mod_x0[module], mod_y0[module], mod_z0[module] );
		
		
		// TRotation Rmod;
		// Rmod.RotateX( mod_ax[module] );
		// Rmod.RotateY( mod_ay[module] );
		// Rmod.RotateZ( mod_az[module] );
		
		// TVector3 trackpos_local = Rmod.Inverse()*(trackpos_global - modcenter_global);
		
		// double utracktemp = trackpos_local.X()*mod_Pxu[module] + trackpos_local.Y()*mod_Pyu[module];
		// double vtracktemp = trackpos_local.X()*mod_Pxv[module] + trackpos_local.Y()*mod_Pyv[module];
		
		//double utracktemp = (BestTrack[0] + BestTrack[1]*zhittemp)*mod_Pxu[module] + (BestTrack[2] + BestTrack[3])*mod_Pyu[module];
		//double vtracktemp = (BestTrack[0] + BestTrack[1]*zhittemp)*mod_Pxv[module] + (BestTrack[2] + BestTrack[3])*mod_Pyv[module];
		
		modlisttemp.push_back( module );
		hitlisttemp.push_back( iclust );
		// residxtemp.push_back( xhittemp - (BestTrack[0] + BestTrack[1]*zhittemp) );
		// residytemp.push_back( yhittemp - (BestTrack[2] + BestTrack[3]*zhittemp) );
		residxtemp.push_back( bestxresid[layer] );
		residytemp.push_back( bestyresid[layer] );
		
		mod_clusters[module].itrack_clust2D[iclust] = trackdata.ntracks;
	      }
	    }
	    
	    trackdata.modlist_track.push_back( modlisttemp );
	    trackdata.hitlist_track.push_back( hitlisttemp );
	    trackdata.residx_hits.push_back( residxtemp );
	    trackdata.residy_hits.push_back( residytemp );
	    trackdata.eresidx_hits.push_back( residxtemp );
	    trackdata.eresidy_hits.push_back( residytemp );
	    
	    trackdata.Xtrack.push_back( BestTrack[0] );
	    trackdata.Xptrack.push_back( BestTrack[1] );
	    trackdata.Ytrack.push_back( BestTrack[2] );
	    trackdata.Yptrack.push_back( BestTrack[3] );
	    
	    trackdata.Chi2NDFtrack.push_back( bestchi2 );
	    
	    trackdata.ntracks++;
	  }
	  
	} //end check on "layers with free hits >= nhits required"
	  //if( !foundtrack ) break; //prevents getting stuck in infinite loop if we don't find any tracks on first iteration:

	//if we fail to find a track at the current hit requirement; we reduce the number of hits required:
	//If the number of hits required falls below 3, we exit the loop:
	if( !foundtrack ) nhitsrequired--;
	
      } //while( nhitsrequired >= 3 )
	    
	//cout << "nhitsrequired = " << nhitsrequired << endl;
    } //while( foundtrack )
  }
}

void GEM_reconstruct( const char *filename, const char *configfilename, const char *outfilename="temp.root" ){

  //Initialize walk correction parameters:
  double walkcor_mean_params[3] = {walkcor_mean_const, walkcor_mean_ADC0, walkcor_mean_exp};
  double walkcor_sigma_params[3] = {walkcor_sigma_const, walkcor_sigma_ADC0, walkcor_mean_exp};
  walkcor_mean_func->SetParameters(walkcor_mean_params);
  walkcor_sigma_func->SetParameters(walkcor_sigma_params);
  
  gROOT->ProcessLine(".x ~/rootlogon.C");
  gStyle->SetPalette(kRainBow);
  
  gStyle->SetOptStat(0);
  
  TFile *fout = new TFile(outfilename,"RECREATE");

  TTree *Tout = new TTree("Tout","INFN GEM 4-layer cosmic tracks");

  //What branches do we need in our ROOT tree?
  //We need the "local" and "global" hit positions of hits on tracks
  //We need the track parameters
  //We need the track residuals
  //We need the module and layer information.

  //Tree design: Separate track and hit arrays. Since 


  long NMAX = -1; //-1, analyze all events in the file; >= 0 = stop at NMAX
  
  ifstream configfile(configfilename);

  int eventdisplaymode=0; //configure event display mode or "analysis" mode
  
  //configuration parameters: read in default values of x,y,z coordinates of centers of modules, and
  // rotation angles.
  if( configfile ){
    TString currentline;
    
    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig")){
      if( !currentline.BeginsWith("#") ){
	TObjArray *tokens = currentline.Tokenize(" ");

	int ntokens = tokens->GetEntries();

	if( ntokens >= 2 ){
	  TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	  if( skey == "NMAX" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    NMAX = stemp.Atoi();
	  }
	  
	  if( skey == "nlayers" ){
	    TString snlayers = ( (TObjString*) (*tokens)[1] )->GetString();
	    nlayers = snlayers.Atoi();
	  }

	  if( skey == "nmodules" ){
	    TString snmodules = ( (TObjString*) (*tokens)[1] )->GetString();
	    nmodules = snmodules.Atoi();
	  }
	
	  if( skey == "mod_x0" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString smodx = ( (TObjString*) (*tokens)[i] )->GetString();
	      
	      mod_x0[i-1] = smodx.Atof();
	    }
	  }

	  if( skey == "mod_y0" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString smody = ( (TObjString*) (*tokens)[i] )->GetString();
	      
	      mod_y0[i-1] = smody.Atof();
	    }
	  }

	  if( skey == "mod_z0" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString smodz = ( (TObjString*) (*tokens)[i] )->GetString();
	      
	      mod_z0[i-1] = smodz.Atof();
	    }
	  }

	  if( skey == "mod_ax" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString smodax = ( (TObjString*) (*tokens)[i] )->GetString();
	      
	      mod_ax[i-1] = smodax.Atof();
	    }
	  }

	  if( skey == "mod_ay" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString smoday = ( (TObjString*) (*tokens)[i] )->GetString();
	      
	      mod_ay[i-1] = smoday.Atof();
	    }
	  }

	  if( skey == "mod_az" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString smodaz = ( (TObjString*) (*tokens)[i] )->GetString();
	      
	      mod_az[i-1] = smodaz.Atof();
	    }
	  }

	  if( skey == "mod_layer" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString smodlayer = ( (TObjString*) (*tokens)[i] )->GetString();

	      mod_layer[i-1] = smodlayer.Atoi();
	    }  
	  }
	  //Number of strips "U"
	  if( skey == "mod_nstripsu" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      mod_nstripsu[i-1] = stemp.Atoi();
	    }
	  }
	  //Number of strips "V"
	  if( skey == "mod_nstripsv" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      mod_nstripsv[i-1] = stemp.Atoi();
	    }
	  }
	  //Strip pitch "U"
	  if( skey == "mod_ustrip_pitch" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      mod_ustrip_pitch[i-1] = stemp.Atof();
	    }
	  }

	  //Strip pitch "V"
	  if( skey == "mod_vstrip_pitch" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      mod_vstrip_pitch[i-1] = stemp.Atof();
	    }
	  }

	  //These angles will be assumed to refer to the strip orientations! The coordinates they measure will be orthogonal to
	  // their orientations
	  //Angle relative to X axis of "U" strips (assumed to be given in degrees)"
	  if( skey == "mod_uangle" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      mod_uangle[i-1] = stemp.Atof()*PI/180.0;
	      mod_Pxu[i-1] = cos(mod_uangle[i-1]);
	      mod_Pyu[i-1] = sin(mod_uangle[i-1]);
	    }
	  }

	  //Angle relative to X axis of "V" strips (assumed to be given in degrees)"
	  if( skey == "mod_vangle" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      mod_vangle[i-1] = stemp.Atof()*PI/180.0;
	      mod_Pxv[i-1] = cos(mod_vangle[i-1]);
	      mod_Pyv[i-1] = sin(mod_vangle[i-1]);
	    }
	  }

	  //Strip dimensions: width of active area along X
	  if( skey == "mod_Lx" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      mod_Lx[i-1] = stemp.Atof();
	    }
	  }

	  //Strip dimensions: width of active area along Y
	  if( skey == "mod_Ly" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      mod_Ly[i-1] = stemp.Atof();
	    }
	  }

	  if( skey == "mod_uplaneID" && ntokens >= nmodules+1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      mod_uplaneID[i-1] = stemp.Atoi();
	    }
	  }

	  if( skey == "mod_vplaneID" && ntokens >= nmodules+1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      mod_vplaneID[i-1] = stemp.Atoi();
	    }
	  }
	  
	  if( skey == "eventdisplay" && ntokens >= 2 ){
	    TString sevdisplay = ( (TObjString*) (*tokens)[1] )->GetString();

	    eventdisplaymode = sevdisplay.Atoi();
	  }

	  // if( skey == "varyclustersigma" && ntokens >= 2 ){
	  //   TString sclustersigmaflag = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   varyclustersigma = sclustersigmaflag.Atoi();
	  // }

	  // if( skey == "varyclustertau" && ntokens >= 2 ){
	  //   TString sclustertauflag = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   varyclustertau = sclustertauflag.Atoi();
	  // }

	  // if( skey == "maxhitspercluster" && ntokens >= 2 ){ //largest number of hits to be considered as part of same cluster:
	  //   TString smaxhits =  ( (TObjString*) (*tokens)[1] )->GetString();
	  //   maxnhitspercluster = smaxhits.Atoi();
	  // }

	  if( skey == "maxstripsperclusterX" && ntokens >= 2 ){ //largest number of hits to be considered as part of same cluster:
	    TString smaxstrips =  ( (TObjString*) (*tokens)[1] )->GetString();
	    maxnstripXpercluster = smaxstrips.Atoi();
	  }

	  if( skey == "maxstripsperclusterY" && ntokens >= 2 ){ //largest number of hits to be considered as part of same cluster:
	    TString smaxstrips =  ( (TObjString*) (*tokens)[1] )->GetString();
	    maxnstripYpercluster = smaxstrips.Atoi();
	  }
	  
	  // if( skey == "clustersigma" && ntokens >= 2 ){
	  //   TString sclustersigma =  ( (TObjString*) (*tokens)[1] )->GetString();
	  //   clustersigma = sclustersigma.Atof();
	  // }

	  // if( skey == "clustertau" && ntokens >= 2 ){
	  //   TString sclustertau =  ( (TObjString*) (*tokens)[1] )->GetString();
	  //   clustertau = sclustertau.Atof();
	  // }

	  // if( skey == "localmaxthreshold" && ntokens >= 2 ){
	  //   TString sthreshold = ( (TObjString*) (*tokens)[1] )->GetString();

	  //   localmaxthreshold_nsigma=sthreshold.Atof();
	  // }

	 
	  if( skey == "maxcor_threshold" && ntokens >= 2 ){ //threshold on MAX correlation coefficient:
	    TString sthreshold = ( (TObjString*) (*tokens)[1] )->GetString();
	    maxstripcorthreshold = sthreshold.Atof();
	  }
	  
	  if( skey == "stripcor_threshold" && ntokens >= 2 ){
	    TString sthreshold = ( (TObjString*) (*tokens)[1] )->GetString();
	    stripcorthreshold = sthreshold.Atof();
	  }

	  if( skey == "clustcor_threshold" && ntokens >= 2 ){
	    TString sthreshold = ( (TObjString*) (*tokens)[1] )->GetString();
	    clustcorthreshold = sthreshold.Atof();
	  }

	  if( skey == "clust2D_ADCasymcut" && ntokens >= 2 ){
	    TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	    cluster2Dmatch_asymcut = scut.Atof();
	  }

	  if( skey == "clust2D_dTcut" && ntokens >= 2 ){
	    TString scut = ( (TObjString*) (*tokens)[1] )->GetString();

	    cluster2Dmatch_tcut = scut.Atof();
	  }

	  if( skey == "threshold_maxsample" && ntokens >= 2 ){
	    TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	    thresh_maxsample = scut.Atof();
	  }

	  if( skey == "threshold_stripsum" && ntokens >= 2 ){
	    TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	    thresh_stripsum = scut.Atof();
	  }

	  if( skey == "threshold_clustersum" && ntokens >= 2 ){
	    TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	    thresh_clustersum = scut.Atof();
	  }

	  if( skey == "trackchi2cut" && ntokens >= 2 ){
	    TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	    TrackChi2Cut = scut.Atof();
	  }

	  if( skey == "trackmaxradius" && ntokens >= 2 ){
	    TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	    TrackFindingMaxRadius = scut.Atof();
	  }

	  if( skey == "trackmaxslopeX" && ntokens >= 2 ){
	    TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	    TrackMaxSlopeX = scut.Atof();
	  }

	  if( skey == "trackmaxslopeY" && ntokens >= 2 ){
	    TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	    TrackMaxSlopeY = scut.Atof();
	  }
	  
	  // if( skey == "maxADCXYthreshold" && ntokens >= 2 ){
	  //   TString scut = ( (TObjString*) (*tokens)[1] )->GetString();

	  //   maxADCXYthreshold = scut.Atof();
	  // }

	  if( skey == "walkcor_mean_t0" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    walkcor_mean_const = stemp.Atof();
	    walkcor_mean_func->SetParameter( 0, stemp.Atof() );
	  }

	  if( skey == "walkcor_mean_ADC0" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    walkcor_mean_ADC0 = stemp.Atof();
	    walkcor_mean_func->SetParameter( 1, stemp.Atof() );
	  }

	  if( skey == "walkcor_mean_exp" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    walkcor_mean_exp = stemp.Atof();
	    walkcor_mean_func->SetParameter( 2, stemp.Atof() );
	  }

	  if( skey == "walkcor_sigma_t0" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    walkcor_sigma_const = stemp.Atof();
	    walkcor_sigma_func->SetParameter( 0, stemp.Atof() );
	  }

	  if( skey == "walkcor_sigma_ADC0" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    walkcor_sigma_ADC0 = stemp.Atof();
	    walkcor_sigma_func->SetParameter( 1, stemp.Atof() );
	  }

	  if( skey == "walkcor_sigma_exp" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    walkcor_sigma_exp = stemp.Atof();
	    walkcor_sigma_func->SetParameter( 2, stemp.Atof() );
	  }

	  if( skey == "tstripcut_nsigma" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    tstripcut_nsigma = stemp.Atof();
	  }

	  if( skey == "clustmap_fname" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();

	    TFile *ftemp = new TFile(stemp,"READ");

	    TH2D *hXtemp;
	    TH2D *hYtemp;
	    
	    ftemp->GetObject( "hClust2D_Xmom_int_vs_NstripX", hXtemp );
	    ftemp->GetObject( "hClust2D_Ymom_int_vs_NstripY", hYtemp );

	    fout->cd();
	    
	    hXmom_int = new TH2D( *hXtemp );
	    hYmom_int = new TH2D( *hYtemp );	    
	    
	    hXmom_int->Print();
	    hYmom_int->Print();
	    
	    ftemp->Close();
	    ftemp->Delete();

	    cout << "closed hit map file" << endl;

	    hXmom_int->Print();
	    hYmom_int->Print();
	    
	    usehitmaps=true;
	  }
	      
	}
      }
    }
  }

  zavg_layer.resize(nlayers);
  int nmod_layer[nlayers];
  for( int ilayer=0; ilayer<nlayers; ilayer++ ){
    nmod_layer[ilayer] = 0;
    zavg_layer[ilayer] = 0.0;
  }
  for( int imod=0; imod<nmodules; imod++ ){
    int layer = mod_layer[imod];
    nmod_layer[layer]++;
    zavg_layer[layer] += mod_z0[imod];
  }
  for( int ilayer=0; ilayer<nlayers; ilayer++ ){
    zavg_layer[ilayer] /= double(nmod_layer[ilayer]);
    cout << "ilayer, zavg = " << ilayer << ", " << zavg_layer[ilayer] << endl;
  }

  //populate the array of all possible combinations of layers:
  //for( int nhitsrequired=3; nhitsrequired<=nlayers; nhitsrequired++ ){
  for( int icombo=0; icombo<pow(2,nlayers); icombo++ ){

    vector<bool> onoff(nlayers);
    vector<int> layercombo;

    int nlayersoncombo=0;
    for( int ilayer=0; ilayer<nlayers; ilayer++ ){
      int testbit = pow(2,ilayer);
      onoff[ilayer] = ( (testbit & icombo) != 0 ); //bitwise AND of icombo and 2^ilayer non-zero:
      if( onoff[ilayer] ){
	nlayersoncombo++;
	layercombo.push_back( ilayer );
      }
    }

    if( nlayersoncombo >= 3 ){
      layercombos[nlayersoncombo].push_back( layercombo );
    }
  }
  
  //Output ROOT tree structure for alignment code, tracking diagnostics, etc:
  
  //Look at one best track per event (usually the first track found):
  double TrackXp,TrackYp,TrackX,TrackY,TrackChi2NDF;
  int TrackNhits; //Since a track is allowed to be associated with at MOST one hit per layer, we can use a fixed-size array for the hit list:
  //int TrackHitList[nlayers];

  int Ntracks;
  //  int TrackLayers[nlayers];

  int Nclustperlayer[nlayers];
  
  int HitModule[nlayers];
  int HitLayer[nlayers];
  
  double HitXlocal[nlayers];
  double HitYlocal[nlayers];
  double HitXglobal[nlayers];
  double HitYglobal[nlayers];
  double HitZglobal[nlayers];
  double HitXresid[nlayers];
  double HitYresid[nlayers];
  double HitXresidE[nlayers];
  double HitYresidE[nlayers];
  double HitSigX[nlayers]; 
  double HitSigY[nlayers];
  double HitXmom[nlayers];
  double HitYmom[nlayers];
  double HitADCX[nlayers]; 
  double HitADCY[nlayers];
  double HitADCasym[nlayers];
  double HitTmean[nlayers];
  double HitdT[nlayers];
  double HitCorrCoeff[nlayers];
  double StripMaxCorrCoeff[nlayers];
  int HitNstripX[nlayers];
  int HitNstripY[nlayers];

  Tout->Branch("Ntracks",&Ntracks,"Ntracks/I");
  Tout->Branch("TrackXp",&TrackXp,"TrackXp/D");
  Tout->Branch("TrackYp",&TrackYp,"TrackYp/D");
  Tout->Branch("TrackX",&TrackX,"TrackX/D");
  Tout->Branch("TrackY",&TrackY,"TrackY/D");
  Tout->Branch("TrackChi2NDF",&TrackChi2NDF,"TrackChi2NDF/D");
  Tout->Branch("TrackNhits",&TrackNhits,"TrackNhits/I");
  Tout->Branch("Nlayers",&nlayers,"Nlayers/I");
  Tout->Branch("Ncluster",Nclustperlayer,"Ncluster[Nlayers]/I");
  Tout->Branch("HitModule",HitModule,"HitModule[TrackNhits]/I");
  Tout->Branch("HitLayer",HitLayer,"HitLayer[TrackNhits]/I");
  Tout->Branch("HitXlocal",HitXlocal,"HitXlocal[TrackNhits]/D");
  Tout->Branch("HitYlocal",HitYlocal,"HitYlocal[TrackNhits]/D");
  Tout->Branch("HitXglobal",HitXglobal,"HitXglobal[TrackNhits]/D");
  Tout->Branch("HitYglobal",HitYglobal,"HitYglobal[TrackNhits]/D");
  Tout->Branch("HitZglobal",HitZglobal,"HitZglobal[TrackNhits]/D");
  Tout->Branch("HitXresid",HitXresid,"HitXresid[TrackNhits]/D");
  Tout->Branch("HitYresid",HitYresid,"HitYresid[TrackNhits]/D");
  Tout->Branch("HitXresidE",HitXresidE,"HitXresidE[TrackNhits]/D");
  Tout->Branch("HitYresidE",HitYresidE,"HitYresidE[TrackNhits]/D");
  Tout->Branch("HitSigX",HitSigX,"HitSigX[TrackNhits]/D");
  Tout->Branch("HitSigY",HitSigY,"HitSigY[TrackNhits]/D");
  Tout->Branch("HitXmom",HitXmom,"HitXmom[TrackNhits]/D");
  Tout->Branch("HitYmom",HitYmom,"HitYmom[TrackNhits]/D");
  Tout->Branch("HitADCX",HitADCX,"HitADCX[TrackNhits]/D");
  Tout->Branch("HitADCY",HitADCY,"HitADCY[TrackNhits]/D");
  Tout->Branch("HitADCasym",HitADCasym,"HitADCasym[TrackNhits]/D");
  Tout->Branch("HitTmean",HitTmean,"HitTmean[TrackNhits]/D");
  Tout->Branch("HitdT",HitdT,"HitdT[TrackNhits]/D");
  Tout->Branch("HitCorrCoeff",HitCorrCoeff,"HitCorrCoeff[TrackNhits]/D");
  Tout->Branch("StripMaxCorrCoeff",StripMaxCorrCoeff,"StripMaxCorrCoeff[TrackNhits]/D");
  Tout->Branch("HitNstripX",HitNstripX,"HitNstripX[TrackNhits]/I");
  Tout->Branch("HitNstripY",HitNstripY,"HitNstripY[TrackNhits]/I");

  double xgmin_all=1e9, xgmax_all=-1e9, ygmin_all=1e9, ygmax_all=-1e9;
  map<int,double> xgmin_layer, xgmax_layer, ygmin_layer, ygmax_layer;

  //compute upper and lower limits for global X, Y coordinates within each layer for purposes of histogram definitions for output and event display:
  for( int imodule=0; imodule<nmodules; imodule++ ){
    int layer = mod_layer[imodule];
    
    double xlo_module = mod_x0[imodule] - mod_Lx[imodule]/2.0;
    double xhi_module = mod_x0[imodule] + mod_Lx[imodule]/2.0;
    double ylo_module = mod_y0[imodule] - mod_Ly[imodule]/2.0;
    double yhi_module = mod_y0[imodule] + mod_Ly[imodule]/2.0;

    if( xgmin_layer.find( layer ) == xgmin_layer.end() ){ //first time in this layer:
      xgmin_layer[layer] = xlo_module;
      xgmax_layer[layer] = xhi_module;
      ygmin_layer[layer] = ylo_module;
      ygmax_layer[layer] = yhi_module;
    } else {
      xgmin_layer[layer] = ( xlo_module < xgmin_layer[layer] ) ? xlo_module : xgmin_layer[layer];
      xgmax_layer[layer] = ( xhi_module > xgmax_layer[layer] ) ? xhi_module : xgmax_layer[layer];
      ygmin_layer[layer] = ( ylo_module < ygmin_layer[layer] ) ? ylo_module : ygmin_layer[layer];
      ygmax_layer[layer] = ( yhi_module > ygmax_layer[layer] ) ? yhi_module : ygmax_layer[layer];
    }

    xgmin_all = ( xlo_module < xgmin_all ) ? xlo_module : xgmin_all;
    xgmax_all = ( xhi_module > xgmax_all ) ? xhi_module : xgmax_all;
    
    ygmin_all = ( ylo_module < ygmin_all ) ? ylo_module : ygmin_all;
    ygmax_all = ( yhi_module > ygmax_all ) ? yhi_module : ygmax_all;
    
  }
  
  
  TCanvas *c1 = new TCanvas("c1","c1",1600,1500);
  c1->Divide( nlayers,1,.001,.001);

  TClonesArray *hframe_layers = new TClonesArray("TH2D",nlayers);
  
  
  // TH2D *hframe1 = new TH2D("hframe1","",nstripsy,-0.5,nstripsy-0.5,3*nstripsx,-0.5,3*nstripsx-0.5);
  // TH2D *hframe2 = new TH2D("hframe2","",nstripsy,-0.5,nstripsy-0.5,3*nstripsx,-0.5,3*nstripsx-0.5);
  // TH2D *hframe3 = new TH2D("hframe3","",nstripsy,-0.5,nstripsy-0.5,3*nstripsx,-0.5,3*nstripsx-0.5);
  // TH2D *hframe4 = new TH2D("hframe4","",nstripsy,-0.5,nstripsy-0.5,3*nstripsx,-0.5,3*nstripsx-0.5);

  // double strip_pitch = 0.4; //mm
  // double ymin = -220.0, ymax=220.0;
  // double xmin = -800.0, xmax=800.0;

  // int nbinsy = (ymax-ymin)/strip_pitch;
  // int nbinsx = (xmax-xmin)/strip_pitch;
  
  // TH2D *hframe1 = new TH2D("hframe1","",nbinsy,ymin,ymax,nbinsx,xmin,xmax);
  // TH2D *hframe2 = new TH2D("hframe2","",nbinsy,ymin,ymax,nbinsx,xmin,xmax);
  // TH2D *hframe3 = new TH2D("hframe3","",nbinsy,ymin,ymax,nbinsx,xmin,xmax);
  // TH2D *hframe4 = new TH2D("hframe4","",nbinsy,ymin,ymax,nbinsx,xmin,xmax);
  
  // if( eventdisplaymode != 0 ){

   
  //   //c1 = new TCanvas("c1","c1",400,1500);
  //   //c1->Divide(2,1,.001,.001);

  //   // hframe1 = new TH2D("hframe1","",nstripsy,-0.5,nstripsy-0.5,3*nstripsx,-0.5,3*nstripsx-0.5);
  //   // hframe2 = new TH2D("hframe2","",nstripsy,-0.5,nstripsy-0.5,3*nstripsx,-0.5,3*nstripsx-0.5);
  //   // hframe3 = new TH2D("hframe3","",nstripsy,-0.5,nstripsy-0.5,3*nstripsx,-0.5,3*nstripsx-0.5);
  //   // hframe4 = new TH2D("hframe4","",nstripsy,-0.5,nstripsy-0.5,3*nstripsx,-0.5,3*nstripsx-0.5);
  
  //   c1->cd(1);
  //   hframe1->Draw("colz");
  //   c1->cd(2);
  //   hframe2->Draw("colz");
  //   c1->cd(3);
  //   hframe3->Draw("colz");
  //   c1->cd(4);
  //   hframe4->Draw("colz");
  //   c1->Update();
  // } //else {
  //   //c1->Delete();
  //   //}
    
  
  TChain *C = new TChain("GEMHit");

  C->Add(filename);

  // int nch;
  // int *strip;
  // int *adc0;

  // C->SetBranchAddress("nch",&nch);
  // C->SetBranchAddress("strip",strip);
  // C->SetBranchAddress("adc0",adc0);

  GEMHit_tree_UVA_EEL *T = new GEMHit_tree_UVA_EEL(C);

  cout << "Total events = " << C->GetEntries() << endl;
  
  long nevent=0;

  fout->cd();
  
  TH1D *hNlayers_hit = new TH1D("hNlayers_hit","",nlayers+1,-0.5,nlayers+0.5);
  TH1D *hNlayers_hitX = new TH1D("hNlayers_hitX","",nlayers+1,-0.5,nlayers+0.5);
  TH1D *hNlayers_hitY = new TH1D("hNlayers_hitY","",nlayers+1,-0.5,nlayers+0.5);
  TH1D *hNlayers_hitXY = new TH1D("hNlayers_hitXY","",nlayers+1,-0.5,nlayers+0.5);
  TH2D *hNlayers_hitXvsY = new TH2D("hNlayers_hitXvsY","",nlayers+1,-0.5,nlayers+0.5,nlayers+1,-0.5,nlayers+0.5);

  TH2D *hNstripsX_layer = new TH2D("hNstripsX_layer","N strips X per layer",101,-0.5,100.5,nlayers,-0.5,nlayers-0.5);
  TH2D *hNstripsY_layer = new TH2D("hNstripsY_layer","N strips Y per layer",101,-0.5,100.5,nlayers,-0.5,nlayers-0.5);
  
  TH2D *hNclustX_layer = new TH2D("hNclustX_layer","N clusters X per layer",101,-0.5,100.5,nlayers,-0.5,nlayers-0.5);
  TH2D *hNclustY_layer = new TH2D("hNclustY_layer","N clusters Y per layer",101,-0.5,100.5,nlayers,-0.5,nlayers-0.5);

  TH1D *hNlayers_2Dclust = new TH1D("hNlayers_2Dclust","N layers with 2D cluster reconstructed",nlayers+1,-0.5,nlayers+0.5);

  TH2D *hNstripsX_module = new TH2D("hNstripsX_module","N strips X per module",101,-0.5,100.5,nmodules,-0.5,nmodules-0.5);
  TH2D *hNstripsY_module = new TH2D("hNstripsY_module","N strips Y per module",101,-0.5,100.5,nmodules,-0.5,nmodules-0.5);
  
  TH2D *hNclust_module = new TH2D("hNclust_module","N clusters 2D per module",101,-0.5,100.5,nmodules,-0.5,nmodules-0.5);
  //TH2D *hNclustY_module = new TH2D("hNclustY_module","N clusters Y per layer",101,-0.5,100.5,nmodules,-0.5,nmodules-0.5);
  
  TH1D *hclustwidth_x = new TH1D("hclustwidth_x","Width of X clusters in strips",100,0.5,100.5);
  TH1D *hclustwidth_y = new TH1D("hclustwidth_y","Width of Y clusters in strips",100,0.5,100.5);
  TH1D *hclustADCsum_x = new TH1D("hclustADCsum_x","ADC sum in X strips",200,0.0,100000.0);
  TH1D *hclustADCsum_y = new TH1D("hclustADCsum_y","ADC sum in Y strips",200,0.0,100000.0);

  TH2D *hNclustX_vs_NclustY = new TH2D("hNclustX_vs_NclustY","Nclusters X vs Y", 101,-0.5,100.5,101,-0.5,100.5);
  // TH2D *hNhitX_vs_NhitY = new TH2D("hNhitX_vs_NhitY","N hits X vs Y",101,-0.5,100.5,101,-0.5,100.5);
  
  // TH1D *hHitAxfit = new TH1D("hHitAx","Max ADC of fitted X hits",200,0.0,10000.0);
  // TH1D *hHitTX0fit = new TH1D("hHitTX0fit","T_0 of fitted X hits",200,-150.0,300.0);
  // TH1D *hHitXChi2NDF = new TH1D("hHitXChi2NDF","Chi2/NDF of fitted X hits",200,0.0,100.0);
  // TH1D *hHitXSigmaXfit = new TH1D("hHitXSigmaXfit", "sigma of Lorentzian fit to X cluster",200,0.0,5.0);
  // TH1D *hHitXTaufit = new TH1D("hHitXTaufit","Tau of fit to pulse shape in time-domain, X",200,0.0,100.0);
  

  // TH1D *hHitAyfit = new TH1D("hHitAy","Max ADC of fitted Y hits",200,0.0,10000.0);
  // TH1D *hHitTY0fit = new TH1D("hHitTY0fit","T_0 of fitted Y hits",200,-150.0,300.0);
  // TH1D *hHitYChi2NDF = new TH1D("hHitYChi2NDF","Chi2/NDF of fitted Y hits",200,0.0,100.0);
  // TH1D *hHitYSigmaYfit = new TH1D("hHitYSigmaYfit", "sigma of Lorentzian fit to Y cluster",200,0.0,5.0);
  // TH1D *hHitYTaufit = new TH1D("hHitYTaufit","Tau of fit to pulse shape in time-domain, Y",200,0.0,100.0);

  // TH2D *hHit2D_ADCXvsY = new TH2D("hHit2D_ADCXvsY","Amax X vs Y, 2D matches",250,0.0,4000.0,250,0.0,4000.0);
  // TH2D *hHit2D_T0XvsY = new TH2D("hHit2D_T0XvsY","T0 X vs Y, 2D matches",250,-150.0,300.0,250,-150.0,300.0);
  // TH1D *hHit2D_ADCdiff = new TH1D("hHit2D_ADCdiff","Amax X - Y, 2D matches",250,-1000.0,1000.0);
  // TH1D *hHit2D_T0diff = new TH1D("hHit2D_T0diff","T0 X - Y, 2D matches",250,-100.0,100.0);

  // TH1D *hHit2D_Chi2Match = new TH1D("hHit2D_Chi2Match","Chi^2 of 2D match, best 2D hit per module",250,0.0,150.0);

  // TH1D *hHit2D_ADCasym = new TH1D("hHit2D_ADCasym","(ADCX-ADCy)/(ADCX+ADCY)",250,-1.5,1.5);
  
  TH2D *hADCvsSampleAllStrips = new TH2D("hADCvsSampleAllStrips","",6,-0.5,5.5,250,-500,4500.0);
  TH2D *hADCvsSampleMaxStrip = new TH2D("hADCvsSampleMaxStrip","",6,-0.5,5.5,250,-500,4500.0);

  TH1D *hADCsumAllStrips = new TH1D("hADCsumAllStrips","",250,-2000.0,20000.0);
  TH1D *hStripTmean = new TH1D("hStripTmean","",250,0.0,150.0);
  TH1D *hStripSigmaT = new TH1D("hStripSigmaT","",250,0.0,500.0);

  TH1D *hStripTcorr = new TH1D("hStripTcorr","",250,0.0,150.0);
  
  TH1D *hN2Dmatch = new TH1D("hN2Dmatch","N 2D matches per event",101,-0.5,100.5);
  
  TH2D *hClust2D_ADCXvsY = new TH2D("hClust2D_ADCXvsY","ADC sum X vs Y, 2D cluster match",500,0.0,50000.0,500,0.0,5e4);
  TH2D *hClust2D_T0XvsY = new TH2D("hClust2D_T0XvsY", "T0 X vs T, 2D cluster match",250,0.0,150.0,250,0.,150.);

  TH1D *hClust2D_ADCdiff = new TH1D("hClust2D_ADCdiff","ADC X - Y, 2D match",500,-10000.0,10000.0);
  TH1D *hClust2D_Tdiff = new TH1D("hClust2D_Tdiff","T0 X - Y, 2D match",500,-100.,100.);
  TH1D *hClust2D_ADCasym = new TH1D("hClust2D_ADCasym","(ADC X - ADC Y)/(ADC X + ADC Y), 2D match",500,-1.5,1.5);
  TH2D *hClust2D_ADCasym_vs_ADCavg = new TH2D("hClust2D_ADCasym_vs_ADCavg","ADC X-Y asym vs X-y avg",500,0.0,5e4,500,-1.5,1.5);
  TH2D *hClust2D_ADCdiff_vs_ADCavg = new TH2D("hClust2D_ADCdiff_vs_ADCavg","ADC X-Y diff vs X-Y avg",500,0.0,5e4,500,-10000.,10000.);

  TH2D *hClust2D_NstripXvsNstripY = new TH2D("hClust2D_NstripX_vs_NstripY","Nstrips X vs Y",10,0.5,10.5,10,0.5,10.5);

  TH2D *hClust2D_Xmom_vs_NstripX = new TH2D("hClust2D_Xmom_vs_NstripX","x - xstripmax vs nstripx",maxnstripXpercluster,0.5,maxnstripXpercluster+0.5,600,-3,3);
  TH2D *hClust2D_Ymom_vs_NstripY = new TH2D("hClust2D_Ymom_vs_NstripY","y - ystripmax vs nstripy",maxnstripYpercluster,0.5,maxnstripYpercluster+0.5,600,-3,3);
  
  TH1D *hADCsumXstrip_max = new TH1D("hADCsumXstrip_max","ADC sum for max X strip, cluster on track only",1000,0.0,1.5e4);
  TH1D *hADCsumYstrip_max = new TH1D("hADCsumYstrip_max","ADC sum for max Y strip, cluster on track only",1000,0.0,1.5e4);

  TH2D *hADCsampXstrip_max = new TH2D("hADCsampXstrip_max","ADC samples for max strip, clusters on track only",6,-0.5,5.5,500,0.0,4000.0);
  TH2D *hADCsampYstrip_max = new TH2D("hADCsampYstrip_max","ADC samples for max strip, clusters on track only",6,-0.5,5.5,500,0.0,4000.0);

  TH2D *hADCsampXclust = new TH2D("hADCsampXclust","Cluster-summed X ADC samples, cluster on track",6,-0.5,5.5,500,0.0,4000.0);
  TH2D *hADCsampYclust = new TH2D("hADCsampYclust","Cluster-summed Y ADC samples, cluster on track",6,-0.5,5.5,500,0.0,4000.0);
  
  TH1D *hADCsampmax_Xstrip = new TH1D("hADCsampmax_Xstrip","max ADC sample for max strip, clusters on track",1000,0.0,3000.0);
  TH1D *hADCsampmax_Ystrip = new TH1D("hADCsampmax_Ystrip","max ADC sample for max strip, clusters on track",1000,0.0,3000.0);

  TH1D *hADCprodXYstrip_max = new TH1D("hADCprodXYstrip_max","sqrt(ADCx*ADCy) sums for max x,y strips",1000,0.0,1.5e4);
  TH1D *hADCprodXYsamp_max = new TH1D("hADCprodXYsamp_max","sqrt(ADCx*ADCy) max samples for max x,y strips",1000,0.0,3000.0);
  TH1D *hADCprodXYcluster = new TH1D("hADCprodXYcluster","sqrt(ADCX*ADCY) cluster sums",1000.0,0.0,5.0e4);
  
  TH1D *hADCsumXclust_max = new TH1D("hADCsumXclust_max","Cluster ADC sum X, clusters on tracks",1000,0.0,5e4);
  TH1D *hADCsumYclust_max = new TH1D("hADCsumYclust_max","Cluster ADC sum Y, clusters on tracks",1000,0.0,5e4);
  
  TH1D *hNtracks_found = new TH1D("hNtracks_found","Number of tracks found",11,-0.5,10.5);
  TH1D *hNhitspertrack = new TH1D("hNhitspertrack","Number of hits per track",nlayers+1,-0.5,nlayers+0.5);

  //(resid^2)/sigma^2 = chi^2: Max resid = sigma*sqrt(chi^2)

  double maxresid = 2.0*sigma_hitpos*sqrt(TrackChi2Cut);
  
  TH1D *hTrackChi2NDF = new TH1D("hTrackChi2NDF","Track chi^2/NDF",1000,0.0,TrackChi2Cut);
  TH2D *hTrackXresid_vs_layer = new TH2D("hTrackXresid_vs_layer","Track X residuals by layer",nlayers,-0.5,nlayers-0.5,200,-maxresid,maxresid); //units are mm
  TH2D *hTrackYresid_vs_layer = new TH2D("hTrackYresid_vs_layer","Track Y residuals by layer",nlayers,-0.5,nlayers-0.5,200,-maxresid,maxresid); //units are mm
  TH2D *hTrackXresid_vs_module = new TH2D("hTrackXresid_vs_module", "Track X residuals by module",nmodules,-0.5,nmodules-0.5,200,-maxresid,maxresid);
  TH2D *hTrackYresid_vs_module = new TH2D("hTrackYresid_vs_module", "Track Y residuals by module",nmodules,-0.5,nmodules-0.5,200,-maxresid,maxresid);

  TH2D *hTrackXeresid_vs_layer = new TH2D("hTrackXeresid_vs_layer","Track X eresiduals by layer",nlayers,-0.5,nlayers-0.5,200,-maxresid,maxresid); //units are mm
  TH2D *hTrackYeresid_vs_layer = new TH2D("hTrackYeresid_vs_layer","Track Y eresiduals by layer",nlayers,-0.5,nlayers-0.5,200,-maxresid,maxresid); //units are mm
  TH2D *hTrackXeresid_vs_module = new TH2D("hTrackXeresid_vs_module", "Track X eresiduals by module",nmodules,-0.5,nmodules-0.5,200,-maxresid,maxresid);
  TH2D *hTrackYeresid_vs_module = new TH2D("hTrackYeresid_vs_module", "Track Y eresiduals by module",nmodules,-0.5,nmodules-0.5,200,-maxresid,maxresid);
  
  TH2D *hTrackXY = new TH2D("hTrackXY","Track X fit vs Y fit",1000,ygmin_all-25.0, ygmax_all+25.0, 1000, xgmin_all-25.0, xgmax_all+25.0);
  TH1D *hTrackXp = new TH1D("hTrackXp","Track dx/dz fit", 1000,-1.0,1.0);
  TH1D *hTrackYp = new TH1D("hTrackYp","Track dy/dz fit", 1000,-1.0,1.0);

  TH1D *hClust_corr = new TH1D("hClust_corr","Cluster correlation coefficient, clusters on tracks",1000,-1.1,1.1);
  TH1D *hStrip_maxcor = new TH1D("hStrip_maxcor","Corr. coeff, max X and Y strips",1000,-1.1,1.1);

  TH1D *hStrip_dT = new TH1D("hStrip_dT","t_{strip}-t_{clust} (ns), strips in clusters on tracks",1000,-100.0,100.0);
  TH2D *hStrip_dT_vs_ADCsum = new TH2D("hStrip_dT_vs_ADCsum","t_{strip}-t_{clust} (ns) vs strip ADC sum",250,0.0,1.5e4,250,-100.0,100.0);

  TH1D *hStrip_dTwalkcor = new TH1D("hStrip_dTwalkcor","t_{strip}-t_{clust} (ns) w/walk correction",1000,-100,100);
  TH2D *hStrip_dTwalkcor_vs_ADC = new TH2D("hStrip_dTwalkcor_vs_ADC","t_{strip}-t_{clust} (ns) vs ADC, walk cor.",250,0.0,1.5e4,250,-100,100);
  
  TClonesArray *hdidhit_layer = new TClonesArray("TH2D",nlayers);
  TClonesArray *hshouldhit_layer = new TClonesArray("TH2D",nlayers);

  TClonesArray *hxyhit_layer = new TClonesArray("TH2D",nlayers);
  //TClonesArray *hxytrack_layer = new TClonesArray("TH2D",nlayers);
  
  long ntotal = C->GetEntries();

  //guesstimate binning for efficiency histos so that we can have ~few hundred events/bin
  double nbins_eff = double( ntotal )/400.0;

  //3.75*n^2 = ntot
  double nbinsy_eff = sqrt( nbins_eff/3.75 );
  double nbinsx_eff = 3.75*nbinsy_eff;
  
  for( int ilayer=0; ilayer<nlayers; ilayer++ ){
    TString hnametemp;
    
    new( (*hdidhit_layer)[ilayer] ) TH2D( hnametemp.Format("hdidhit_layer%d",ilayer), "Hit on track in layer", TMath::Nint(nbinsy_eff), ygmin_layer[ilayer]-10.0, ygmax_layer[ilayer]+10.0, TMath::Nint(nbinsx_eff), xgmin_layer[ilayer]-25.0,xgmax_layer[ilayer]+25.0 );
    new( (*hshouldhit_layer)[ilayer] ) TH2D( hnametemp.Format("hshouldhit_layer%d",ilayer), "Track passed through", TMath::Nint(nbinsy_eff), ygmin_layer[ilayer]-10.0, ygmax_layer[ilayer]+10.0, TMath::Nint(nbinsx_eff), xgmin_layer[ilayer]-25.0,xgmax_layer[ilayer]+25.0);

    //add a few-mm "buffer zone" at the edges of this histogram:
    new( (*hxyhit_layer)[ilayer] ) TH2D( hnametemp.Format("hxyhit_layer%d",ilayer), "X vs. Y of hits on tracks",
					 280, ygmin_layer[ilayer]-10.0, ygmax_layer[ilayer]+10.0,
					 1050, xgmin_layer[ilayer]-25.0, xgmax_layer[ilayer]+25.0 );


    new( (*hframe_layers)[ilayer] ) TH2D( hnametemp.Format("hframe_layer%d", ilayer), hnametemp.Format("Layer %d", ilayer),
					  200, ygmin_layer[ilayer]-10.0, ygmax_layer[ilayer]+10.0,
					  200, xgmin_layer[ilayer]-25.0, xgmax_layer[ilayer]+25.0 );
  }

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();


  for( int imodule=0; imodule<nmodules; imodule++ ){
    TRotation Rtemp;
    Rtemp.RotateX(mod_ax[imodule]);
    Rtemp.RotateY(mod_ay[imodule]);
    Rtemp.RotateZ(mod_az[imodule]);

    mod_Rot[imodule] = Rtemp;
    mod_Rotinv[imodule] = Rtemp.Inverse();
  }
  
  while( T->GetEntry(nevent++) && (NMAX < 0 || nevent < NMAX ) ){
    if( nevent % 1000 == 0 ) cout << nevent << endl;

    //Clustering and hit reconstruction:

    set<int> modules_hit;
    //list of unique strips fired in X and Y directions, by module:
    set<int> layers_hit;
    set<int> layers_hitX;
    set<int> layers_hitY;
    set<int> layers_hitXY;

    int NstripX_layer[nlayers];
    int NstripY_layer[nlayers];

    int Nclust2D_layer[nlayers];
    int NclustX_layer[nlayers];
    int NclustY_layer[nlayers];
    
    for( int ilay=0; ilay<nlayers; ilay++ ){
      NstripX_layer[ilay] = 0;
      NstripY_layer[ilay] = 0;
      NclustX_layer[ilay] = 0;
      NclustY_layer[ilay] = 0;
      Nclust2D_layer[ilay] = 0;
    }

    //let's consolidate all this into the "moduledata_t" data structure to avoid overhead of copying all these complicated arrays into the c struct data members:

    map<int,moduledata_t> ModData;
    
    // map<int,set<int> > mod_xstrips_hit;
    // map<int,set<int> > mod_ystrips_hit;

    // //how to store strip amplitudes: mapping of ADC sums by module, strip ID:
    // map<int,map<int,double> > ADCsum_xstrips;
    // map<int,map<int,double> > ADCsum_ystrips;

    // //how to store strip ADC samples: mapping of individual samples by module, strip ID:
    // map<int,map<int,vector<double> > > ADCsamp_xstrips;
    // map<int,map<int,vector<double> > > ADCsamp_ystrips;

    // map<int,map<int,double> > ADCmax_xstrips; //max ADC sample mapped by module, strip:
    // map<int,map<int,double> > ADCmax_ystrips; //max ADC sample mapped by module, strip:

    // map<int,map<int,int> > isampmax_xstrips;
    // map<int,map<int,int> > isampmax_ystrips;
    
    // map<int,map<int,double> > Tmean_xstrips;
    // map<int,map<int,double> > Tsigma_xstrips;
    // map<int,map<int,double> > Tmean_ystrips;
    // map<int,map<int,double> > Tsigma_ystrips;

    // map<int,map<int,double> > Tmean_xstrips_walkcor;
    // map<int,map<int,double> > Tmean_ystrips_walkcor;
    // map<int,map<int,double> > Tsigma_xstrips_walkcor;
    // map<int,map<int,double> > Tsigma_ystrips_walkcor;
    
    map<int,int> mod_maxstripX;
    map<int,double> mod_ADCmaxX;

    map<int,int> mod_maxstripY;
    map<int,double> mod_ADCmaxY;
    
    for( int ich=0; ich<T->nch; ich++ ){
      //Modified for UVA EEL GEM setup:
      int strip = (T->strip)[ich];
      int plane = (T->axis)[ich];
      //int module = (T->detID)[ich];
      int module = (T->moduleID)[ich] + 4*( (T->planeID)[ich] - 1 );
      //int layer = mod_layer[module];
      int layer = (T->planeID)[ich] - 1;
      
      modules_hit.insert( module );
      
      layers_hit.insert(layer);
      
      int ADCsamples[6];
      ADCsamples[0] = (T->adc0)[ich];
      ADCsamples[1] = (T->adc1)[ich];
      ADCsamples[2] = (T->adc2)[ich];
      ADCsamples[3] = (T->adc3)[ich];
      ADCsamples[4] = (T->adc4)[ich];
      ADCsamples[5] = (T->adc5)[ich];

      bool keepstrip=false;
      double maxsamp=0.0,sumsamp=0.0;

      int isamp_max = -1;
      
      for( int isamp=0; isamp<6; isamp++ ){
	sumsamp += ADCsamples[isamp];
	//maxsamp = (ADCsamples[isamp] > maxsamp ) ? ADCsamples[isamp] : maxsamp;
	if( isamp_max < 0 || ADCsamples[isamp] > maxsamp ){
	  maxsamp = ADCsamples[isamp];
	  isamp_max = isamp;
	}
      }

      if( maxsamp >= thresh_maxsample && sumsamp >= thresh_stripsum ) keepstrip = true;
      
      double tsum = 0.0;
      double tsum2 = 0.0;

      double tmean = 0.0;
      double tsigma = 0.0;

      double tcorr;
      
      if( plane == mod_uplaneID[module] && keepstrip ){ //x (vertical/long) axis
	//mod_xstrips_hit[module].insert( strip );

	ModData[module].xstrips.insert( strip );
	ModData[module].ADCsum_xstrips[strip] = 0.0;
	ModData[module].ADCsamp_xstrips[strip].resize(6);
	
	// ADCsum_xstrips[module][strip] = 0.0;
	// ADCsamp_xstrips[module][strip].resize(6);
	for( int isamp=0; isamp<6; isamp++ ){

	  ModData[module].ADCsamp_xstrips[strip][isamp] = ADCsamples[isamp];
	  ModData[module].ADCsum_xstrips[strip] += ADCsamples[isamp];
	  
	  // ADCsamp_xstrips[module][strip][isamp] = ADCsamples[isamp];
	  // ADCsum_xstrips[module][strip] += ADCsamples[isamp];

	  hADCvsSampleAllStrips->Fill( isamp, ADCsamples[isamp] );
	  double tsample = 12.5+25.0*isamp;

	  tsum += ADCsamples[isamp]*tsample;
	  tsum2 += ADCsamples[isamp]*pow(tsample,2);

	  if( isamp == 0 || ADCsamples[isamp] > ModData[module].ADCmax_xstrips[strip] ){
	    //    ADCmax_xstrips[module][strip] = ADCsamples[isamp];

	    ModData[module].ADCmax_xstrips[strip] = ADCsamples[isamp];
	    ModData[module].isampmax_xstrips[strip] = isamp;
	  }
	}

	tmean = tsum/ModData[module].ADCsum_xstrips[strip];
	tsigma = sqrt(fabs(tsum2/ModData[module].ADCsum_xstrips[strip]-pow(tmean,2)));
	
	tcorr = tmean - walkcor_mean_func->Eval( ModData[module].ADCsum_xstrips[strip] );
	
	ModData[module].Tmean_xstrips[strip] = tmean;
	ModData[module].Tsigma_xstrips[strip] = tsigma;

	ModData[module].Tmean_xstrips_walkcor[strip] = tcorr;
	ModData[module].Tsigma_xstrips_walkcor[strip] = walkcor_sigma_func->Eval( ModData[module].ADCsum_xstrips[strip] );
	
	hADCsumAllStrips->Fill( ModData[module].ADCsum_xstrips[strip] );

	layers_hitX.insert(layer);

	NstripX_layer[layer] += 1;

	if( mod_maxstripX.find( module ) == mod_maxstripX.end() || ModData[module].ADCsum_xstrips[strip] > mod_ADCmaxX[module] ){
	  mod_ADCmaxX[module] = ModData[module].ADCsum_xstrips[strip];
	  mod_maxstripX[module] = strip;
	}
	
      } else if( plane == mod_vplaneID[module] && keepstrip ){ //y (horizontal/short) axis
	ModData[module].ystrips.insert( strip );

	ModData[module].ADCsum_ystrips[strip] = 0.0;
	ModData[module].ADCsamp_ystrips[strip].resize(6);
	for( int isamp=0; isamp<6; isamp++ ){
	  ModData[module].ADCsamp_ystrips[strip][isamp] = ADCsamples[isamp];
	  ModData[module].ADCsum_ystrips[strip] += ADCsamples[isamp];

	  hADCvsSampleAllStrips->Fill( isamp, ADCsamples[isamp] );

	  double tsample = 12.5+25.0*isamp;
	  
	  tsum += ADCsamples[isamp]*tsample;
	  tsum2 += ADCsamples[isamp]*pow(tsample,2);

	  if( isamp == 0 || ADCsamples[isamp] > ModData[module].ADCmax_ystrips[strip] ){
	    ModData[module].ADCmax_ystrips[strip] = ADCsamples[isamp];
	    ModData[module].isampmax_ystrips[strip] = isamp;
	  }
	  
	}

	tmean = tsum/ModData[module].ADCsum_ystrips[strip];
	tsigma = sqrt(fabs(tsum2/ModData[module].ADCsum_ystrips[strip]-pow(tmean,2)));

	tcorr = tmean - walkcor_mean_func->Eval( ModData[module].ADCsum_ystrips[strip] );
	
	ModData[module].Tmean_ystrips[strip] = tmean;
	ModData[module].Tsigma_ystrips[strip] = tsigma;

	ModData[module].Tmean_ystrips_walkcor[strip] = tcorr;
	ModData[module].Tsigma_ystrips_walkcor[strip] = walkcor_sigma_func->Eval( ModData[module].ADCsum_ystrips[strip] );
	
	hADCsumAllStrips->Fill( ModData[module].ADCsum_ystrips[strip] );
	
	layers_hitY.insert(layer);

	NstripY_layer[layer] += 1;

	if( mod_maxstripY.find( module ) == mod_maxstripY.end() || ModData[module].ADCsum_ystrips[strip] > mod_ADCmaxY[module] ){
	  mod_ADCmaxY[module] = ModData[module].ADCsum_ystrips[strip];
	  mod_maxstripY[module] = strip;
	}
      }

      if( keepstrip ){
	hStripTmean->Fill(tmean);
	hStripSigmaT->Fill(tsigma);

	hStripTcorr->Fill( tcorr );
	
      }
    }


    for( set<int>::iterator ilay=layers_hitX.begin(); ilay != layers_hitX.end(); ++ilay ){
      if( layers_hitY.find( *ilay ) != layers_hitY.end() ){
	layers_hitXY.insert( *ilay );
      }
    }
    
    hNlayers_hit->Fill( layers_hit.size() );
    hNlayers_hitX->Fill( layers_hitX.size() );
    hNlayers_hitY->Fill( layers_hitY.size() );
    hNlayers_hitXY->Fill( layers_hitXY.size() );
    hNlayers_hitXvsY->Fill( layers_hitX.size(), layers_hitY.size() );

    for( int ilay=0; ilay<nlayers; ilay++ ){
      hNstripsX_layer->Fill( NstripX_layer[ilay],ilay );
      hNstripsY_layer->Fill( NstripY_layer[ilay],ilay );
    }
    
    if( layers_hitXY.size() >= 3 ){ //enough layers hit to (possibly) form a track: only bother with clustering and attempted track finding if this is the case: 

      if( eventdisplaymode != 0 ){
	for( int ilayer=0; ilayer<nlayers; ilayer++ ){
	  ( (TH2D*) (*hframe_layers)[ilayer] )->Reset();
	}
      }

      //we should probably restructure the earlier part of the code to eliminate the overhead of copying all this info:
      
      map<int,clusterdata_t> mod_clusters; //mapping of found clusters by module;

      //map<int,int> N2D_match_layer;
      
      for(set<int>::iterator imod=modules_hit.begin(); imod != modules_hit.end(); ++imod ){
	//moduledata_t datatemp;
	clusterdata_t clusttemp;

	// ModData[*imod].modindex = *imod;
	// datatemp.modindex = *imod;
	// datatemp.layerindex = mod_layer[*imod];
	
	// datatemp.xstrips = mod_xstrips_hit[*imod];
	// datatemp.ystrips = mod_ystrips_hit[*imod];
	// datatemp.ADCsum_xstrips = ADCsum_xstrips[*imod];
	// datatemp.ADCsum_ystrips = ADCsum_ystrips[*imod];
	// datatemp.ADCsamp_xstrips = ADCsamp_xstrips[*imod];
	// datatemp.ADCsamp_ystrips = ADCsamp_ystrips[*imod];
	// datatemp.ADCmax_xstrips = ADCmax_xstrips[*imod];
	// datatemp.ADCmax_ystrips = ADCmax_ystrips[*imod];
	// datatemp.isampmax_xstrips = isampmax_xstrips[*imod];
	// datatemp.isampmax_ystrips = isampmax_ystrips[*imod];
	// datatemp.Tmean_xstrips = Tmean_xstrips[*imod];
	// datatemp.Tsigma_xstrips = Tsigma_xstrips[*imod];
	// datatemp.Tmean_ystrips = Tmean_ystrips[*imod];
	// datatemp.Tsigma_ystrips = Tsigma_ystrips[*imod];
	// datatemp.Tmean_xstrips_walkcor = Tmean_xstrips_walkcor[*imod];
	// datatemp.Tsigma_xstrips_walkcor = Tsigma_xstrips_walkcor[*imod];
	// datatemp.Tmean_ystrips_walkcor = Tmean_ystrips_walkcor[*imod];
	// datatemp.Tsigma_ystrips_walkcor = Tsigma_ystrips_walkcor[*imod];
	
	//for( map<int,int>::iterator jmod=mod_maxstripX.begin(); jmod!=mod_maxstripX.end(); ++jmod ){
	int module = *imod;
	int layer = mod_layer[module];

	hNstripsX_module->Fill( ModData[module].xstrips.size(), module );
	hNstripsY_module->Fill( ModData[module].ystrips.size(), module );
	
	ModData[module].modindex = module;
	ModData[module].layerindex = layer;
	int strip = mod_maxstripX[module];
	if( fabs( ModData[module].Tmean_xstrips[strip] - 64.3 ) <= 2.5*16.3 &&
	    fabs( ModData[module].Tsigma_xstrips[strip] -36.25) <= 2.5*5.7 ){
	  for( int isamp=0; isamp<6; isamp++ ){
	    
	    hADCvsSampleMaxStrip->Fill( isamp, ModData[module].ADCsamp_xstrips[strip][isamp] );
	  }
	}
	
	strip = mod_maxstripY[module];
	//for( map<int,int>::iterator imod=mod_maxstripY.begin(); imod!=mod_maxstripY.end(); ++imod ){
	//int module = imod->first;
	//int strip = imod->second;
	if( fabs( ModData[module].Tmean_ystrips[strip] - 64.3 ) <= 2.5*16.3 &&
	    fabs( ModData[module].Tsigma_ystrips[strip] -36.25) <= 2.5*5.7 ){
	  for( int isamp=0; isamp<6; isamp++ ){
	    hADCvsSampleMaxStrip->Fill( isamp, ModData[module].ADCsamp_ystrips[strip][isamp] );
	  }
	}
	

	//	cout << "starting cluster finding, event " << T->evtID << ", module " << module << ". (nstrip X, nstrip Y) = (" << datatemp.xstrips.size()
	//    << ", " << datatemp.ystrips.size() << ")" << endl;
	int clusterflag = find_clusters_by_module( ModData[module], clusttemp );

	if( clusterflag != 0 ){
	  cout << "event " << T->evtID << ", module " << module << " too noisy, gave up" << endl;
	}
	
	//	cout << "ending cluster finding " << endl;

	hNclust_module->Fill(clusttemp.nclust2D, module );
	
	mod_clusters[module] = clusttemp;

	NclustX_layer[mod_layer[module]]+=clusttemp.nclustx;
	NclustY_layer[mod_layer[module]]+=clusttemp.nclusty;

	Nclust2D_layer[mod_layer[module]]+=clusttemp.nclust2D;
	
	for( int iclustx=0; iclustx<clusttemp.nclustx; iclustx++ ){
	  hclustwidth_x->Fill( clusttemp.nstripx[iclustx] );
	  hclustADCsum_x->Fill( clusttemp.totalchargex[iclustx] );
	}

	for( int iclusty=0; iclusty<clusttemp.nclusty; iclusty++ ){
	  hclustwidth_y->Fill( clusttemp.nstripy[iclusty] );
	  hclustADCsum_y->Fill( clusttemp.totalchargey[iclusty] );
	}

	hNclustX_vs_NclustY->Fill( clusttemp.nclusty, clusttemp.nclustx );

	//	hNhitX_vs_NhitY->Fill( clusttemp.nhity1D, clusttemp.nhitx1D );
	
	// for( int ihitx=0; ihitx<clusttemp.nhitx1D; ihitx++ ){
	//   if( clusttemp.nstripx[clusttemp.clustidx_xhit[ihitx]]>2 ){
	//     hHitAxfit->Fill( clusttemp.Axhit1D[ihitx] );
	//     hHitXChi2NDF->Fill( clusttemp.chi2ndf_hitx1D[ihitx] );

	//     if( clusttemp.Axhit1D[ihitx] < 2200.0 ){
	//       hHitTX0fit->Fill( clusttemp.tx0hit1D[ihitx] );
	   
	//       //if( clusttemp.nstripx[clusttemp.clustidx_xhit[ihitx]] > 1 )
	//       hHitXTaufit->Fill( clusttemp.taux0hit1D[ihitx] );
	//       hHitXSigmaXfit->Fill( clusttemp.sigx0hit1D[ihitx] );
	//     }
	//   }
	// }
	
	// for( int ihity=0; ihity<clusttemp.nhity1D; ihity++ ){
	//   if( clusttemp.nstripy[clusttemp.clustidx_yhit[ihity]] > 2 ){
	//     hHitYChi2NDF->Fill( clusttemp.chi2ndf_hity1D[ihity] );
	//     hHitAyfit->Fill( clusttemp.Ayhit1D[ihity] );
	    
	//     if( clusttemp.Ayhit1D[ihity] < 2200.0 ){ //cut out "saturation" events from these hits:
	//       hHitTY0fit->Fill( clusttemp.ty0hit1D[ihity] );	      
	//       hHitYSigmaYfit->Fill( clusttemp.sigy0hit1D[ihity] );
	//       hHitYTaufit->Fill( clusttemp.tauy0hit1D[ihity] );
	//     }
	//   }
	// }

	//	for( int ihit2D=0; ihit2D<clusttemp.nhit2D; ihit2D++ ){
	// for( int ihit2D=0; ihit2D<TMath::Min(1,clusttemp.nhit2D); ihit2D++ ){

	//   if( clusttemp.nstripx[clusttemp.clustidx_xhit[clusttemp.ixhit2D[ihit2D]]] > 2 &&
	//       clusttemp.nstripy[clusttemp.clustidx_yhit[clusttemp.iyhit2D[ihit2D]]] > 2 ){
	  
	//     hHit2D_ADCXvsY->Fill( clusttemp.Axhit1D[clusttemp.ixhit2D[ihit2D]],
	// 			  clusttemp.Ayhit1D[clusttemp.iyhit2D[ihit2D]] );
	//     hHit2D_T0XvsY->Fill( clusttemp.tx0hit1D[clusttemp.ixhit2D[ihit2D]],
	// 			 clusttemp.ty0hit1D[clusttemp.iyhit2D[ihit2D]] );
	    
	//     hHit2D_ADCdiff->Fill( clusttemp.dAhit2D[ihit2D] );
	//     hHit2D_T0diff->Fill( clusttemp.dthit2D[ihit2D] );
	    
	//     hHit2D_Chi2Match->Fill( clusttemp.chi2match2D[ihit2D] );
	    
	//     hHit2D_ADCasym->Fill( clusttemp.dAhit2D[ihit2D]/(2.*clusttemp.Ahit2D[ihit2D]) );
	//   }
	// }

	//cout << "nclust2D = " << clusttemp.nclust2D << endl;
	
	// for( int iclust2D=0; iclust2D<clusttemp.nclust2D; iclust2D++ ){
	//   //if( clusttemp.nstripx2D[iclust2D]>1&&clusttemp.nstripy2D[iclust2D]>1 ){
	//   if( clusttemp.itrack_clust2D[iclust2D] >= 0 ){ //only fill this information if the cluster is on a found track:
	//     hClust2D_ADCXvsY->Fill( clusttemp.totalchargex[clusttemp.ixclust2D[iclust2D]],
	// 			    clusttemp.totalchargey[clusttemp.iyclust2D[iclust2D]] );
	//     hClust2D_T0XvsY->Fill( clusttemp.txmean[clusttemp.ixclust2D[iclust2D]],
	// 			   clusttemp.tymean[clusttemp.iyclust2D[iclust2D]] );
	//     hClust2D_ADCdiff->Fill( clusttemp.dEclust2D[iclust2D] );
	//     hClust2D_Tdiff->Fill( clusttemp.dtclust2D[iclust2D] );
	//     hClust2D_ADCasym->Fill( clusttemp.dEclust2D[iclust2D]/(2.0*clusttemp.Eclust2D[iclust2D]) );
	//     hClust2D_ADCasym_vs_ADCavg->Fill( clusttemp.Eclust2D[iclust2D], clusttemp.dEclust2D[iclust2D]/(2.0*clusttemp.Eclust2D[iclust2D]) );
	//     hClust2D_ADCdiff_vs_ADCavg->Fill( clusttemp.Eclust2D[iclust2D], clusttemp.dEclust2D[iclust2D] );
	  
	//     hStrip_maxcor->Fill( clusttemp.CorrCoeff2D[iclust2D] );

	//     hClust2D_NstripXvsNstripY->Fill( clusttemp.nstripy2D, clusttemp.nstripx2D );
	//   }
	//   //}
	// }
      }

      int nlayers_with_2Dclust = 0;
      for( int ilay=0; ilay<nlayers; ilay++ ){
	if( Nclust2D_layer[ilay] > 0 ) nlayers_with_2Dclust++;
      }

      hNlayers_2Dclust->Fill( nlayers_with_2Dclust );
      
      trackdata_t tracktemp;

      tracktemp.ntracks = 0;

      //      cout << "Finding tracks, event..." << T->evtID << endl;

      if( nlayers_with_2Dclust >= 3 ){
      
	find_tracks( mod_clusters, tracktemp );
	
	//	cout << "track finding successful..." << endl;
	
	hNtracks_found->Fill( tracktemp.ntracks );
	
	for( int ilay=0; ilay<nlayers; ilay++ ){
	  Nclustperlayer[ilay] = Nclust2D_layer[ilay];
	}

	Ntracks = tracktemp.ntracks;
      
	//for( int itrack=0; itrack<tracktemp.ntracks; itrack++ ){
	for( int itrack=0; itrack<TMath::Min(tracktemp.ntracks,1); itrack++ ){
	  hNhitspertrack->Fill( tracktemp.nhitsontrack[itrack] );
	  hTrackChi2NDF->Fill( tracktemp.Chi2NDFtrack[itrack] );

	  if( itrack == 0 ){ 
	    TrackXp = tracktemp.Xptrack[itrack];
	    TrackYp = tracktemp.Yptrack[itrack];
	    TrackX = tracktemp.Xtrack[itrack];
	    TrackY = tracktemp.Ytrack[itrack];

	    TrackNhits = tracktemp.nhitsontrack[itrack];
	    TrackChi2NDF = tracktemp.Chi2NDFtrack[itrack];

	    for( int ilay=0; ilay<nlayers; ilay++ ){
	      double xtracktemp = TrackX + zavg_layer[ilay]*TrackXp;
	      double ytracktemp = TrackY + zavg_layer[ilay]*TrackYp;

	      ( (TH2D*) (*hshouldhit_layer)[ilay] )->Fill( ytracktemp, xtracktemp );
	    }
	  }

	  
	
	  //Only fill residuals if we have all four layers firing (may help clarify the situation):
	  //if( tracktemp.nhitsontrack[itrack] == 4 ){
	  //&& tracktemp.Chi2NDFtrack[itrack] < 1000.0 ){
	  for( int ihit=0; ihit<tracktemp.nhitsontrack[itrack]; ihit++ ){
	    int module= tracktemp.modlist_track[itrack][ihit];
	    int layer = mod_layer[module];

	    int iclust2D = tracktemp.hitlist_track[itrack][ihit];
	    
	    //	  if( tracktemp.nhitsontrack[itrack] == 4 && Nclustperlayer[layer] == 1 ){
	    hTrackXresid_vs_layer->Fill( layer, -tracktemp.residx_hits[itrack][ihit] );
	    hTrackYresid_vs_layer->Fill( layer, -tracktemp.residy_hits[itrack][ihit] );
	    
	    hTrackXresid_vs_module->Fill( module, -tracktemp.residx_hits[itrack][ihit] );
	    hTrackYresid_vs_module->Fill( module, -tracktemp.residy_hits[itrack][ihit] );
	    //}
	    clusterdata_t clusttemp = mod_clusters[module];
	    //Fill 2D cluster properties histograms IFF they are on tracks:

	    //if we have at least 4 hits on the track, then we can compute "exclusive residuals":
	    
	    if( tracktemp.nhitsontrack[itrack] > 3 ){
	      
	      double sumxhits=0.0,sumyhits=0.0,sumzhits=0.0,sumxzhits=0.0,sumyzhits=0.0,sumz2hits=0.0;
	      
	      double xtrtemp,ytrtemp,ztrtemp,xptrtemp,yptrtemp; //temporary storage for track parameters: 

	      double varx,vary,varxp,varyp,covxxp,covyyp;

	      //loop over all hits OTHER than ihit, compute sums needed for track fit:
	      for( int jhit=0; jhit<tracktemp.nhitsontrack[itrack]; jhit++ ){
		if( jhit != ihit ){

		  int modj = tracktemp.modlist_track[itrack][jhit];
		  int layj = mod_layer[modj];
		  
		  double xhittemp=mod_clusters[modj].xglobal2D[tracktemp.hitlist_track[itrack][jhit]];
		  double yhittemp=mod_clusters[modj].yglobal2D[tracktemp.hitlist_track[itrack][jhit]];
		  double zhittemp=mod_clusters[modj].zglobal2D[tracktemp.hitlist_track[itrack][jhit]];
		  
		  sumxhits += xhittemp;
		  sumyhits += yhittemp;
		  sumzhits += zhittemp;
		  sumxzhits += xhittemp*zhittemp;
		  sumyzhits += yhittemp*zhittemp;
		  sumz2hits += pow(zhittemp,2);
		}
	      }

	      int nhittemp = tracktemp.nhitsontrack[itrack]-1;

	      double denom = (sumz2hits*nhittemp - pow(sumzhits,2));
	      xptrtemp = (nhittemp*sumxzhits-sumxhits*sumzhits)/denom;
	      yptrtemp = (nhittemp*sumyzhits-sumyhits*sumzhits)/denom;
	      xtrtemp = (sumz2hits*sumxhits-sumzhits*sumxzhits)/denom;
	      ytrtemp = (sumz2hits*sumyhits-sumzhits*sumyzhits)/denom;

	      //now compute track position at layer in local coordinates:
	      double uhittemp = clusttemp.xclust2Dcorr[iclust2D];
	      double vhittemp = clusttemp.yclust2Dcorr[iclust2D];
	      double zhittemp = clusttemp.zglobal2D[iclust2D];
	      
	      TVector3 trackpos_global( xtrtemp+xptrtemp*zhittemp, ytrtemp+yptrtemp*zhittemp, zhittemp );
	      TVector3 modcenter_global( mod_x0[module], mod_y0[module], mod_z0[module] );
		   
		    
	      TVector3 trackpos_local = mod_Rotinv[module]*(trackpos_global - modcenter_global);

	      //Compute the local "u" and "v" coordinates of the track within the module, to allow for the
	      //possibility of different strip orientations than X/Y:
	      double utracktemp = trackpos_local.X()*mod_Pxu[module] + trackpos_local.Y()*mod_Pyu[module];
	      double vtracktemp = trackpos_local.X()*mod_Pxv[module] + trackpos_local.Y()*mod_Pyv[module];

	      tracktemp.eresidx_hits[itrack][ihit] = uhittemp - utracktemp;
	      tracktemp.eresidy_hits[itrack][ihit] = vhittemp - vtracktemp;

	      hTrackXeresid_vs_layer->Fill( layer, utracktemp - uhittemp );
	      hTrackYeresid_vs_layer->Fill( layer, vtracktemp - vhittemp );

	      hTrackXeresid_vs_module->Fill( module, utracktemp - uhittemp );
	      hTrackYeresid_vs_module->Fill( module, vtracktemp - vhittemp );
	      
	    }
	    
	  
	    hClust2D_ADCXvsY->Fill( clusttemp.totalchargex[clusttemp.ixclust2D[iclust2D]],
				    clusttemp.totalchargey[clusttemp.iyclust2D[iclust2D]] );
	    hClust2D_T0XvsY->Fill( clusttemp.txmean[clusttemp.ixclust2D[iclust2D]],
				   clusttemp.tymean[clusttemp.iyclust2D[iclust2D]] );
	    hClust2D_ADCdiff->Fill( clusttemp.dEclust2D[iclust2D] );
	    hClust2D_Tdiff->Fill( clusttemp.dtclust2D[iclust2D] );
	    hClust2D_ADCasym->Fill( clusttemp.dEclust2D[iclust2D]/(2.0*clusttemp.Eclust2D[iclust2D]) );
	    hClust2D_ADCasym_vs_ADCavg->Fill( sqrt( clusttemp.totalchargex[clusttemp.ixclust2D[iclust2D]]*
						    clusttemp.totalchargey[clusttemp.iyclust2D[iclust2D]] ),
					      clusttemp.dEclust2D[iclust2D]/(2.0*clusttemp.Eclust2D[iclust2D]) );
	    hClust2D_ADCdiff_vs_ADCavg->Fill( clusttemp.Eclust2D[iclust2D], clusttemp.dEclust2D[iclust2D] );
	    
	    hClust_corr->Fill( clusttemp.CorrCoeff2D[iclust2D] );
	    hStrip_maxcor->Fill( clusttemp.CorrCoeffMaxStrips[iclust2D] );
	    
	    hClust2D_NstripXvsNstripY->Fill( clusttemp.nstripy2D[iclust2D], clusttemp.nstripx2D[iclust2D] );

	    hADCsumXclust_max->Fill( clusttemp.totalchargex[clusttemp.ixclust2D[iclust2D]] );
	    hADCsumYclust_max->Fill( clusttemp.totalchargey[clusttemp.iyclust2D[iclust2D]] );

	    hADCsumXstrip_max->Fill( ModData[module].ADCsum_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]] );
	    hADCsumYstrip_max->Fill( ModData[module].ADCsum_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]] );
	    
	    for( int isamp=0; isamp<6; isamp++ ){
	      hADCsampXstrip_max->Fill( isamp, ModData[module].ADCsamp_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]][isamp] );
	      hADCsampYstrip_max->Fill( isamp, ModData[module].ADCsamp_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]][isamp] );

	      hADCsampXclust->Fill( isamp, clusttemp.ADCsamp_xclust[iclust2D][isamp] );
	      hADCsampYclust->Fill( isamp, clusttemp.ADCsamp_yclust[iclust2D][isamp] );
	    }

	    hADCsampmax_Xstrip->Fill( ModData[module].ADCsamp_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]][ModData[module].isampmax_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]]] );
	    hADCsampmax_Ystrip->Fill( ModData[module].ADCsamp_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]][ModData[module].isampmax_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]]] );

	    hADCprodXYstrip_max->Fill( sqrt( ModData[module].ADCsum_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]] * 
					     ModData[module].ADCsum_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]] ) );
	    hADCprodXYsamp_max->Fill( sqrt( ModData[module].ADCsamp_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]][ModData[module].isampmax_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]]] *
					    ModData[module].ADCsamp_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]][ModData[module].isampmax_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]]] ) );

	    hADCprodXYcluster->Fill( sqrt( clusttemp.totalchargex[clusttemp.ixclust2D[iclust2D]]*
					   clusttemp.totalchargey[clusttemp.iyclust2D[iclust2D]] ) );

	    
	    
	    //loop over strips in cluster and fill strip deltaT histograms:
	    for( int ixstrip=clusttemp.ixstriplo[clusttemp.ixclust2D[iclust2D]]; ixstrip<=clusttemp.ixstriphi[clusttemp.ixclust2D[iclust2D]]; ixstrip++ ){
	      if( ModData[module].Tmean_xstrips.find( ixstrip ) != ModData[module].Tmean_xstrips.end() ){
		hStrip_dT->Fill( ModData[module].Tmean_xstrips[ixstrip] - clusttemp.tclust2D[iclust2D] );
		hStrip_dT_vs_ADCsum->Fill( ModData[module].ADCsum_xstrips[ixstrip], ModData[module].Tmean_xstrips[ixstrip] - clusttemp.tclust2D[iclust2D] );

		hStrip_dTwalkcor->Fill( ModData[module].Tmean_xstrips_walkcor[ixstrip] - clusttemp.tclust2Dwalkcorr[iclust2D] );
		hStrip_dTwalkcor_vs_ADC->Fill( ModData[module].ADCsum_xstrips[ixstrip], ModData[module].Tmean_xstrips_walkcor[ixstrip] - clusttemp.tclust2Dwalkcorr[iclust2D] );
		
	      }
	    }

	    //loop over strips in cluster and fill strip deltaT histograms:
	    for( int iystrip=clusttemp.iystriplo[clusttemp.iyclust2D[iclust2D]]; iystrip<=clusttemp.iystriphi[clusttemp.iyclust2D[iclust2D]]; iystrip++ ){
	      if( ModData[module].Tmean_ystrips.find( iystrip ) != ModData[module].Tmean_ystrips.end() ){
		hStrip_dT->Fill( ModData[module].Tmean_ystrips[iystrip] - clusttemp.tclust2D[iclust2D] );
		hStrip_dT_vs_ADCsum->Fill( ModData[module].ADCsum_ystrips[iystrip], ModData[module].Tmean_ystrips[iystrip] - clusttemp.tclust2D[iclust2D] );

		hStrip_dTwalkcor->Fill( ModData[module].Tmean_ystrips_walkcor[iystrip] - clusttemp.tclust2Dwalkcorr[iclust2D] );
		hStrip_dTwalkcor_vs_ADC->Fill( ModData[module].ADCsum_ystrips[iystrip], ModData[module].Tmean_ystrips_walkcor[iystrip] - clusttemp.tclust2Dwalkcorr[iclust2D] );
	      }
	    }
	    
	    if( itrack == 0 ){ //set tree output variables:

	      HitModule[ihit] = module;
	      HitLayer[ihit] = layer;
	    
	      HitXlocal[ihit] = clusttemp.xclust2Dcorr[iclust2D];
	      HitYlocal[ihit] = clusttemp.yclust2Dcorr[iclust2D];
	      HitXresid[ihit] = tracktemp.residx_hits[itrack][ihit];
	      HitYresid[ihit] = tracktemp.residy_hits[itrack][ihit];
	      HitXresidE[ihit] = tracktemp.eresidx_hits[itrack][ihit];
	      HitYresidE[ihit] = tracktemp.eresidy_hits[itrack][ihit];

	      // //re-compute global coordinates:
	      // TVector3 hitpos_local( HitXlocal[ihit], HitYlocal[ihit], 0.0 );

	      // TRotation R;
	      // //Order of the rotations affects results somewhat because rotations don't commute. But
	      // // if they are sufficiently small, they are almost independent:
	      // R.RotateX(mod_ax[module]);
	      // R.RotateY(mod_ay[module]);
	      // R.RotateZ(mod_az[module]);

	      // TVector3 modulecenter_global( mod_x0[module], mod_y0[module], mod_z0[module] );
	      // TVector3 hitpos_global = R * hitpos_local + modulecenter_global;

	      TVector3 hitpos_global( clusttemp.xglobal2D[iclust2D],
				      clusttemp.yglobal2D[iclust2D],
				      clusttemp.zglobal2D[iclust2D] );
	      
	      HitXglobal[ihit] = hitpos_global.X();
	      HitYglobal[ihit] = hitpos_global.Y();
	      HitZglobal[ihit] = hitpos_global.Z();

	      HitXmom[ihit] = clusttemp.xmom2D[iclust2D];
	      HitYmom[ihit] = clusttemp.ymom2D[iclust2D];
	      
	      //if there IS a hit: fill "did hit" histos with TRACK coordinates, not hit coordinates, to avoid bin migration effects:
	      
	      //( (TH2D*) (*hdidhit_layer)[layer] )->Fill( HitYglobal[ihit], HitXglobal[ihit] );

	      //Need to fill this with identical coordinates to how "should hit" is done to avoid bin migration effects:
	      ( (TH2D*) (*hdidhit_layer)[layer] )->Fill( TrackY + zavg_layer[layer]*TrackYp,
							TrackX + zavg_layer[layer]*TrackXp );

	      ( (TH2D*) (*hxyhit_layer)[layer] )->Fill( HitYglobal[ihit], HitXglobal[ihit] );
	      
	      HitSigX[ihit] = clusttemp.xsigma[clusttemp.ixclust2D[iclust2D]];
	      HitSigY[ihit] = clusttemp.ysigma[clusttemp.iyclust2D[iclust2D]];
	      HitADCX[ihit] = clusttemp.totalchargex[clusttemp.ixclust2D[iclust2D]];
	      HitADCY[ihit] = clusttemp.totalchargey[clusttemp.iyclust2D[iclust2D]];
	      HitADCasym[ihit] = clusttemp.dEclust2D[iclust2D]/(2.0*clusttemp.Eclust2D[iclust2D]);
	      HitTmean[ihit] = clusttemp.tclust2D[iclust2D];
	      HitdT[ihit] = clusttemp.dtclust2D[iclust2D];
	      HitCorrCoeff[ihit] = clusttemp.CorrCoeff2D[iclust2D];
	      StripMaxCorrCoeff[ihit] = clusttemp.CorrCoeffMaxStrips[iclust2D];
	      HitNstripX[ihit] = clusttemp.nstripx2D[iclust2D];
	      HitNstripY[ihit] = clusttemp.nstripy2D[iclust2D];

	      hClust2D_Xmom_vs_NstripX->Fill( HitNstripX[ihit], HitXmom[ihit] );
	      hClust2D_Ymom_vs_NstripY->Fill( HitNstripY[ihit], HitYmom[ihit] );
					    
	    }
	  }
      

	  hTrackXY->Fill( tracktemp.Ytrack[itrack],tracktemp.Xtrack[itrack] );
	  hTrackXp->Fill( tracktemp.Xptrack[itrack] );
	  hTrackYp->Fill( tracktemp.Yptrack[itrack] );
	}
      }

      if( tracktemp.ntracks > 0 ){
	Tout->Fill();
      }
      
      if( eventdisplaymode != 0 ){

	TLine Ltemp;

	Ltemp.SetLineWidth(1);
	
	double stripADCmax=1.2e4;
	int ncolors = gStyle->GetNumberOfColors();

	for( int ilayer=0; ilayer<nlayers; ilayer++ ){
	  c1->cd( ilayer+1 );

	  ( (TH2D*) (*hframe_layers)[ilayer] )->Draw();
	}
	
	for( set<int>::iterator imod=modules_hit.begin(); imod != modules_hit.end(); ++imod ){
	  
	  int layer = mod_layer[*imod];
	  int module = *imod;

	  set<int> xstrips = ModData[module].xstrips;
	  set<int> ystrips = ModData[module].ystrips;
	  
	  c1->cd(layer+1);
	  
	  //clusterdata_t clusttemp = mod_clusters[module];
	  
	  for( set<int>::iterator istrip=xstrips.begin(); istrip != xstrips.end(); ++istrip ){
	    int strip = *istrip;
	    //int binx = strip + 1280*( module % 3 );

	    //double xlocal = (strip+0.5 - 0.5*nstripsx)*0.4;
	    //Local U coordinate of the center of the strip according to the convention:
	    double ulocal = (strip + 0.5 - 0.5*mod_nstripsu[module])*mod_ustrip_pitch[module];

	    // The X strips are along the line of constant "U" =
	    TVector3 uhat( mod_Pxu[module], mod_Pyu[module], 0.0);
	    TVector3 zaxis(0,0,1);

	    //Unit vector perp. to U; i.e., ALONG strip:
	    TVector3 uperphat = zaxis.Cross(uhat).Unit();
	    
	    TVector3 strip_center_pos = ulocal*uhat;
	    
	    if( fabs( strip_center_pos.X() ) <= mod_Lx[module]/2.0 &&
		fabs( strip_center_pos.Y() ) <= mod_Ly[module]/2.0 ){
	      //Then we can draw this strip:

	      if( uhat.X() != 0 ){ //strip is not along Y, compute upper and lower limits in X:
		double xmin_strip = strip_center_pos.X() + uperphat.X()/uperphat.Y() * (-mod_Ly[module]/2.0-strip_center_pos.Y());
		double xmax_strip = strip_center_pos.X() + uperphat.X()/uperphat.Y() * (mod_Ly[module]/2.-strip_center_pos.Y());

		double ymin_strip = -mod_Ly[module]/2.0;
		double ymax_strip = mod_Ly[module]/2.0;
		
		xmin_strip = (xmin_strip < -mod_Lx[module]/2.0) ? -mod_Lx[module]/2.0 : xmin_strip;
		xmax_strip = (xmax_strip > mod_Lx[module]/2.0) ? mod_Lx[module]/2.0 : xmax_strip;

		if( xmin_strip == -mod_Lx[module]/2.0 ) ymin_strip = strip_center_pos.Y() + uperphat.Y()/uperphat.X() * ( xmin_strip - strip_center_pos.X() );
		if( xmax_strip == mod_Lx[module]/2.0 ) ymax_strip = strip_center_pos.Y() + uperphat.Y()/uperphat.X() * ( xmax_strip - strip_center_pos.X() );

		//These are local coordinates: need to convert to global coordinates; it should be sufficient to convert the two points:
		TVector3 point1local(xmin_strip, ymin_strip, 0.0);
		TVector3 point2local(xmax_strip,ymax_strip, 0.0 );
		TVector3 modcenter(mod_x0[module],mod_y0[module],mod_z0[module]);
		TVector3 point1global = mod_Rot[module] * point1local + modcenter;
		TVector3 point2global = mod_Rot[module] * point2local + modcenter;

		//How to implement a logarithmic color scale: ADCstrip over ADCmax varies from threshold/max up to 1:
		double logADCmin = log(thresh_stripsum/stripADCmax);
		double logADCmax = 0.0;

		int logADCbin = int( (log(ModData[module].ADCsum_xstrips[strip]/stripADCmax)-logADCmin)/(logADCmax-logADCmin)*double(ncolors) );

		//cout << "logADCmin, ADC, logADCbin = " << logADCmin << ", " << ADCsum_xstrips[module][strip] << ", " << logADCbin << ", color = " << gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,logADCbin))) << endl;
		int ADCbin = int( (ModData[module].ADCsum_xstrips[strip]-thresh_stripsum)/(stripADCmax-thresh_stripsum)*double(ncolors) );
		
		Ltemp.SetLineColor( gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,ADCbin))));

		Ltemp.DrawLine( point1global.Y(), point1global.X(), point2global.Y(), point2global.X() );
		
	      } else { //strip IS along Y:
		double xmax_strip = strip_center_pos.X();
		double xmin_strip = strip_center_pos.X();
		double ymax_strip = mod_Ly[module]/2.0;
		double ymin_strip = -mod_Ly[module]/2.0;

		//These are local coordinates: need to convert to global coordinates; it should be sufficient to convert the two points:
		TVector3 point1local(xmin_strip, ymin_strip, 0.0);
		TVector3 point2local(xmax_strip,ymax_strip, 0.0 );
		TVector3 modcenter(mod_x0[module],mod_y0[module],mod_z0[module]);
		TVector3 point1global = mod_Rot[module] * point1local + modcenter;
		TVector3 point2global = mod_Rot[module] * point2local + modcenter;

		//How to implement a logarithmic color scale: ADCstrip over ADCmax varies from threshold/max up to 1:
		double logADCmin = log(thresh_stripsum/stripADCmax);
		double logADCmax = 0.0;
		
		int logADCbin = int( (log(ModData[module].ADCsum_xstrips[strip]/stripADCmax)-logADCmin)/(logADCmax-logADCmin)*double(ncolors) );

		//cout << "logADCmin, ADC, logADCbin = " << logADCmin << ", " << ADCsum_xstrips[module][strip] << ", " << logADCbin << ", color = " << gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,logADCbin))) << endl;
		int ADCbin = int( (ModData[module].ADCsum_xstrips[strip]-thresh_stripsum)/(stripADCmax-thresh_stripsum)*double(ncolors) );
		
		Ltemp.SetLineColor( gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,ADCbin))));
		
	    

		Ltemp.DrawLine( point1global.Y(), point1global.X(), point2global.Y(), point2global.X() );
	      }
	      
	    }
	    
	  }
	    
	  for( set<int>::iterator istrip=ystrips.begin(); istrip != ystrips.end(); ++istrip ){
	    int strip = *istrip;
	    //int binx = strip + 1280*( module % 3 );

	    //double xlocal = (strip+0.5 - 0.5*nstripsx)*0.4;
	    //Local U coordinate of the center of the strip according to the convention:
	    double vlocal = (strip + 0.5 - 0.5*mod_nstripsv[module])*mod_vstrip_pitch[module];

	    // The X strips are along the line of constant "U" =
	    TVector3 vhat( mod_Pxv[module], mod_Pyv[module], 0.0);
	    TVector3 zaxis(0,0,1);

	    //Unit vector perp. to U; i.e., ALONG strip:
	    TVector3 vperphat = zaxis.Cross(vhat).Unit();
	    
	    TVector3 strip_center_pos = vlocal*vhat;
	    
	    if( fabs( strip_center_pos.X() ) <= mod_Lx[module]/2.0 &&
		fabs( strip_center_pos.Y() ) <= mod_Ly[module]/2.0 ){
	      //Then we can draw this strip:

	      //I think this covers all possible scenarios, SHOULD avoid divide-by-zero errors:
	      if( vhat.X() != 0 ){ //strip is not along Y, compute upper and lower limits in X:
		double xmin_strip = strip_center_pos.X() + vperphat.X()/vperphat.Y() * (-mod_Ly[module]/2.0-strip_center_pos.Y());
		double xmax_strip = strip_center_pos.X() + vperphat.X()/vperphat.Y() * (mod_Ly[module]/2.-strip_center_pos.Y());

		double ymin_strip = -mod_Ly[module]/2.0;
		double ymax_strip = mod_Ly[module]/2.0;
		
		xmin_strip = (xmin_strip < -mod_Lx[module]/2.0) ? -mod_Lx[module]/2.0 : xmin_strip;
		xmax_strip = (xmax_strip > mod_Lx[module]/2.0) ? mod_Lx[module]/2.0 : xmax_strip;

		//These conditions shouldn't be able to be satisfied if vperphat.X == 0; thus avoiding divide-by-zero errors:
		if( xmin_strip == -mod_Lx[module]/2.0 ) ymin_strip = strip_center_pos.Y() + vperphat.Y()/vperphat.X() * ( xmin_strip - strip_center_pos.X() );
		if( xmax_strip == mod_Lx[module]/2.0 ) ymax_strip = strip_center_pos.Y() + vperphat.Y()/vperphat.X() * ( xmax_strip - strip_center_pos.X() );

		//These are local coordinates: need to convert to global coordinates; it should be sufficient to convert the two points:
		TVector3 point1local(xmin_strip, ymin_strip, 0.0);
		TVector3 point2local(xmax_strip,ymax_strip, 0.0 );
		TVector3 modcenter(mod_x0[module],mod_y0[module],mod_z0[module]);
		TVector3 point1global = mod_Rot[module] * point1local + modcenter;
		TVector3 point2global = mod_Rot[module] * point2local + modcenter;

		//How to implement a logarithmic color scale: ADCstrip over ADCmax varies from threshold/max up to 1:
		double logADCmin = log(thresh_stripsum/stripADCmax);
		double logADCmax = 0.0;

		int logADCbin = int( (log(ModData[module].ADCsum_ystrips[strip]/stripADCmax)-logADCmin)/(logADCmax-logADCmin)*double(ncolors) );

		//cout << "logADCmin, ADC, logADCbin = " << logADCmin << ", " << ADCsum_ystrips[module][strip] << ", " << logADCbin << ", color = " << gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,logADCbin))) << endl;
		int ADCbin = int( (ModData[module].ADCsum_ystrips[strip]-thresh_stripsum)/(stripADCmax-thresh_stripsum)*double(ncolors) );
		
		Ltemp.SetLineColor( gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,ADCbin))));
		
		//Ltemp.SetLineColor( gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,TMath::Nint(ADCsum_xstrips[module][strip]/stripADCmax*(ncolors-1))))));

		Ltemp.DrawLine( point1global.Y(), point1global.X(), point2global.Y(), point2global.X() );
		
	      } else { //strip IS along Y:
		double xmax_strip = strip_center_pos.X();
		double xmin_strip = strip_center_pos.X();
		double ymax_strip = mod_Ly[module]/2.0;
		double ymin_strip = -mod_Ly[module]/2.0;

		//These are local coordinates: need to convert to global coordinates; it should be sufficient to convert the two points:
		TVector3 point1local(xmin_strip, ymin_strip, 0.0);
		TVector3 point2local(xmax_strip,ymax_strip, 0.0 );
		TVector3 modcenter(mod_x0[module],mod_y0[module],mod_z0[module]);
		TVector3 point1global = mod_Rot[module] * point1local + modcenter;
		TVector3 point2global = mod_Rot[module] * point2local + modcenter;

		//How to implement a logarithmic color scale: ADCstrip over ADCmax varies from threshold/max up to 1:
		double logADCmin = log(thresh_stripsum/stripADCmax);
		double logADCmax = 0.0;

		int logADCbin = int( (log(ModData[module].ADCsum_ystrips[strip]/stripADCmax)-logADCmin)/(logADCmax-logADCmin)*double(ncolors) );

		//	cout << "logADCmin, ADC, logADCbin = " << logADCmin << ", " << ADCsum_ystrips[module][strip] << ", " << logADCbin << ", color = " << gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,logADCbin))) << endl;
		int ADCbin = int( (ModData[module].ADCsum_ystrips[strip]-thresh_stripsum)/(stripADCmax-thresh_stripsum)*double(ncolors) );
		
		Ltemp.SetLineColor( gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,ADCbin))));
		
		//Ltemp.SetLineColor( gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,TMath::Nint(ADCsum_xstrips[module][strip]/stripADCmax*(ncolors-1))))));

		Ltemp.DrawLine( point1global.Y(), point1global.X(), point2global.Y(), point2global.X() );
	      }
	      
	    }
	    
	  }
	}
	// hframe1->SetMinimum(0);
	// hframe2->SetMinimum(0);
	// hframe3->SetMinimum(0);
	// hframe4->SetMinimum(0);

	// hframe1->SetMaximum(20000);
	// hframe2->SetMaximum(20000);
	// hframe3->SetMaximum(20000);
	// hframe4->SetMaximum(20000);
	
	// c1->SetGrid();
	// c1->cd(1)->Clear();
	// hframe1->Draw("colz");
	// gPad->Modified();
	// gPad->Update();
	// c1->cd(2)->Clear();
	// hframe2->Draw("colz");
	// gPad->Modified();
	// gPad->Update();
	// c1->cd(3)->Clear();
	// hframe3->Draw("colz");
	// gPad->Modified();
	// gPad->Update();
	// c1->cd(4)->Clear();
	// hframe4->Draw("colz");
	gPad->Modified();
	gPad->Update();
	c1->Update();

	TMarker Mhit;
	TMarker Mtrack;
	//	TMarker Mhit_notrack;
	
	Mhit.SetMarkerStyle(5);
	Mhit.SetMarkerSize(3.0);
	Mhit.SetMarkerColor(6);
	//Mhit.SetLineColor(6);
	Mtrack.SetMarkerStyle(45);
	Mtrack.SetMarkerSize(3.0);
	Mtrack.SetMarkerColor(3);
	//Mtrack.SetLineColor(3);

	int mstyle_track[] = {45, 44, 43, 42, 49, 34, 28, 41};
	int mcolor_track[] = {1,2,3,4,5,6,7,8};
	
	for( set<int>::iterator imod=modules_hit.begin(); imod != modules_hit.end(); ++imod ){
	  //clusterdata_t clusttemp = mod_clusters[*imod];
	  int module = *imod;
	  int layer = mod_layer[module];
	  clusterdata_t clusttemp = mod_clusters[module];

	  for( int iclust=0; iclust<clusttemp.nclust2D; iclust++ ){
	    

	    //local hit position:
	    double xhit = clusttemp.xglobal2D[iclust];
	    double yhit = clusttemp.yglobal2D[iclust];
	    
	    

	    c1->cd(layer+1);

	    Mhit.DrawMarker( yhit, xhit );

	    gPad->Modified();
	    gPad->Update();
	    //gSystem->ProcessEvents();
	    c1->Update();
	  }
	}
	
	if( tracktemp.ntracks > 0 && tracktemp.ntracks<9){ //Draw track information for each layer:
	  for( int itrack=0; itrack<tracktemp.ntracks; itrack++ ){
	    for( int ihit=0; ihit<tracktemp.nhitsontrack[itrack]; ihit++ ){
	      int module=tracktemp.modlist_track[itrack][ihit];
	      int layer=mod_layer[module];

	      int iclust = tracktemp.hitlist_track[itrack][ihit];
	      
	      clusterdata_t clusttemp = mod_clusters[module];

	      
	      TVector3 hitpos_global(clusttemp.xglobal2D[iclust],clusttemp.yglobal2D[iclust],clusttemp.zglobal2D[iclust]);

	      //track position at module:
	      double xtrack = tracktemp.Xtrack[itrack]+tracktemp.Xptrack[itrack]*hitpos_global.Z();
	      double ytrack = tracktemp.Ytrack[itrack]+tracktemp.Yptrack[itrack]*hitpos_global.Z();

	      TVector3 trackpos_global( xtrack, ytrack, hitpos_global.Z() );
	     
	      
	      c1->cd(layer+1);

	      // Mhit.DrawMarker( yhit_local_stripcoord, xhit_local_stripcoord+nstripsx*(module%3) );

	      Mtrack.SetMarkerStyle( mstyle_track[itrack] );
	      Mtrack.SetMarkerColor( mcolor_track[itrack] );
	      Mtrack.DrawMarker( ytrack, xtrack );
	      //Mhit.Draw("SAME");
	      //Mtrack.Draw("SAME");
		  
	      gPad->Modified();
	      gPad->Update();
	      //gSystem->ProcessEvents();
	      c1->Update();
	      
	      
	    }
	  }
	  
	}
	gSystem->ProcessEvents();

	cout << "press any key to continue (q to quit):" << endl;
	TString reply;
	reply.ReadLine(cin,kFALSE);

	if( reply.BeginsWith("q") ) break;
	
      }
	  
	  
      //gSystem->Sleep(2500);

      //	  c1->Update();
	  
	  
    }
  

    for( int ilayer=0; ilayer<nlayers; ilayer++ ){
      // if( NstripX_layer[ilayer] > 0 ){
      hNclustX_layer->Fill( NclustX_layer[ilayer], ilayer );
	//}
      //if( NstripY_layer[ilayer] > 0 ){
      hNclustY_layer->Fill( NclustY_layer[ilayer], ilayer );
	//}
    }

    
      

    // if( eventdisplaymode != 0 ){

      
    // }
  }

  //form efficiency histograms:
  TClonesArray *heff_layer = new TClonesArray("TH2D",nlayers);
  
  for( int ilay=0; ilay<nlayers; ilay++ ){
    TString histname;

    TH2D *htemp = ( (TH2D*) (*hdidhit_layer)[ilay] );
    new( (*heff_layer)[ilay] ) TH2D( *htemp );

    ( (TH2D*) (*heff_layer)[ilay] )->SetName( histname.Format( "heff_layer%d", ilay ) );
    ( (TH2D*) (*heff_layer)[ilay] )->SetTitle( histname.Format( "Track-based efficiency, layer %d", ilay ) );
    
    ( (TH2D*) (*heff_layer)[ilay] )->Divide( (TH2D*) (*hshouldhit_layer)[ilay] );
  }	 
  
  //  hframe_layers->Delete();

  fout->Write();
  
  
}




// void find_clusters_by_module_old( moduledata_t mod_data, clusterdata_t &clust_data ){

//   filter_strips_by_module( mod_data, 0.3 );
  
//   //Build clusters from contiguous groups of x and y strips:

//   //Look for local maxima within contiguous groups of strips and fit template functions:
  
//   int module = mod_data.modindex;
//   int layer = mod_data.layerindex;
  
//   int nlocalmaxima=0;
//   int ncontiguous=0, laststrip=-1;
//   int striplo=-1, striphi=-1;
//   vector<int> stripmax;
//   vector<int> sampmax;
//   //vector<int> localmaxstrips;
  
//   //int maxstrip=-1;
//   double ADCmaxclust = 0.0;
//   int stripmaxclust = -1;
//   int sampmaxclust = -1;
  
  
//   int nclustx=0,nclusty=0;

//   //vectors of cluster properties:
//   vector<int> nstripxclust,nstripyclust,ixstriplo,ixstriphi,iystriplo,iystriphi,ixstripmax,iystripmax;
//   vector<double> xmean, xsigma, totalchargex, txmean, txsigma;
//   vector<double> ymean, ysigma, totalchargey, tymean, tysigma;

//   //  vector<int> nhitclust;
//   //"Hit" fit results:
//   int nhitx1D=0;
//   vector<int> clustidx_xhit;
//   vector<double> x0hit1D, tx0hit1D, Axhit1D, chi2ndf_hitX1D, sigx0hit1D, taux0hit1D;

//   int nhity1D=0;
//   vector<int> clustidx_yhit;
//   vector<double> y0hit1D, ty0hit1D, Ayhit1D, chi2ndf_hitY1D, sigy0hit1D, tauy0hit1D;

//   map<int,vector<double> > ADCsample_clust; //mapped by strip
  
//   double xsum=0.0, x2sum=0.0, tsum=0.0, t2sum=0.0, ADCsum=0.0;
//   double ysum=0.0, y2sum=0.0;

//   map<int,map<int,bool> > islocalmax; //mapped by strip/sample:
//   map<int,bool> isclustermax; 
  
//   for( set<int>::iterator ixstrip=mod_data.xstrips_filtered.begin(); ixstrip != mod_data.xstrips_filtered.end(); ++ixstrip ){
//     int strip=*ixstrip;

//     // if( fabs( mod_data.Tmean_xstrips[strip] - 64.3 )<=3.0*16.3 &&
//     // 	fabs( mod_data.Tsigma_xstrips[strip] - 36.3) <=3.0*5.8 ){
  
//     double xcoord_local = (strip+0.5 - 0.5*nstripsx)*0.4;
//     double ADCval = mod_data.ADCsum_xstrips[strip];
//     double tval  = mod_data.Tmean_xstrips[strip];

//     //  int stripleft=strip-1,stripright=strip+1;
//     // isclustermax[strip] = false;
   
//     if( laststrip < 0 || strip > laststrip+1 ||
// 	(strip-striplo+1 > maxnstripspercluster && strip-stripmaxclust > 1 ) ){ //new cluster:
//       if( laststrip >= 0 ){ //not the first cluster: add current cluster to the array before resetting sums:
// 	nstripxclust.push_back( ncontiguous );
// 	ixstriplo.push_back( striplo );
// 	ixstriphi.push_back( striphi );
// 	ixstripmax.push_back( stripmaxclust );
// 	xmean.push_back( xsum/ADCsum );
// 	xsigma.push_back( sqrt(fabs( x2sum/ADCsum - pow(xsum/ADCsum,2))) );
// 	totalchargex.push_back( ADCsum );
// 	txmean.push_back( tsum/ADCsum );
// 	txsigma.push_back( sqrt(fabs( t2sum/ADCsum - pow(tsum/ADCsum,2))));

// 	//set up and execute hit template fit:
// 	vector<double> xtemp,ytemp,ztemp;
// 	vector<double> extemp,eytemp,eztemp;
// 	//double Amaxtemp=0.0, tmaxtemp;
// 	map<int,double> Amaxstriptemp;
// 	map<int,double> tmaxstriptemp;

// 	//	  nlocalmaxima = 1;
	  
// 	nlocalmaxima=0;
	 
// 	//	  islocalmax.clear();
// 	stripmax.clear();
// 	sampmax.clear();

// 	// if( isclustermax[istripclust] ){
// 	stripmax.push_back( stripmaxclust );
// 	sampmax.push_back( sampmaxclust );
	 
// 	for( int istripclust=striplo; istripclust<=striphi; istripclust++ ){
		
// 	  //    Amaxstriptemp[istripclust] = 0.0;
// 	  // tmaxstriptemp[istripclust] = 0.0;

// 	  // if ANY local maximum occurs in this strip for any sample, increment the number of local maxima/hits:
// 	  //   if( islocalmax.find(istripclust) != islocalmax.end() ){
// 	  //  nlocalmaxima++;
// 	  //  stripmax.push_back( istripclust );
// 	  //}

	      
// 	  for( int isamp=0; isamp<6; isamp++ ){
// 	    xtemp.push_back( (istripclust+0.5 - 0.5*nstripsx)*0.4 ); //local X coordinate
// 	    ytemp.push_back( 12.5+25.0*isamp ); //center of time bin
// 	    ztemp.push_back( mod_data.ADCsamp_xstrips[istripclust][isamp] ); //ADC value for this sample:

// 	    //"Errors": 
// 	    extemp.push_back( 0.0 ); //0.25 * strip pitch (rough guess)
// 	    eytemp.push_back( 0.0 ); //0.25 * time bin width (rough guess)
// 	    eztemp.push_back( sigmasample ); //individual sample noise width (rough guess).

	      
	      
// 	    // if( mod_data.ADCsamp_xstrips[istripclust][isamp] > Amaxstriptemp[istripclust] ){
		
// 	    // 	Amaxstriptemp[istripclust] = mod_data.ADCsamp_xstrips[istripclust][isamp];
// 	    // 	tmaxstriptemp[istripclust] = 12.5+25.0*isamp;

// 	    // 	//if( Amaxstriptemp[istripclust] > ADCmaxtemp ) {
// 	    // 	//  ADCmaxtemp = Amaxstriptemp[istripclust];
// 	    // 	//  istripmax[nlocalmaxima-1] = istripclust;
// 	    // 	//}
// 	    // }

// 	    // if( mod_data.ADCsamp_xstrips[istripclust][isamp] > ADCmaxtemp ){
// 	    // 	ADCmaxtemp = mod_data.ADCsamp_xstrips[istripclust][isamp];
// 	    // 	stripmaxtemp = istripclust;
// 	    // }

// 	    double ADCleft = 0.0, ADCright=0.0, ADCnext=0.0, ADCprev=0.0;
// 	    if( istripclust > striplo ) ADCleft = mod_data.ADCsamp_xstrips[istripclust-1][isamp];
// 	    if( istripclust < striphi ) ADCright = mod_data.ADCsamp_xstrips[istripclust+1][isamp];
// 	    if( isamp > 0 ) ADCprev = mod_data.ADCsamp_xstrips[istripclust][isamp-1];
// 	    if( isamp < 5 ) ADCnext = mod_data.ADCsamp_xstrips[istripclust][isamp+1];
	      
// 	    //check for an additional local maximum in any of these samples with nsigma "prominence":
// 	    if( mod_data.ADCsamp_xstrips[istripclust][isamp] > ADCleft + sigmasample*localmaxthreshold_nsigma &&
// 		mod_data.ADCsamp_xstrips[istripclust][isamp] > ADCright + sigmasample*localmaxthreshold_nsigma &&
// 		mod_data.ADCsamp_xstrips[istripclust][isamp] > ADCprev + sigmasample*localmaxthreshold_nsigma &&
// 		mod_data.ADCsamp_xstrips[istripclust][isamp] > ADCnext + sigmasample*localmaxthreshold_nsigma &&
// 		( abs( istripclust - stripmaxclust ) > 1 || abs( isamp - sampmaxclust ) > 1 ) &&
// 		stripmax.size() < maxnhitspercluster ){
// 	      stripmax.push_back( istripclust );
// 	      sampmax.push_back( isamp );
// 	    }
	      
// 	  } 
// 	}
	  
// 	nlocalmaxima = stripmax.size();
	  
// 	TGraph2DErrors *gtemp = new TGraph2DErrors( xtemp.size(),
// 						    &(xtemp[0]), &(ytemp[0]), &(ztemp[0]),
// 						    &(extemp[0]), &(eytemp[0]), &(eztemp[0]) );

// 	int nfreeparamperhit = 5;
// 	if( varyclustersigma == 0 ) nfreeparamperhit--;
// 	if( varyclustertau   == 0 ) nfreeparamperhit--;
	  
// 	//We assume that the number of hits involved in the cluster is equal to the number of local maxima in the cluster occuring in unique strips:
// 	//We probably need to be smarter than this:
// 	int ndftemp = xtemp.size() - nfreeparamperhit*nlocalmaxima;

// 	if( ndftemp > 0 && fittemplatefunctions != 0 ){
// 	  HitFunc->FixParameter(0,TMath::Min(double(maxnhitspercluster),double(nlocalmaxima)));
// 	  for( int ihittemp=0; ihittemp<TMath::Min(nlocalmaxima,maxnhitspercluster); ihittemp++ ){
// 	    HitFunc->ReleaseParameter(5*ihittemp+1);
// 	    HitFunc->ReleaseParameter(5*ihittemp+2);
// 	    HitFunc->ReleaseParameter(5*ihittemp+3);
// 	    HitFunc->ReleaseParameter(5*ihittemp+4);
// 	    HitFunc->ReleaseParameter(5*ihittemp+5);

// 	    //	      for( int ihittemp=0; ihittemp<maxnhitspercluster; ihittemp++ ){
// 	    HitFunc->SetParLimits(5*ihittemp+1,0.0,100000.0);
// 	    HitFunc->SetParLimits(5*ihittemp+3,-300.0,450.0);
// 	    HitFunc->SetParLimits(5*ihittemp+4,0.0,100.0); //mm
// 	    HitFunc->SetParLimits(5*ihittemp+5,0.0,300.0);

	      
// 	    HitFunc->SetParameter(5*ihittemp+1, mod_data.ADCmax_xstrips[stripmax[ihittemp]] );
// 	    HitFunc->SetParameter(5*ihittemp+2, (stripmax[ihittemp]+0.5-0.5*nstripsx)*0.4 );
// 	    if( striplo == striphi ){ //only one strip; fix position and "sigma":
// 	      HitFunc->FixParameter(5*ihittemp+2, (stripmax[ihittemp]+0.5-0.5*nstripsx)*0.4 );
// 	      HitFunc->FixParameter(5*ihittemp+4, clustersigma );
// 	    }
// 	    HitFunc->SetParameter(5*ihittemp+3, 12.5+25.0*sampmax[ihittemp]);
// 	    HitFunc->SetParameter(5*ihittemp+4, clustersigma );
// 	    HitFunc->SetParameter(5*ihittemp+5, clustertau );
// 	    if( varyclustersigma == 0 ) HitFunc->FixParameter( 5*ihittemp+4,clustersigma );
// 	    if( varyclustertau   == 0 ) HitFunc->FixParameter( 5*ihittemp+5,clustertau   ); 
// 	  }

// 	  for( int ihittemp=nlocalmaxima; ihittemp<maxnhitspercluster; ihittemp++ ){
// 	    HitFunc->FixParameter( 5*ihittemp+1, 0.0 );
// 	    HitFunc->FixParameter( 5*ihittemp+2, 0.0 );
// 	    HitFunc->FixParameter( 5*ihittemp+3, 0.0 );
// 	    HitFunc->FixParameter( 5*ihittemp+4, 1.0 );
// 	    HitFunc->FixParameter( 5*ihittemp+5, 1.0 );
// 	  }

// 	  TFitResultPtr fitresult = gtemp->Fit( HitFunc, "S0Q", "" );

// 	  if( fitresult->IsValid() ){
	    
// 	    // cout << "successful fit to hit template function, nlocalmaxima = " << nlocalmaxima << endl;

// 	    //cout << "successful fit of hit template function to X strip data, nhits = " << nlocalmaxima << endl;


// 	    int ihitglobal=nhitx1D;
	      
// 	    nhitx1D += nlocalmaxima;
	      
// 	    for( int ihittemp=0; ihittemp<nlocalmaxima; ihittemp++ ){
// 	      x0hit1D.push_back( HitFunc->GetParameter( 5*ihittemp+2 ) );
// 	      tx0hit1D.push_back( HitFunc->GetParameter( 5*ihittemp+3 ) );
// 	      Axhit1D.push_back( HitFunc->GetParameter( 5*ihittemp+1 ) );
// 	      chi2ndf_hitX1D.push_back( fitresult->Chi2()/fitresult->Ndf() );
// 	      sigx0hit1D.push_back( HitFunc->GetParameter( 5*ihittemp+4 ) );
// 	      taux0hit1D.push_back( HitFunc->GetParameter( 5*ihittemp+5 ) );

// 	      clustidx_xhit.push_back( nclustx );
		
// 	      // cout << "strip lo, strip hi = " << striplo << ", " << striphi << endl;
// 	      // cout << "Hit " << ihittemp << ", (A,x0,t0,sigma,tau)=("
// 	      //      << Axhit1D[ihitglobal+ihittemp] << ", "
// 	      //      << x0hit1D[ihitglobal+ihittemp] << ", "
// 	      //      << tx0hit1D[ihitglobal+ihittemp] << ", "
// 	      //      << sigx0hit1D[ihitglobal+ihittemp] << ", "
// 	      //      << taux0hit1D[ihitglobal+ihittemp] << "), chi2/ndf = "
// 	      //      << chi2ndf_hitX1D[ihitglobal+ihittemp] << endl;
		  
// 	    }
	      
// 	  }

// 	}

// 	gtemp->Delete();
	  
// 	// cout << "Adding new x cluster, xmean = "
// 	//      << xmean[nclustx] << " mm, xsigma = "
// 	//      << xsigma[nclustx] << " mm, tmean = "
// 	//      << txmean[nclustx] << " ns, tsigma = "
// 	//      << txsigma[nclustx] << " ns, ADC sum = "
// 	//      << totalchargex[nclustx] << ", nstrips = "
// 	//      << nstripxclust[nclustx] << ", nlocalmaxima = "
// 	//      << nlocalmaxima << endl;

// 	//ADCsample_clust.clear();
// 	stripmax.clear();
// 	//islocalmax.clear();
	  
// 	nclustx++;
//       }

//       //copy ADC samples to the array for hit template fitting:
//       //ADCsample_clust[strip] = mod_data.ADCsamp_xstrips[strip];
	
//       ADCsum = ADCval;
//       xsum = xcoord_local*ADCval;
//       x2sum = pow(xcoord_local,2)*ADCval;
//       tsum = tval*ADCval;
//       t2sum = pow(tval,2)*ADCval;
	
//       ncontiguous = 1;

//       nlocalmaxima = 1;
//       striplo = strip;
//       striphi = strip;
//       //isclustermax[strip] = true;

//       stripmaxclust = strip;
//       ADCmaxclust = mod_data.ADCmax_xstrips[strip];
//       sampmaxclust = mod_data.isampmax_xstrips[strip];
//     } else { //strip = laststrip+1: Add to existing cluster:
	
//       ADCsum += ADCval;
//       xsum += xcoord_local*ADCval;
//       x2sum += pow(xcoord_local,2)*ADCval;
//       tsum += tval*ADCval;
//       t2sum += pow(tval,2)*ADCval;

//       // if( ADCval > ADCmax ){
//       //   if( strip > stripmax[nlocalmaxima-1]+1 ){
//       //     nlocalmaxima++;
//       //     stripmax.push_back( strip );
//       //   } else {
//       //     stripmax[nlocalmaxima-1] = strip;
//       //   }
//       //   ADCmax = ADCval;
//       // }

//       // if( ADCval > ADCmax ){
//       //   ADCmax = ADCval;
//       //   stripmax[nlocalmaxima-1] = strip;
//       // }

//       if( mod_data.ADCmax_xstrips[strip] > ADCmaxclust ){
// 	ADCmaxclust = mod_data.ADCmax_xstrips[strip];
// 	sampmaxclust = mod_data.isampmax_xstrips[strip];
// 	stripmaxclust = strip;
// 	//isclustermax[strip]=true;
// 	//isclustermax[laststrip] = false;
//       }
	
//       striphi=strip;
	
//       ncontiguous++;
//     }
//     laststrip = strip;
//   }
  

//   clust_data.nclustx = nclustx;
//   clust_data.nstripx = nstripxclust;
//   clust_data.ixstriplo = ixstriplo;
//   clust_data.ixstriphi = ixstriphi;
//   clust_data.ixstripmax = ixstripmax;
//   clust_data.xmean = xmean;
//   clust_data.xsigma = xsigma;
//   clust_data.totalchargex = totalchargex;
//   clust_data.txmean = txmean;
//   clust_data.txsigma = txsigma;

//   clust_data.nhitx1D = nhitx1D;
//   clust_data.clustidx_xhit = clustidx_xhit;
//   clust_data.x0hit1D = x0hit1D;
//   clust_data.tx0hit1D = tx0hit1D;
//   clust_data.Axhit1D = Axhit1D;
//   clust_data.chi2ndf_hitx1D = chi2ndf_hitX1D;
//   clust_data.sigx0hit1D = sigx0hit1D;
//   clust_data.taux0hit1D = taux0hit1D;
  
//   laststrip = -1;

//   ncontiguous = 0;

//   striplo=-1; striphi=-1;

//   ADCmaxclust = 0.0;
//   stripmaxclust = -1;
//   sampmaxclust = -1;
  
//   tsum=0.0;
//   t2sum=0.0;
//   ADCsum=0.0;
//   ysum=0.0;
//   y2sum=0.0;

//   islocalmax.clear();
  
//   for( set<int>::iterator iystrip=mod_data.ystrips_filtered.begin(); iystrip != mod_data.ystrips_filtered.end(); ++iystrip ){
//     int strip=*iystrip;

//     // if( fabs( mod_data.Tmean_ystrips[strip] - 64.3 )<=3.0*16.3 &&
//     // 	fabs( mod_data.Tsigma_ystrips[strip] - 36.3) <=3.0*5.8 ){

//     double ycoord_local = (strip+0.5 - 0.5*nstripsy)*0.4;
//     double ADCval = mod_data.ADCsum_ystrips[strip];
//     double tval  = mod_data.Tmean_ystrips[strip];

//     int stripleft = strip-1;
//     int stripright = strip+1;

//     for(int isamp=0; isamp<6; isamp++ ){
//       double ADCtemp = mod_data.ADCsamp_ystrips[strip][isamp];

//       double ADCleft=0.0, ADCright=0.0;
//       if( mod_data.ystrips.find(stripleft) != mod_data.ystrips.end() ){
// 	ADCleft = mod_data.ADCsamp_ystrips[stripleft][isamp];
//       }
//       if( mod_data.ystrips.find(stripright) != mod_data.ystrips.end() ){
// 	ADCright = mod_data.ADCsamp_ystrips[stripright][isamp];
//       }
	
//       if( ADCval > ADCleft && ADCval > ADCright ){
// 	islocalmax[strip][isamp] = true;
//       } else {
// 	islocalmax[strip][isamp] = false;
//       }
//     }
      
      
      
//     if( laststrip < 0 || strip > laststrip+1 ||
// 	(strip-striplo+1 > maxnstripspercluster && strip > stripmaxclust+1 ) ){ //new cluster:
//       if( laststrip >= 0 ){ //not the first cluster: add current cluster to the array before resetting sums:
// 	nstripyclust.push_back( ncontiguous );
// 	iystriplo.push_back( striplo );
// 	iystriphi.push_back( striphi );
// 	iystripmax.push_back( stripmaxclust );
// 	ymean.push_back( ysum/ADCsum );
// 	ysigma.push_back( sqrt(fabs( y2sum/ADCsum - pow(ysum/ADCsum,2))) );
// 	totalchargey.push_back( ADCsum );
// 	tymean.push_back( tsum/ADCsum );
// 	tysigma.push_back( sqrt(fabs( t2sum/ADCsum - pow(tsum/ADCsum,2))));

// 	// cout << "Adding new y cluster, ymean = "
// 	//      << ymean[nclusty] << " mm, ysigma = "
// 	//      << ysigma[nclusty] << " mm, tmean = "
// 	//      << tymean[nclusty] << " ns, tsigma = "
// 	//      << tysigma[nclusty] << " ns, ADC sum = "
// 	//      << totalchargey[nclusty] << ", nstrips = "
// 	//      << nstripyclust[nclusty] << endl;

// 	//set up and execute hit template fit:
// 	vector<double> xtemp,ytemp,ztemp;
// 	vector<double> extemp,eytemp,eztemp;
// 	//double Amaxtemp=0.0, tmaxtemp;
// 	map<int,double> Amaxstriptemp;
// 	map<int,double> tmaxstriptemp;

// 	//	  nlocalmaxima = 1;

// 	nlocalmaxima=0;
// 	stripmax.clear();
// 	sampmax.clear();

// 	stripmax.push_back( stripmaxclust );
// 	sampmax.push_back( sampmaxclust );
	  
// 	for( int istripclust=striplo; istripclust<=striphi; istripclust++ ){
	  
// 	  //ADCmaxtemp = 0.0;
// 	  for(int isamp=0; isamp<6; isamp++ ){    
// 	    xtemp.push_back( (istripclust+0.5 - 0.5*nstripsy)*0.4 );
// 	    ytemp.push_back( 12.0+25.0*isamp );
// 	    ztemp.push_back( mod_data.ADCsamp_ystrips[istripclust][isamp] );

// 	    extemp.push_back( 0.0 );
// 	    eytemp.push_back( 0.0 );
// 	    eztemp.push_back( sigmasample );

// 	    double ADCleft = 0.0, ADCright=0.0, ADCnext=0.0, ADCprev=0.0;
// 	    if( istripclust > striplo ) ADCleft = mod_data.ADCsamp_ystrips[istripclust-1][isamp];
// 	    if( istripclust < striphi ) ADCright = mod_data.ADCsamp_ystrips[istripclust+1][isamp];
// 	    if( isamp > 0 ) ADCprev = mod_data.ADCsamp_ystrips[istripclust][isamp-1];
// 	    if( isamp < 5 ) ADCnext = mod_data.ADCsamp_ystrips[istripclust][isamp+1];

// 	    //check for an additional local maximum in any of these samples with nsigma "prominence":
// 	    if( mod_data.ADCsamp_ystrips[istripclust][isamp] > ADCleft + sigmasample*localmaxthreshold_nsigma &&
// 		mod_data.ADCsamp_ystrips[istripclust][isamp] > ADCright + sigmasample*localmaxthreshold_nsigma &&
// 		mod_data.ADCsamp_ystrips[istripclust][isamp] > ADCprev + sigmasample*localmaxthreshold_nsigma &&
// 		mod_data.ADCsamp_ystrips[istripclust][isamp] > ADCnext + sigmasample*localmaxthreshold_nsigma &&
// 		( abs( istripclust - stripmaxclust ) > 1 || abs( isamp - sampmaxclust ) > 1 ) &&
// 		stripmax.size() < maxnhitspercluster ){
// 	      stripmax.push_back( istripclust );
// 	      sampmax.push_back( isamp );
// 	    }
	      
// 	  }
// 	}

// 	nlocalmaxima = stripmax.size();
	  
// 	TGraph2DErrors *gtemp = new TGraph2DErrors( xtemp.size(),
// 						    &(xtemp[0]), &(ytemp[0]), &(ztemp[0]),
// 						    &(extemp[0]), &(eytemp[0]), &(eztemp[0]) );

// 	int nfreeparamperhit = 5;
// 	if( varyclustersigma == 0 ) nfreeparamperhit--;
// 	if( varyclustertau   == 0 ) nfreeparamperhit--;
	  
// 	int ndftemp = xtemp.size() - nfreeparamperhit*nlocalmaxima;

// 	if( ndftemp > 0 && fittemplatefunctions != 0 ){ //then we can attempt a fit:
// 	  HitFunc->FixParameter(0,TMath::Min(double(maxnhitspercluster),double(nlocalmaxima)));
// 	  for( int ihittemp=0; ihittemp<TMath::Min(nlocalmaxima,maxnhitspercluster); ihittemp++ ){
// 	    HitFunc->ReleaseParameter(5*ihittemp+1);
// 	    HitFunc->ReleaseParameter(5*ihittemp+2);
// 	    HitFunc->ReleaseParameter(5*ihittemp+3);
// 	    HitFunc->ReleaseParameter(5*ihittemp+4);
// 	    HitFunc->ReleaseParameter(5*ihittemp+5);

// 	    HitFunc->SetParLimits(5*ihittemp+1,0.0,100000.0);
// 	    HitFunc->SetParLimits(5*ihittemp+3,-300.0,450.0);
// 	    HitFunc->SetParLimits(5*ihittemp+4,0.0,100.0); //mm
// 	    HitFunc->SetParLimits(5*ihittemp+5,0.0,300.0);
	      
// 	    HitFunc->SetParameter(5*ihittemp+1, mod_data.ADCmax_ystrips[stripmax[ihittemp]] );
// 	    HitFunc->SetParameter(5*ihittemp+2, (stripmax[ihittemp]+0.5-0.5*nstripsy)*0.4 );
// 	    if( striplo == striphi ){ //only one strip; fix position and "sigma":
// 	      HitFunc->FixParameter(5*ihittemp+2, (stripmax[ihittemp]+0.5-0.5*nstripsy)*0.4 );
// 	      HitFunc->FixParameter(5*ihittemp+4, clustersigma );
// 	    }
// 	    HitFunc->SetParameter(5*ihittemp+3, 12.5+25.0*sampmax[ihittemp]);
// 	    HitFunc->SetParameter(5*ihittemp+4, clustersigma );
// 	    HitFunc->SetParameter(5*ihittemp+5, clustertau );
// 	    if( varyclustersigma == 0 ) HitFunc->FixParameter( 5*ihittemp+4,clustersigma );
// 	    if( varyclustertau   == 0 ) HitFunc->FixParameter( 5*ihittemp+5,clustertau   ); 
// 	  }

// 	  for( int ihittemp=nlocalmaxima; ihittemp<maxnhitspercluster; ihittemp++ ){
// 	    HitFunc->FixParameter( 5*ihittemp+1, 0.0 );
// 	    HitFunc->FixParameter( 5*ihittemp+2, 0.0 );
// 	    HitFunc->FixParameter( 5*ihittemp+3, 0.0 );
// 	    HitFunc->FixParameter( 5*ihittemp+4, 1.0 );
// 	    HitFunc->FixParameter( 5*ihittemp+5, 1.0 );
// 	  }

// 	  TFitResultPtr fitresult = gtemp->Fit( HitFunc, "S0Q", "" );

// 	  if( fitresult->IsValid() ){
// 	    // cout << "successful fit of hit template function to Y strip data, nhits = " << nlocalmaxima << endl;


// 	    int ihitglobal=nhity1D;
	      
// 	    nhity1D += nlocalmaxima;
	      
// 	    for( int ihittemp=0; ihittemp<nlocalmaxima; ihittemp++ ){
// 	      y0hit1D.push_back( HitFunc->GetParameter( 5*ihittemp+2 ) );
// 	      ty0hit1D.push_back( HitFunc->GetParameter( 5*ihittemp+3 ) );
// 	      Ayhit1D.push_back( HitFunc->GetParameter( 5*ihittemp+1 ) );
// 	      chi2ndf_hitY1D.push_back( fitresult->Chi2()/fitresult->Ndf() );
// 	      sigy0hit1D.push_back( HitFunc->GetParameter( 5*ihittemp+4 ) );
// 	      tauy0hit1D.push_back( HitFunc->GetParameter( 5*ihittemp+5 ) );

// 	      clustidx_yhit.push_back( nclusty );
		
// 	      // cout << "strip lo, strip hi = " << striplo << ", " << striphi << endl;
// 	      // cout << "Hit " << ihittemp << ", (A,y0,t0,sigma,tau)=("
// 	      //      << Ayhit1D[ihitglobal+ihittemp] << ", "
// 	      //      << y0hit1D[ihitglobal+ihittemp] << ", "
// 	      //      << ty0hit1D[ihitglobal+ihittemp] << ", "
// 	      //      << sigy0hit1D[ihitglobal+ihittemp] << ", "
// 	      //      << tauy0hit1D[ihitglobal+ihittemp] << "), chi2/ndf = "
// 	      //      << chi2ndf_hitY1D[ihitglobal+ihittemp] << endl;
		  
// 	    }
	      
// 	  }
// 	}

// 	gtemp->Delete();

// 	// cout << "Adding new y cluster, ymean = "
// 	//      << ymean[nclusty] << " mm, ysigma = "
// 	//      << ysigma[nclusty] << " mm, tmean = "
// 	//      << tymean[nclusty] << " ns, tsigma = "
// 	//      << tysigma[nclusty] << " ns, ADC sum = "
// 	//      << totalchargey[nclusty] << ", nstrips = "
// 	//      << nstripyclust[nclusty] << ", nlocalmaxima = "
// 	//      << nlocalmaxima << endl;

// 	stripmax.clear();
// 	sampmax.clear();
	  
// 	nclusty++;
//       }

//       //stripmax.clear();
	
	
//       ADCsum = ADCval;
//       ysum = ycoord_local*ADCval;
//       y2sum = pow(ycoord_local,2)*ADCval;
//       tsum = tval*ADCval;
//       t2sum = pow(tval,2)*ADCval;
	
//       ncontiguous = 1;
//       striplo = strip;
//       striphi = strip;
//       //stripmax.push_back(strip);

//       stripmaxclust = strip;
//       ADCmaxclust = mod_data.ADCmax_ystrips[strip];
//       sampmaxclust = mod_data.isampmax_ystrips[strip];
	
//       //	ADCmax = ADCval;
//     } else { //strip = laststrip+1: Add to existing cluster:
//       ADCsum += ADCval;
//       ysum += ycoord_local*ADCval;
//       y2sum += pow(ycoord_local,2)*ADCval;
//       tsum += tval*ADCval;
//       t2sum += pow(tval,2)*ADCval;

//       if( mod_data.ADCmax_ystrips[strip] > ADCmaxclust ){
// 	ADCmaxclust = mod_data.ADCmax_ystrips[strip];
// 	sampmaxclust = mod_data.isampmax_ystrips[strip];
// 	stripmaxclust = strip;
// 	//isclustermax[strip]=true;
// 	//isclustermax[laststrip] = false;
//       }
	
//       striphi = strip;
	
//       ncontiguous++;
//     }
//     laststrip = strip;
//   }
  

//   // cout << "module " << module << " layer " << layer << " cluster statistics:" << endl;
//   // cout << "nclust x = " << nclustx << " nclust y = " << nclusty << endl;

//   clust_data.nclusty = nclusty;
//   clust_data.nstripy = nstripyclust;
//   clust_data.iystriplo = iystriplo;
//   clust_data.iystriphi = iystriphi;
//   clust_data.iystripmax = iystripmax;
//   clust_data.ymean = ymean;
//   clust_data.ysigma = ysigma;
//   clust_data.totalchargey = totalchargey;
//   clust_data.tymean = tymean;
//   clust_data.tysigma = tysigma;

//   clust_data.nhity1D = nhity1D;
//   clust_data.clustidx_yhit = clustidx_yhit;
//   clust_data.y0hit1D = y0hit1D;
//   clust_data.ty0hit1D = ty0hit1D;
//   clust_data.Ayhit1D = Ayhit1D;
//   clust_data.chi2ndf_hity1D = chi2ndf_hitY1D;
//   clust_data.sigy0hit1D = sigy0hit1D;
//   clust_data.tauy0hit1D = tauy0hit1D;

//   clust_data.nhit2D = 0;

//   clust_data.nclust2D = 0;
  
//   ////////////// 2D matching: try based on cluster information with no fit, and hit template fit results: /////////////////
//   ////now loop over X and Y clusters and try to match each x cluster with one y cluster:
//   //// Let's compute a correlation coefficient a la Evaristo:

  
  
//   if( clust_data.nclustx > 0 && clust_data.nclusty > 0 ){
//     clust_data.nclust2D = 0;
    
//     vector<bool> clustx_used(clust_data.nclustx);
//     vector<bool> clusty_used(clust_data.nclusty);

//     for( int ix=0; ix<clust_data.nclustx; ix++ ){
//       clustx_used[ix] = false;
//     }

//     for( int iy=0; iy<clust_data.nclusty; iy++ ){
//       clusty_used[iy] = false;
//     }

    
    
//     bool foundmatch = true;
//     while( foundmatch ){
//       foundmatch = false;
//       int ntest=0;
//       double bestchi2match=-1000.0;
//       int ixbest=-1, iybest=-1;

//       double bestmaxcor=-1.1;
      
//       for( int ix=0; ix<clust_data.nclustx; ix++ ){
// 	for( int iy=0; iy<clust_data.nclusty; iy++ ){
// 	  if( !clustx_used[ix] && !clusty_used[iy] ){

	   
	      
	    
	    
// 	    double ASYMtemp = (clust_data.totalchargex[ix]-clust_data.totalchargey[iy])/(clust_data.totalchargex[ix]+clust_data.totalchargey[iy]);
// 	    double chi2match = pow( ASYMtemp/0.1 ,2)+
// 	      pow( (clust_data.txmean[ix]-clust_data.tymean[iy])/5.1,2)+
// 	      pow( (clust_data.totalchargex[ix]-clust_data.totalchargey[iy])/20.0,2);

// 	    if( chi2match < bestchi2match || ntest == 0 ){
// 	      //if( ntest == 0 || cormax > bestmaxcor ){
// 	      bestchi2match = chi2match;
// 	      //bestmaxcor = cormax;
// 	      ixbest = ix;
// 	      iybest = iy;
// 	    }
// 	    ntest++;
// 	  }
// 	}
//       }
      
//       if( ntest > 0 && ixbest >= 0 && iybest >= 0 ){
// 	foundmatch = true;
// 	clust_data.nclust2D++;

// 	clust_data.ixclust2D.push_back( ixbest );
// 	clust_data.iyclust2D.push_back( iybest );
// 	clust_data.nstripx2D.push_back( clust_data.nstripx[ixbest] );
// 	clust_data.nstripy2D.push_back( clust_data.nstripy[iybest] );
// 	clust_data.xclust2D.push_back( clust_data.xmean[ixbest] );
// 	clust_data.yclust2D.push_back( clust_data.ymean[iybest] );
// 	clust_data.Eclust2D.push_back( 0.5*(clust_data.totalchargex[ixbest]+clust_data.totalchargey[iybest]) );
// 	clust_data.tclust2D.push_back( 0.5*(clust_data.txmean[ixbest]+clust_data.tymean[iybest]) );
// 	clust_data.dEclust2D.push_back( clust_data.totalchargex[ixbest]-clust_data.totalchargey[iybest] );
// 	clust_data.dtclust2D.push_back( clust_data.txmean[ixbest]-clust_data.tymean[iybest] );
// 	clust_data.CorrCoeff2D.push_back( bestmaxcor );

// 	clustx_used[ixbest] = true;
// 	clusty_used[iybest] = true;
//       }
//     }
    
//     if( clust_data.nhitx1D > 0 && clust_data.nhity1D > 0 ){

//       clust_data.nhit2D = 0;
//       //map<int, bool> hitx_used;
//       //map<int, bool> hity_used;
//       vector<bool> hitx_used(nhitx1D);
//       vector<bool> hity_used(nhity1D);

//       for( int ix=0; ix<nhitx1D; ix++ ){
// 	hitx_used[ix] = false;
//       }
//       for( int iy=0; iy<nhity1D; iy++ ){
// 	hity_used[iy] = false;
//       }

//       foundmatch = true;
    
//       while( foundmatch ){
// 	foundmatch = false;
// 	int ntest=0;
// 	double bestchi2match=-1000.0;

// 	int ixbestmatch=-1,iybestmatch=-1;
      
// 	for( int ix=0; ix<clust_data.nhitx1D; ix++ ){
// 	  for( int iy=0; iy<clust_data.nhity1D; iy++ ){

// 	    if( !hitx_used[ix] && !hity_used[iy] ){
// 	      double chi2match = pow( (tx0hit1D[ix]-ty0hit1D[iy])/sigmaT_2Dmatch, 2 ) +
// 		pow( (Axhit1D[ix]-Ayhit1D[iy])/sigmaADC_2Dmatch, 2 );
	    
// 	      if( chi2match < bestchi2match || ntest == 0 ){
// 		bestchi2match = chi2match;
// 		ixbestmatch=ix;
// 		iybestmatch=iy;
// 	      }
	    
// 	      ntest++;
// 	    }
// 	  }
// 	}

// 	if( ntest > 0 && ixbestmatch >= 0 && iybestmatch >= 0 ){
// 	  foundmatch = true;
// 	  clust_data.nhit2D++;

// 	  clust_data.xhit2D.push_back( x0hit1D[ixbestmatch] );
// 	  clust_data.yhit2D.push_back( y0hit1D[iybestmatch] );
// 	  clust_data.thit2D.push_back( 0.5*(tx0hit1D[ixbestmatch]+ty0hit1D[iybestmatch]) );
// 	  clust_data.Ahit2D.push_back( 0.5*(Axhit1D[ixbestmatch]+Ayhit1D[iybestmatch]) );
// 	  clust_data.dthit2D.push_back( tx0hit1D[ixbestmatch]-ty0hit1D[iybestmatch] );
// 	  clust_data.dAhit2D.push_back( Axhit1D[ixbestmatch]-Ayhit1D[iybestmatch] );
// 	  clust_data.ixhit2D.push_back( ixbestmatch );
// 	  clust_data.iyhit2D.push_back( iybestmatch );
// 	  clust_data.chi2match2D.push_back( bestchi2match );

// 	  hitx_used[ixbestmatch] = true;
// 	  hity_used[iybestmatch] = true;
// 	}
      
//       }
//     }
//   }
//   return;
// }
