#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include "TString.h"
//#include "TClonesArray.h"
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
#include "TStyle.h"
//#include "TApplication.h"

//contains struct and variables needed for the new track finder function --WX
#include "TrackFindingContainer.h"

//used for time measurement for the code --WX
#include <chrono> 
using namespace std::chrono;
auto totalTime = duration_cast<nanoseconds>(high_resolution_clock::now() - high_resolution_clock::now());

//TODO: don't hard-code the number of APV25 samples. But for this we need a better decoded hit format anyway

//Note that our definition of "X" and "Y" will differ from simprox.cpp method. What we are calling the "X" direction is the "vertical" (long) axis of the
//layer, with more strips (planeID == 1)
//What we are calling "Y" is the horizontal dimension with fewer strips (planeID == 0).

//We need to modify this so that we no longer assume that every module has the same numbers of strips or strip orientations:
//The only assumption we will retain is the assumption that each module has two non-parallel strip orientations, generically denoted "U" and "V"
//We will also make the strip pitch along each direction a configurable parameter:

const int MAXNCH = 60000;
const int MAXNFADC = 100;
const int MAXNTDC  = 100;

const int TOTAL_REQUIRED_HIT = 3;

int TrackingAlgorithmFlag = 0; //0 = OLD, "brute force" algorithm, 1 = NEW, "grid container" algorithm"

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

int maxneighborsX=4; //+/-4 around local maximum:
int maxneighborsY=3; //+/-3 around local maximum:

//double clustersigma=0.4; //mm
//double clustertau=50.0; //ns
//double sigmasample=20.0; //individual sample noise width
//double localmaxthreshold_nsigma=3.0; //number of sigmas to create a new local max/hit:
//double sigmaT_2Dmatch=26.0; //guess
//double sigmaADC_2Dmatch=18.0; //guess

TString DataFormat="INFN";

double cluster2Dmatch_asymcut=0.5; //dimensionless
double cluster2Dmatch_tcut = 25.0; //ns
double maxstripcorthreshold = 0.0;
//double maxADCXYthreshold = 3000.0; 
double stripcorthreshold = 0.0; //correlation coefficient must be positive:
double clustcorthreshold = 0.0;

double thresh_maxsample = 150.0; //max ADC sample on a strip must exceed this to be included
double thresh_stripsum =  300.0; //threshold on the sum of all samples on a strip:
double thresh_clustersum = 1000.0; //threshold on the summed ADCs of all strips in a cluster:

double sigma_sample = 20.0; //Sigma of ADC sample noise; used to get threshold for overlapping cluster condition
double thresh_2ndmax_nsigma = 5.0; //Number of sigmas above noise level to flag a second maximum in a contiguous grouping of strips. For this condition to occur, the 
double thresh_2ndmax_fraction = 0.25; //threshold to flag 2nd maximum as a fraction of first maximum
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

//flag to allow reuse of strips in multiple 2D clusters: this will require significant retooling of
//significant parts of the code
int reusestripsflag=0; 

//Default values for walk correction parameters:
//tstrip_cor = tstrip + pow( ADC0/ADCstrip, exp )
double walkcor_mean_const = 1.87334;
double walkcor_mean_ADC0=2818.0;
double walkcor_mean_exp =1.133;

double walkcor_sigma_const = 1.725;
double walkcor_sigma_ADC0=6448.0;
double walkcor_sigma_exp=0.78896;

double tau_pulseshape = 50.6; //ns 
//Width of strip timing cut for addition to cluster, in standard deviations:
double tstripcut_nsigma = 5.0;

//Walk correction for strip times:
TF1 *walkcor_mean_func = new TF1("walkcor_mean_func","[0]-pow([1]/x,[2])",250.0,2e4);
TF1 *walkcor_sigma_func = new TF1("walkcor_sigma_func","[0]+pow([1]/x,[2])",250.0,2e4);

//TF1 to fit strip times: Let's see how computationally expensive this is:
TF1 *PulseShape = new TF1("PulseShape","TMath::Max(0.0,[0]*(x-[1])*exp(-(x-[1])/[2]))",0.0,150.0);

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

  //Strip timing fit results:
  map<int,double> Tfit_xstrips;
  map<int,double> dTfit_xstrips;
  map<int,double> Afit_xstrips;
  map<int,double> dAfit_xstrips;
  map<int,double> taufit_xstrips;
  map<int,double> dtaufit_xstrips;
  
  map<int,double> Tfit_Chi2NDF_xstrips;

  map<int,double> Tfit_ystrips;
  map<int,double> dTfit_ystrips;
  map<int,double> Afit_ystrips;
  map<int,double> dAfit_ystrips;
  map<int,double> taufit_ystrips;
  map<int,double> dtaufit_ystrips;
  map<int,double> Tfit_Chi2NDF_ystrips;
  //map<int,double> Taufit_xstrips;
  //map<int,double> dTa

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

  int nkeep2D; 
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
  vector<vector<double> > ADCsamp_xclust;
  vector<vector<double> > ADCsamp_yclust;
  
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
//EXPERIMENTAL: Ratio of "y" gain to "x" gain: use to center ADC asymmetry at zero and improve correlation coefficient determination:
map<int,double> mod_RYX; //Ratio of Y gain to X gain: should set to 1 by default
map<int,double> mod_thresh_maxsample; //max sample threshold by module:
map<int,double> mod_thresh_stripsum;  //sum of ADC samples on strip
map<int,double> mod_thresh_clustersum; //cluster sum (technically sqrt(ADCX*ADCY) minimum value)
map<int,double> mod_stripcorthreshold; //XY correlation coefficient of max X and Y strips
map<int,double> mod_clustcorthreshold; //XY correlation coefficient of cluster-summed ADC samples
map<int,double> mod_ADCasymcut; //ADC asymmetry cut for cluster X and Y sums
map<int,double> mod_dTcut;      //tmeanx - tmeany cut for cluster summed ADC-weighted cluster mean times

//Also useful to define number of modules per layer and list of modules by layer:
map<int,int> nmodules_layer;
map<int,set<int> > modlist_layer;

// Tracking layer combinatorics: populate these arrays once so we don't do it every event:
//For each possible number of layers from 3 up to the total number of layers, we list all possible combinations of n layers
map<int,vector<vector<int> > > layercombos;


vector<double> zavg_layer;

vector<GridHitContainer> gridContainer;

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

      for( int isamp=0; isamp<nADCsamples; isamp++ ){
	double ADCx = mod_data.ADCsamp_xstrips[ixstrip][isamp]; //Here ADCX has already been multiplied by RYX which represents the gain ratio Y/X for the total cluster charge. But there is also a ratio of "strip gains": The average cluster charge is divided between about 4 X strips but only 3 Y strips. This means that for STRIP correlation coefficients, we should actually multiply the ADCx here by the ratio of average strip multiplicity between Y and X, before computing the correlation coefficient of the samples:
	double ADCy = mod_data.ADCsamp_ystrips[iystrip][isamp];
	
	sumx += ADCx;
	sumy += ADCy;
	sumx2 += pow(ADCx,2);
	sumy2 += pow(ADCy,2);
	sumxy += ADCx*ADCy;
      }

      double nSAMP = double(nADCsamples);
      
      double mux = sumx/nSAMP;
      double muy = sumy/nSAMP;
      double varx = sumx2/nSAMP-pow(mux,2);
      double vary = sumy2/nSAMP-pow(muy,2);
      double sigx = sqrt(varx);
      double sigy = sqrt(vary);

      double ccor = ( sumxy - nSAMP*mux*muy )/( nSAMP*sigx*sigy );

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

  clusttemp.nkeep2D = ngood;
}


int find_clusters_by_module_new( moduledata_t mod_data, clusterdata_t &clust_data ){
  int module = mod_data.modindex;
  int layer = mod_data.layerindex;
  int nclust=0;
  double cormax=-1.1;
  bool foundclust = true;

  clust_data.modindex=module;
  clust_data.layerindex=layer;
  
  map<int,bool> pixelused;
  map<int,double> pixelcorrcoeff;
  map<int,double> pixelADCXY; //sqrt(ADCX*ADCY)
  map<int,double> pixelADCX;
  map<int,double> pixelADCY;
  
  set<int> localmaxpixels; //ALL local maxima of sqrt(ADCY * ADCY)
  map<int,bool> is2x2; //Flag if local max is part of a cluster of at least 2x2 strips:
  map<int,bool> overlap; //Flag if cluster is part of a contiguous 2D grouping of strips overlapping with another local maximum
  map<int,set<int> > overlaplist; //List of all other local maxima with which this one overlaps; FIT to separate?
  
  bool any2x2 = false;
  //First loop: initialize pixelused flag to false AND calculate correlation coefficients of all pixels (should we do this?)
  for( set<int>::iterator ix=mod_data.xstrips.begin(); ix!=mod_data.xstrips.end(); ++ix ){
    int ixstrip = *ix;
    for( set<int>::iterator iy=mod_data.ystrips.begin(); iy!=mod_data.ystrips.end(); ++iy ){
      int iystrip = *iy;
      
      int pixel = ixstrip + mod_nstripsu[module]*iystrip;

      //  islocalmax[pixel] = false;
      //is2x2
      
      double sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;
	
      for( int isamp=0; isamp<nADCsamples; isamp++ ){
	double ADCx = mod_data.ADCsamp_xstrips[ixstrip][isamp];
	double ADCy = mod_data.ADCsamp_ystrips[iystrip][isamp];
	
	sumx += ADCx;
	sumy += ADCy;
	sumx2 += pow(ADCx,2);
	sumy2 += pow(ADCy,2);
	sumxy += ADCx*ADCy;
      }

      double nSAMP = double(nADCsamples);
      
      double mux = sumx/nSAMP;
      double muy = sumy/nSAMP;
      double varx = sumx2/nSAMP-pow(mux,2);
      double vary = sumy2/nSAMP-pow(muy,2);
      double sigx = sqrt(varx);
      double sigy = sqrt(vary);
      
      double ccor = ( sumxy - nSAMP*mux*muy )/( nSAMP*sigx*sigy );
 
      pixelused[pixel] = false;
      pixelcorrcoeff[pixel] = ccor;
      //Not clear if we actually want to calculate these other quantities yet:
      pixelADCXY[pixel] = sqrt(sumx*sumy);
      pixelADCX[pixel] = sumx;
      pixelADCY[pixel] = sumy;
    }
  }

  //Second loop: finding local maxima:
  for( set<int>::iterator ix=mod_data.xstrips.begin(); ix!=mod_data.xstrips.end(); ++ix ){
    int ixstrip = *ix;
    for( set<int>::iterator iy=mod_data.ystrips.begin(); iy!=mod_data.ystrips.end(); ++iy ){
      int iystrip = *iy;

      int pixel = ixstrip + mod_nstripsu[module]*iystrip;

      int ixlo=ixstrip, ixhi=ixstrip;
      int iylo=iystrip, iyhi=iystrip;

      double ADC = sqrt(mod_data.ADCsum_xstrips[ixstrip]*mod_data.ADCsum_ystrips[iystrip]);
      
      if( mod_data.xstrips.find(ixstrip+1) != mod_data.xstrips.end() ){
	ixhi = ixstrip+1;
      }
      if( mod_data.xstrips.find(ixstrip-1) != mod_data.xstrips.end() ){
	ixlo = ixstrip-1;
      }

      if( mod_data.ystrips.find(iystrip+1) != mod_data.ystrips.end() ){
	iyhi = iystrip+1;
      }

      if( mod_data.ystrips.find(iystrip-1) != mod_data.ystrips.end() ){
	iylo = iystrip-1;
      }

      bool islocalmax = true;
      // bool is2x2 = false;

      //Check all nearest-neighbor strips:
      for( int ixtest=ixlo; ixtest<=ixhi; ixtest++ ){
	for( int iytest=iylo; iytest<=iyhi; iytest++ ){
	  double ADCtemp = sqrt(mod_data.ADCsum_xstrips[ixtest]*mod_data.ADCsum_ystrips[iytest]);
	  if( ADCtemp > ADC && !(ixtest == ixstrip && iytest == iystrip ) ) islocalmax = false;
	}
      }
      
      if( islocalmax && pixelcorrcoeff[pixel] >= stripcorthreshold ){ //This pixel is a local maximum of sqrt(ADCX*ADCY) with a correlation coefficient above the threshold--evaluate as possible seed for clustering:

	localmaxpixels.insert(pixel);
	is2x2[pixel] = false;
	if( ixhi-ixlo+1 >= 2 && iyhi - iylo+1 >= 2 ){
	  is2x2[pixel] = true;
	  any2x2 = true;
	  //initialize "overlap" to false:
	  overlap[pixel] = false;
	}
	
	
      }
      
    }
  }

  //Third loop: create a map to flag overlapping clusters in X and/or Y: ask how many other local maxima does this cluster
  //share with a nearby local maximum: We have to avoid flagging fake hits as overlapping clusters. So what SHOULD we do?
  //If we have another cluster separated by two or more strips but less than 2*maxneighbors, then we check whether the two clusters
  // are contiguous; and if so, we "split" any shared strips in between two local maxima in proportion to the ratio of the two 
  // maxima and/or the distance from either maximum:
  // map<int,set<int> > overlapx; 
  // map<int,set<int> > overlapy; 

  // for( set<int>::iterator ipixel=localmaxpixels.begin(); ipixel!=localmaxpixels.end(); ++ipixel ){
  //   int ip=*ipixel;
  //   int ixp = ip % (mod_nstripsu[module]);
  //   int iyp = ip / (mod_nstripsu[module]);
  //   for( set<int>::iterator jpixel=localmaxpixels.begin(); jpixel!=localmaxpixels.end(); ++jpixel ){
  //     int jp = *jpixel;

  //     if( jp != ip ){

  // 	int jxp = jp % (mod_nstripsu[module]);
  // 	int jyp = jp / (mod_nstripsu[module]);
	
  // 	if( abs(jxp - ixp) > 1 && abs(jxp - ixp) < 2*maxneighborsX ){ // do these clusters overlap or potentially overlap in X?
  // 	  int ixlo=ixp, ixhi=jxp;
  // 	  if( ixlo > ixhi ){
  // 	    ixlo = jxp;
  // 	    ixhi = ixp; 
  // 	  }
  // 	  //Are these clusters contiguous? 
	  
  // 	  bool overlap = true;
  // 	  for( int ixtemp=ixlo+1; ixtemp<ixhi; ixtemp++ ){
  // 	    if( mod_data.xstrips.find(ixtemp) == mod_data.xstrips.end() ) overlap = false;
  // 	  }
	  
  // 	  if( overlap ){ //flag two clusters as overlapping in X
  // 	    overlapx[jp].insert(ip);
  // 	    overlapx[ip].insert(jp);
  // 	  }
  // 	}
	
  // 	if( abs(jyp-iyp) > 1 && abs(jyp-iyp) < 2*maxneighborsY ){ //do these clusters overlap or potentially overlap in Y?
  // 	  int iylo = iyp, iyhi=jyp;
  // 	  if( iylo > iyhi ){
  // 	    iylo = jyp;
  // 	    iyhi = iyp;
  // 	  }

  // 	  //is this pair of clusters contiguous in Y?
  // 	  bool overlap = true;
  // 	  for( int iytemp=iylo+1; iytemp<iyhi; iytemp++ ){
  // 	    if( mod_data.ystrips.find(iytemp) == mod_data.ystrips.end() ) overlap = false;
  // 	  }

  // 	  if( overlap ){ //flag two clusters as overlapping in Y
  // 	    overlapy[jp].insert(ip);
  // 	    overlapy[ip].insert(jp);
  // 	  }
  // 	}
	  
  //     }      
  //   }
  // }
  
  //Fourth loop: build candidate clusters around local maxima:
  for( set<int>::iterator ipixel=localmaxpixels.begin(); ipixel!=localmaxpixels.end(); ++ipixel ){

    int pixel = *ipixel;
    if( is2x2[pixel] || !any2x2 ){ //consider ONLY multi-strip clusters if there are ANY multistrip clusters:

      bool overlap = false;
      
      int ixpixel = pixel % (mod_nstripsu[module]);
      int iypixel = pixel / mod_nstripsu[module];
      set<int> xstriplist,ystriplist; //temporary lists of X and Y strips

      xstriplist.insert(ixpixel);
      ystriplist.insert(iypixel);

      //int ixtemp = ixpixel, iytemp=iypixel;
      int ixlo=ixpixel,ixhi=ixpixel,iylo=iypixel,iyhi=iypixel;

      while( mod_data.xstrips.find( ixlo-1 ) != mod_data.xstrips.end() && ixpixel-(ixlo-1) <=maxneighborsX
	     && xstriplist.size() < maxnstripXpercluster ){
	xstriplist.insert( ixlo-1 );
	ixlo--;
      }

      while( mod_data.xstrips.find( ixhi+1 ) != mod_data.xstrips.end() && ixhi+1 - ixpixel <= maxneighborsX
	     && xstriplist.size() < maxnstripXpercluster ){
	xstriplist.insert( ixhi+1 );
	ixhi++;
      }

      while( mod_data.ystrips.find( iylo-1 ) != mod_data.ystrips.end() && iypixel-(iylo-1) <= maxneighborsY
	     && ystriplist.size() < maxnstripYpercluster ){
	ystriplist.insert( iylo-1 );
	iylo--;
      }

      while( mod_data.ystrips.find( iyhi+1 ) != mod_data.ystrips.end() && iyhi+1 - iypixel <= maxneighborsY
	     && ystriplist.size() < maxnstripYpercluster ){
	ystriplist.insert( iyhi+1 );
	iyhi++;
      }

      //TODO: split overlapping clusters:
      //LATER

      double ADCsumy=0.0;
      double ADCsumx=0.0, xsum=0.0, ysum=0.0, xsum2=0.0, ysum2=0.0, txsum=0.0, txsum2=0.0,tysum=0.0,tysum2=0.0;      

      // double sumxsamp[nADCsamples],sumysamp[nADCsamples];

      //arrays to hold cluster-summed ADC samples:
      vector<double> xADCsamples(nADCsamples, 0.0);
      vector<double> yADCsamples(nADCsamples, 0.0);
      
      // for( int isamp=0; isamp<nADCsamples; isamp++ ){
      // 	sumxsamp[isamp] = 0.0;
      // 	sumysamp[isamp] = 0.0;
      // }

      double maxADCX=0.0;
      double xstripmax = -1.e9;
      int ixmax=-1;
      //Now compute candidate cluster properties: worry about how to deal with overlap later:
      for( int ix=ixlo; ix<=ixhi; ix++ ){
	double ADCtemp = mod_data.ADCsum_xstrips[ix];
	double xlocal = (ix + 0.5 - 0.5*mod_nstripsu[module])*mod_ustrip_pitch[module];
	//uncorrected strip time:
	double tstrip = mod_data.Tmean_xstrips[ix];

	if( ixmax < 0 || ADCtemp > maxADCX ){
	  maxADCX = ADCtemp;
	  ixmax = ix;
	  xstripmax = xlocal;
	}
	xsum += ADCtemp * xlocal;
	xsum2 += ADCtemp * pow(xlocal,2);
	txsum += ADCtemp * tstrip;
	txsum2 += ADCtemp * pow(tstrip,2);

	for( int isamp=0; isamp<nADCsamples; isamp++ ){
	  xADCsamples[isamp] += mod_data.ADCsamp_xstrips[ix][isamp];
	}

	ADCsumx += ADCtemp;
	
      }

      double maxADCY = 0.0;
      double ystripmax = -1.e9;
      int iymax=-1;
      for( int iy=iylo; iy<=iyhi; iy++ ){
	double ADCtemp = mod_data.ADCsum_ystrips[iy];
	double ylocal = (iy+0.5 - 0.5*mod_nstripsv[module])*mod_vstrip_pitch[module];
	double tstrip = mod_data.Tmean_ystrips[iy];

	if( iymax < 0 || ADCtemp > maxADCY ){
	  maxADCY = ADCtemp;
	  iymax = iy;
	  ystripmax = ylocal;
	}

	ysum += ADCtemp * ylocal;
	ysum2 += ADCtemp * pow(ylocal,2);
	tysum += ADCtemp * tstrip;
	tysum2 += ADCtemp * pow(tstrip,2);

	for( int isamp=0; isamp<nADCsamples; isamp++ ){
	  yADCsamples[isamp] += mod_data.ADCsamp_ystrips[iy][isamp];
	}

	ADCsumy += ADCtemp;
      }
      
      double csumxsamp =0.0, csumysamp=0.0, csumx2samp=0.0, csumy2samp=0.0, csumxysamp=0.0;
      for( int isamp=0; isamp<nADCsamples; isamp++ ){
	csumxsamp += xADCsamples[isamp];
	csumysamp += yADCsamples[isamp];
	csumx2samp += pow(xADCsamples[isamp],2);
	csumy2samp += pow(yADCsamples[isamp],2);
	csumxysamp += xADCsamples[isamp]*yADCsamples[isamp];
	// xADCsamples.push_back( sumxsamp[isamp] );
	// yADCsamples.push_back( sumysamp[isamp] );
      }

      double nSAMP = double(nADCsamples);
      
      double mux = csumxsamp/nSAMP;
      double muy = csumysamp/nSAMP;
      double varx = csumx2samp/nSAMP-pow(mux,2);
      double vary = csumy2samp/nSAMP-pow(muy,2);
      double sigx = sqrt(varx);
      double sigy = sqrt(vary);

      double ccor = (csumxysamp - nSAMP*mux*muy)/(nSAMP*sigx*sigy);

      int nstripx = ixhi - ixlo + 1;
      int nstripy = iyhi - iylo + 1;
      
      // if( sqrt(ADCsumx*ADCsumy) >= thresh_clustersum && ccor >= clustcorthreshold &&
      // 	  nstripx <= maxnstripXpercluster && nstripy <= maxnstripYpercluster ){
      if( nstripx <= maxnstripXpercluster && nstripy <= maxnstripYpercluster ){
	double txmean = txsum/ADCsumx;
	double tymean = tysum/ADCsumy;
	double txsigma = txsum2/ADCsumx - pow(txmean,2);
	double tysigma = tysum2/ADCsumy - pow(tymean,2);

	double xmean = xsum/ADCsumx;
	double ymean = ysum/ADCsumy;
	double xsigma = xsum2/ADCsumx - pow(xmean,2);
	double ysigma = ysum2/ADCsumy - pow(ymean,2);

	//1D cluster info:
	clust_data.nstripx.push_back( nstripx );
	clust_data.ixstriplo.push_back( ixlo );
	clust_data.ixstriphi.push_back( ixhi );
	clust_data.ixstripmax.push_back( ixmax );
	clust_data.nstripy.push_back( nstripy );
	clust_data.iystriplo.push_back( iylo );
	clust_data.iystriphi.push_back( iyhi );
	clust_data.iystripmax.push_back( iymax );
	clust_data.xmean.push_back( xmean );
	clust_data.ymean.push_back( ymean );
	clust_data.xsigma.push_back( xsigma );
	clust_data.ysigma.push_back( ysigma );
	clust_data.txmean.push_back( txmean );
	clust_data.tymean.push_back( tymean );
	clust_data.txsigma.push_back( txsigma );
	clust_data.tysigma.push_back( tysigma );
	clust_data.totalchargex.push_back( ADCsumx );
	clust_data.totalchargey.push_back( ADCsumy );
	
	//"2D" cluster info: 
	clust_data.itrack_clust2D.push_back( -1 );
	clust_data.ixclust2D.push_back( nclust );
	clust_data.iyclust2D.push_back( nclust );
	clust_data.nstripx2D.push_back( nstripx );
	clust_data.nstripy2D.push_back( nstripy );
	clust_data.xclust2D.push_back( xmean );
	clust_data.yclust2D.push_back( ymean );
	clust_data.Eclust2D.push_back( 0.5*(ADCsumx+ADCsumy) );
	clust_data.dEclust2D.push_back( ADCsumx - ADCsumy );
	clust_data.tclust2D.push_back( 0.5*(txmean+tymean) );
	clust_data.tclust2Dwalkcorr.push_back( 0.5*(txmean+tymean) );
	clust_data.dtclust2D.push_back( txmean-tymean );
	clust_data.dtclust2Dwalkcorr.push_back( txmean-tymean );
	clust_data.CorrCoeff2D.push_back( ccor );
	clust_data.CorrCoeffMaxStrips.push_back( pixelcorrcoeff[pixel] );
	clust_data.keepclust2D.push_back( true );
	clust_data.ADCsamp_xclust.push_back( xADCsamples );
	clust_data.ADCsamp_yclust.push_back( yADCsamples );
	clust_data.xmom2D.push_back( (xmean-xstripmax)/mod_ustrip_pitch[module] );
	clust_data.ymom2D.push_back( (ymean-ystripmax)/mod_vstrip_pitch[module] );
	clust_data.xclust2Dcorr.push_back( xmean );
	clust_data.yclust2Dcorr.push_back( ymean );

	double Utemp = xmean;
	double Vtemp = ymean;

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

      }
      
    }

  }

  clust_data.nclustx = nclust;
  clust_data.nclusty = nclust;
  clust_data.nclust2D = nclust;
  
  if( nclust > 0 ){
    prune_clusters( clust_data );
  }
  return 0;
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

  //Now we want to get a bit more intelligent: we want to find the local maximum in sqrt(ADCX*ADCY) with the largest correlation coefficient; preferably with a cluster
  //size of at least 2x2:
  map<int,bool> pixelused;
  map<int,bool> stripxused,stripyused;
  //map<int,bool> islocalmax; 
  //for( set<int>::iterator ix=mod_data.xstrips_filtered.begin(); ix!=mod_data.xstrips_filtered.end(); ++ix ){
  //Let's initialize this flag for ALL fired strips, regardless of whether they passed a correlation threshold:
  for( set<int>::iterator ix=mod_data.xstrips.begin(); ix != mod_data.xstrips.end(); ++ix ){
    int ixstrip = *ix;
    stripxused[ixstrip] = false; //these may not be used anymore:
    for( set<int>::iterator iy=mod_data.ystrips.begin(); iy != mod_data.ystrips.end(); ++iy ){
      int iystrip = *iy;
      stripyused[iystrip] = false;

      //The "pixel" array is nstripx*nstripy
      int pixel = ixstrip + mod_nstripsu[module]*iystrip;
      pixelused[pixel] = false;
    }
  }
  

  int returnval = 0;

  //A more intelligent strategy: find all local maxima of sqrt(ADCX*ADCY):
  
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

    bool any2x2 = false;
    
    for( set<int>::iterator ix=mod_data.xstrips_filtered.begin(); ix != mod_data.xstrips_filtered.end(); ++ix ){
      //for( set<int>::iterator iy=mod_data.ystrips_filtered.begin(); iy != mod_data.ystrips_filtered.end(); ++iy ){
      int ixstrip=*ix;
      for( set<int>::iterator iy=mod_data.ystrips_filtered.begin(); iy != mod_data.ystrips_filtered.end(); ++iy ){
	int iystrip = *iy;
	//int pixel = ixstrip+iystrip*nstripsy;
	int pixel = ixstrip + mod_nstripsu[module]*iystrip;
	
	//cluster size:
	if( (!stripxused[ixstrip] && !stripyused[iystrip]) ||
	    (!pixelused[pixel] && reusestripsflag != 0 ) ){ //find next local maximum of ADCX*ADCY
	
	  int sizextemp = 1, sizeytemp = 1;

	  bool is2x2 = false;
	  bool islocalmax = true;
	  //check if adjacent strips fired:
	  //int ixlo=ixstrip, ixhi = ixstrip;
	  //int iylo=iystrip, iyhi = iystrip;

	  int ixlo=ixstrip, ixhi = ixstrip;
	  int iylo=iystrip, iyhi = iystrip;
	
	  double ADCXY = sqrt(mod_data.ADCsum_xstrips[ixstrip]*mod_data.ADCsum_ystrips[iystrip]);
	  //Check nearest-neighbor strips +/-1 in X and Y:
	  if( mod_data.xstrips.find( ixstrip + 1 ) != mod_data.xstrips.end() ){
	    if( !stripxused[ixstrip+1] || reusestripsflag != 0 ){
	      ixhi = ixstrip+1;
	    }
	  }
	  if( mod_data.xstrips.find( ixstrip - 1 ) != mod_data.xstrips.end() ){
	    if( !stripxused[ixstrip-1] || reusestripsflag != 0 ){
	      ixlo = ixstrip-1;
	    }
	  }
	  if( mod_data.ystrips.find( iystrip + 1 ) != mod_data.ystrips.end() ){
	    if( !stripyused[iystrip+1] || reusestripsflag != 0 ){
	      iyhi = iystrip+1;
	    }
	  }
	  if( mod_data.ystrips.find( iystrip - 1 ) != mod_data.ystrips.end() ){
	    if( !stripyused[iystrip-1] || reusestripsflag != 0 ){
	      iylo = iystrip-1;
	    }
	  }

	  for( int ixtest=ixlo; ixtest<=ixhi; ixtest++ ){
	    for( int iytest=iylo; iytest<=iyhi; iytest++ ){
	      double ADCXYtemp = sqrt(mod_data.ADCsum_xstrips[ixtest]*mod_data.ADCsum_ystrips[iytest]);
	      if( ADCXYtemp > ADCXY && !(ixtest==ixstrip&&iytest==iystrip) ) islocalmax = false;
	    }
	  }
	
	  if( ixhi-ixlo+1 > 1 && iyhi-iylo + 1 > 1 && islocalmax ){
	    is2x2 = true; 
	  }
	  //	double ADCXY = mod_data.ADCmax_xstrips[ixstrip] * mod_data.ADCmax_ystrips[iystrip]; 

	  //If there are ANY other 2x2 clusters found and this is NOT one of them, set islocalmax to false
	  if( !is2x2 && any2x2 ) islocalmax = false;
	
	  //Question: do we really need to compute these again? I think maybe not, but to be on the safe side let's do it anyway:
	
	  double sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;
	
	  for( int isamp=0; isamp<nADCsamples; isamp++ ){
	    double ADCx = mod_data.ADCsamp_xstrips[ixstrip][isamp];
	    double ADCy = mod_data.ADCsamp_ystrips[iystrip][isamp];
	  
	    sumx += ADCx;
	    sumy += ADCy;
	    sumx2 += pow(ADCx,2);
	    sumy2 += pow(ADCy,2);
	    sumxy += ADCx*ADCy;
	  }

	  double nSAMP = double(nADCsamples);
	
	  double mux = sumx/nSAMP;
	  double muy = sumy/nSAMP;
	  double varx = sumx2/nSAMP-pow(mux,2);
	  double vary = sumy2/nSAMP-pow(muy,2);
	  double sigx = sqrt(varx);
	  double sigy = sqrt(vary);
	
	  double ccor = ( sumxy - nSAMP*mux*muy )/( nSAMP*sigx*sigy );

	  //double ADCXY = sqrt( mod_data.ADCsum_xstrips[ixstrip]*mod_data.ADCsum_ystrips[iystrip] );

	  // double txcorr = mod_data.Tmean_xstrips[ixstrip]+walkcor_mean_func->Eval( mod_data.ADCsum_xstrips[ixstrip] );
	  // 	double tycorr = mod_data.Tmean_ystrips[iystrip]+walkcor_mean_func->Eval( mod_data.ADCsum_ystrips[iystrip] );
	
	  // //      double tavgcorr = 0.5*(txcorr+tycorr);  
	  // 	double tsigmaxcorr = walkcor_sigma_func->Eval( mod_data.ADCsum_xstrips[ixstrip] );
	  // 	double tsigmaycorr = walkcor_sigma_func->Eval( mod_data.ADCsum_ystrips[iystrip] );

	  double dtcorr = mod_data.Tmean_xstrips_walkcor[ixstrip]-mod_data.Tmean_ystrips_walkcor[iystrip];
	  double sigdt = sqrt(pow(mod_data.Tsigma_xstrips_walkcor[ixstrip],2)+pow(mod_data.Tsigma_ystrips_walkcor[iystrip],2));
	
	  //only seed a cluster from a local maximum if the correlation coefficient is reasonable:
	  if( ((!stripxused[ixstrip] && !stripyused[iystrip])||(reusestripsflag != 0 && !pixelused[pixel]) )&& ccor >= stripcorthreshold && islocalmax ){
	      // fabs( dtcorr ) <= tstripcut_nsigma*sigdt && islocalmax ){
	    if( ixmax < 0 || ccor > maxcor || (is2x2 && !any2x2) ){
	      //the third condition here overrides any previous pixel choice if this happens to be the first local maximum
	      // in a cluster of at least 2x2
	      //if( ixmax < 0 || ADCXY > maxADCXY ){
	      maxcor = ccor;
	      maxADCXY = ADCXY;
	      ixmax = ixstrip;
	      iymax = iystrip;

	      if( is2x2 ) any2x2 = true;
	    }
	  }
	}
      }
    }

    //From here on out the rest is the same:
    
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

      pixelused[ixmax+mod_nstripsu[module]*iymax] = true; 
      
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
	  
	  
	  if( (!stripxused[ixleft] || reusestripsflag != 0 ) && abs( bestymatch - iymax ) < 2000 &&
	      fabs( txcorr - tcorr_mean ) <= tstripcut_nsigma*sigtxcorr &&
	      xstriplistclust.find(ixleft) == xstriplistclust.end() ){
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
	  
	  if( (!stripxused[ixright] || reusestripsflag != 0 ) && abs( bestymatch - iymax ) < 2000 &&
	      fabs( txcorr - tcorr_mean ) <= tstripcut_nsigma*sigtxcorr &&
	      xstriplistclust.find(ixright) == xstriplistclust.end() ){
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
	  
	  if( (!stripyused[iyleft] || reusestripsflag != 0 ) && abs( bestxmatch - ixmax ) < 2000 &&
	      fabs( tycorr - tcorr_mean ) <= tstripcut_nsigma*sigtycorr &&
	      ystriplistclust.find(iyleft) == ystriplistclust.end() ){
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
	  
	  if( (!stripyused[iyright] || reusestripsflag != 0 ) && abs( bestxmatch - ixmax ) < 2000 &&
	      fabs( tycorr - tcorr_mean ) <= tstripcut_nsigma*sigtycorr &&
	      ystriplistclust.find(iyright) == ystriplistclust.end() ){
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
      double sumxsamp[nADCsamples],sumysamp[nADCsamples];

      for( int isamp=0; isamp<nADCsamples; isamp++ ){
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

	for( int isamp=0; isamp<nADCsamples; isamp++ ){
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
	for( int isamp=0; isamp<nADCsamples; isamp++ ){
	  sumysamp[isamp] += mod_data.ADCsamp_ystrips[iystrip][isamp];
	}
      }

      vector<double> xADCsamples,yADCsamples;
      
      double csumxsamp =0.0, csumysamp=0.0, csumx2samp=0.0, csumy2samp=0.0, csumxysamp=0.0;
      for( int isamp=0; isamp<nADCsamples; isamp++ ){
	csumxsamp += sumxsamp[isamp];
	csumysamp += sumysamp[isamp];
	csumx2samp += pow(sumxsamp[isamp],2);
	csumy2samp += pow(sumysamp[isamp],2);
	csumxysamp += sumxsamp[isamp]*sumysamp[isamp];
	xADCsamples.push_back( sumxsamp[isamp] );
	yADCsamples.push_back( sumysamp[isamp] );
      }

      double nSAMP = double(nADCsamples);
      
      double mux = csumxsamp/nSAMP;
      double muy = csumysamp/nSAMP;
      double varx = csumx2samp/nSAMP-pow(mux,2);
      double vary = csumy2samp/nSAMP-pow(muy,2);
      double sigx = sqrt(varx);
      double sigy = sqrt(vary);

      double ccor = (csumxysamp - nSAMP*mux*muy)/(nSAMP*sigx*sigy);

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

void new_find_tracks(map<int,clusterdata_t> mod_clusters, trackdata_t &trackdata, TVector3 fcp, TVector3 bcp)
{
    //a new track finding algorithm for faster execution speed
    
    //compute the constraint slope for the track we are looking for
    double constraintSlopeX =  (fcp.X() - bcp.X())/(fcp.Z() - bcp.Z());
    double constraintSlopeY =  (fcp.Y() - bcp.Y())/(fcp.Z() - bcp.Z());
    
    //copy the hit into the new container that has grid like boxes 
    //for faster hit search
    for( map<int,clusterdata_t>::iterator imod=mod_clusters.begin(); imod != mod_clusters.end(); ++imod ){ 
        int module = imod->first;
        int layer = mod_layer[module];
        
        if (layer >= gridContainer.size()){
            cout<<"no layer "<<layer<<" stored in the grid container"<<endl;
            continue;
        }
        
        if (!gridContainer[layer].HasModule(module)){
            cout<<"grid container "<<layer<<" did not register module "<<module<<endl;
            continue;
        }
        for( int iclust=0; iclust<mod_clusters[module].nclust2D; iclust++ ){
        gridContainer[layer].AddHit(layer, module, iclust, 
                                    mod_clusters[module].xglobal2D[iclust],
                                    mod_clusters[module].yglobal2D[iclust],
                                    mod_clusters[module].zglobal2D[iclust],
                                    sigma_hitpos, sigma_hitpos); 
        }
        //TODO, the resolution should be more precise based on the GEM clustering
    }
    
    //now we figure out what layer combination we need to look at. We only
    //look at combination that has different end layers. In the case of choosing 
    //at lest 4 out of 6, we have 22 different combinations, but only 6 of them 
    //have different end layers: 03, 04, 05, 14, 15, 25.
    
    //when we look at configuration 05, we can allow at most two missing hits
    //when we look at configuration 03, then there will be no missing hit allowed
    //since it is assumed that layer 4 and 5 have no hits. Otherwise the track should 
    //already be picked up from configuration 05
    
    //TODO this part should be move out of this function, appearently there is no
    //reason to calculate this configuration for every event
    vector<pair<int, int>> layerconfig;
    layerconfig.clear();
    //cout<<"size of grid container "<<gridContainer.size()<<endl;
    for (int i = gridContainer.size(); i >= TOTAL_REQUIRED_HIT; i--){
    
        map<int,vector<vector<int> > >::iterator it = layercombos.find(i);
        
        if (it == layercombos.end()){
            cout<<"layercombos does not contain any combo with required hit "<<i<<endl;
            continue;
        }
        
        for (unsigned int j=0; j<layercombos[i].size(); j++){
        
            if (layercombos[i][j].size() == 0) continue;
            
            int front = layercombos[i][j][0];
            int back  = layercombos[i][j][layercombos[i][j].size() - 1];
            assert(front < back);
            bool exist = false;
            
            for (unsigned int k=0; k<layerconfig.size(); k++){
                if (front == layerconfig[k].first && back == layerconfig[k].second)
                exist = true;
            }
            
            if (!exist) layerconfig.push_back(pair<int, int> (front, back));
        }
        
    }
    //NOTE: the above procedure should ensure that layerconfig is sorted, that is
    //we start with configuration that can potentially give tracks with more hits
    
    //the container for all the recon tracks
    vector<SBSReconTrack> allTracks;
    allTracks.clear();
    
    long long int totalCombo = 1;
    
    for (unsigned int conf = 0; conf < layerconfig.size(); conf++){
        int front = layerconfig[conf].first;
        int back  = layerconfig[conf].second;
        
        //calculate the total number of possible combo, skip this config if 
        //there are too many
        totalCombo = 1;
        for (int i=front; i<=back; i++) totalCombo *= gridContainer[i].GetTotalHit();
        if (totalCombo > maxnhitcombinations) continue;
        
        //we start tracking with the front and back trackers (seed), this provides
        //the largest leverage so that we have a good estimate for the track slope
        //at the beginning
        for (std::map<int, std::vector<GridHit*>>::iterator itf = gridContainer[front].hitArray.begin(); 
             itf != gridContainer[front].hitArray.end(); ++itf){
            for (std::map<int, std::vector<GridHit*>>::iterator itb = gridContainer[back].hitArray.begin(); 
                 itb != gridContainer[back].hitArray.end(); ++itb){
                int modf = itf->first;
                int modb = itb->first;
                
                for (unsigned int i = 0; i<gridContainer[front].hitArray[modf].size(); i++){
                    //pass if the hit is used
                    if (gridContainer[front].hitArray[modf][i]->used) continue;
                    for (unsigned int j=0; j<gridContainer[back].hitArray[modb].size(); j++){
                        //pass if the hit is used
                        if (gridContainer[back].hitArray[modb][j]->used) continue;
                        
                        int nmissinghits = gridContainer.size() - (back - front + 1);
                        
                        SBSReconTrack thisTrack(nmissinghits);
                        thisTrack.AddAndFilterWithLLS(gridContainer[back].hitArray[modb][j] );
                        thisTrack.AddAndFilterWithLLS(gridContainer[front].hitArray[modf][i]);
                        
                        if (fabs(thisTrack.trackPara[1] - constraintSlopeX) > TrackMaxSlopeX || 
                            fabs(thisTrack.trackPara[3] - constraintSlopeY) > TrackMaxSlopeY) continue;
                            
                        //calculate the projection on the front and back constaint plane
                        //reject if not satisfy the cuts
                        
                        double projx, projy;
                        thisTrack.GetProjection(fcp.Z(), projx, projy);
                        
                        if (fabs(projx - fcp.X()) > TrackProjPosCut[0] || fabs(projy - fcp.Y()) > TrackProjPosCut[1])
                        continue;
                        
                        thisTrack.GetProjection(bcp.Z(), projx, projy);
                        
                        if (fabs(projx - bcp.X()) > TrackProjPosCut[0] || fabs(projy - bcp.Y()) > TrackProjPosCut[1])
                        continue;
                        
                        //now if the pair passed the selection rules, we will construct a seed
                        
                        vector<SBSReconTrack> tracksystem;
                        tracksystem.push_back(thisTrack);
                        
                        //and then we use the middle layers to build the track
                        
                        assert(back - front > 1); //at lest one middle layer in between
                        for (int mid = front + 1; mid < back; mid++){
                            int ntrack = tracksystem.size();
                            for (unsigned int itrack = 0; itrack < ntrack; itrack++){
                                
                                if (!tracksystem[itrack].status) continue;
                                
                                map<int, vector<vector<GridHit*>*> > hitinrange;
                                int nhitsInRange = 0;
                                
                                for (std::map<int, std::vector<GridHit*>>::iterator itm = gridContainer[mid].hitArray.begin(); 
                                     itm != gridContainer[mid].hitArray.end(); ++itm){
                                     
                                     int modm = itm->first;
                                     double projx, projy;
                                     
                                     tracksystem[itrack].GetProjection(mod_z0[modm], projx, projy);
                                     
                                     int nhitsInMod = 0;
                                     hitinrange[modm] = gridContainer[mid].GetHitInRange(projx, projy, 
                                                        TrackFindingMaxRadius, modm, nhitsInMod);
                                     
                                     nhitsInRange += nhitsInMod;
                                }
                                
                                if (nhitsInRange == 0){
                                    //no hits found within the range
                                    tracksystem[itrack].nmissinghits++;
                                }
                                else{
                                    //hits found in the range, for the first hit
                                    //add it to the original track, there there is more 
                                    //than one hit, make a copy of the track and add it to the 
                                    //new track
                                    bool noHitFound = true;
                                    SBSReconTrack tmpTrack(&tracksystem[itrack]);
                                    map<int, vector<vector<GridHit*>*> >::iterator ii = hitinrange.begin();
                                    for (; ii != hitinrange.end(); ++ii){
                                        int thisMod = ii->first;
                                        for (unsigned int jj=0; jj<hitinrange[thisMod].size(); jj++){
                                            for (unsigned int kk=0; kk<hitinrange[thisMod][jj]->size(); kk++){
                                                if ((hitinrange[thisMod][jj]->at(kk))->used) continue;
                                                
                                                double projx, projy;
                                                tmpTrack.GetProjection((hitinrange[thisMod][jj]->at(kk))->z, projx, projy);
                                                double dist = sqrt(pow(projx - (hitinrange[thisMod][jj]->at(kk))->x, 2) 
                                                                 + pow(projy - (hitinrange[thisMod][jj]->at(kk))->y, 2));
                                                
                                                if (dist > TrackFindingMaxRadius) continue;
                                                
                                                double dz = hitinrange[thisMod][jj]->at(kk)->z - tmpTrack.hits.back()->z;
                                                assert(dz > 0);
                                                double tmpSlopeX = (hitinrange[thisMod][jj]->at(kk)->x - tmpTrack.hits.back()->x)/dz;
                                                double tmpSlopeY = (hitinrange[thisMod][jj]->at(kk)->y - tmpTrack.hits.back()->y)/dz;
                                                
                                                if (fabs(tmpSlopeX - constraintSlopeX) > TrackMaxSlopeX || 
                                                    fabs(tmpSlopeY - constraintSlopeY) > TrackMaxSlopeY) continue;
                                                
                                                if (noHitFound) {
                                                    tracksystem[itrack].AddAndFilterWithLLS(hitinrange[thisMod][jj]->at(kk));
                                                    noHitFound = false;
                                                    
                                                }
                                                else{
                                                    //made a copy of the track first
                                                    tracksystem.push_back(tmpTrack);
                                                    tracksystem[tracksystem.size() - 1].AddAndFilterWithLLS(hitinrange[thisMod][jj]->at(kk));
                                                }
                                            }
                                        }
                                    }
                                    if (noHitFound) tracksystem[itrack].nmissinghits++;
                                    else{
                                        //even if there is hit found we still consider no hit found as one possiblity
                                        tracksystem.push_back(tmpTrack);
                                        tracksystem[tracksystem.size() - 1].nmissinghits++;
                                    }
                                }
                            }
                            //examine all track status before moving towards the next plane
                            for (unsigned int itrack = 0; itrack < tracksystem.size(); itrack++){
                                if (!tracksystem[itrack].status) continue;
                                if (tracksystem[itrack].nmissinghits > gridContainer.size() - TOTAL_REQUIRED_HIT )
                                tracksystem[itrack].status = false; 
                            }
                        }
                        
                        //finish track following with the middle layers, now we check if 
                        //the tracks are good
                        for (unsigned int itrack = 0; itrack < tracksystem.size(); itrack++){
                            //check minimum hit requirement
                            if (tracksystem[itrack].nmissinghits > gridContainer.size() - TOTAL_REQUIRED_HIT ){
                                tracksystem[itrack].status = false; 
                                continue;
                            }
                            //check track slope, can use smaller cuts after fitting all hits
                            if (fabs(tracksystem[itrack].trackPara[1] - constraintSlopeX) > FineTrackMaxSlope[0] || 
                                fabs(tracksystem[itrack].trackPara[3] - constraintSlopeY) > FineTrackMaxSlope[1]){
                                tracksystem[itrack].status = false;
                                continue;
                            }
                            //check projection, can use smaller cuts after fitting all hits
                            double projx, projy;         
                            tracksystem[itrack].GetProjection(fcp.Z(), projx, projy);
                            if (fabs(projx - fcp.X()) > FineTrackProjPosCut[0] || fabs(projy - fcp.Y()) > FineTrackProjPosCut[1]){
                                tracksystem[itrack].status = false;
                                continue;
                            }
                            tracksystem[itrack].GetProjection(bcp.Z(), projx, projy);
                            if (fabs(projx - bcp.X()) > FineTrackProjPosCut[0] || fabs(projy - bcp.Y()) > FineTrackProjPosCut[1]){
                                tracksystem[itrack].status = false;
                                continue;
                            }
                            
                            std::sort(tracksystem[itrack].hits.begin(), tracksystem[itrack].hits.end(), SortHits);
                            tracksystem[itrack].CheckOutlier(TrackFindingMaxRadius);
                            if (tracksystem[itrack].nmissinghits > gridContainer.size() - TOTAL_REQUIRED_HIT ){
                                tracksystem[itrack].status = false; 
                                continue;
                            }
                            
                            if (tracksystem[itrack].ForwardSearchTest(TrackFindingMaxRadius)){
                                tracksystem[itrack].status = false; 
                                continue;
                            }
                            
                            //if the track passed above cuts, calculate chi2/ndf
                            //this is a rather simplified chi2 calculation (assume xy readout and no rotation)
                            //let's use the calculation in the original standalone code
                            //tracksystem[itrack].CalChi2NDF();
                            
                            //the original calculation
                            int ndf = 2*tracksystem[itrack].hits.size()-4;
                            double chi2 = 0.;
                            for (unsigned int ihit = 0; ihit < tracksystem[itrack].hits.size(); ihit++){
                                int module = tracksystem[itrack].hits[ihit]->module;
                                int iclust = tracksystem[itrack].hits[ihit]->index;
                                
                                double zhittemp = mod_clusters[module].zglobal2D[iclust]; 
		                        double uhittemp = mod_clusters[module].xclust2Dcorr[iclust];
		                        double vhittemp = mod_clusters[module].yclust2Dcorr[iclust];

		                        TVector3 trackpos_global( tracksystem[itrack].trackPara[0]+tracksystem[itrack].trackPara[1]*zhittemp, 
		                                                  tracksystem[itrack].trackPara[2]+tracksystem[itrack].trackPara[3]*zhittemp, 
		                                                  zhittemp );
		                        TVector3 modcenter_global( mod_x0[module], mod_y0[module], mod_z0[module] );
		                       
		                        
		                        TVector3 trackpos_local = mod_Rotinv[module]*(trackpos_global - modcenter_global);

		                        //Compute the local "u" and "v" coordinates of the track within the module, to allow for the
		                        //possibility of different strip orientations than X/Y:
		                        double utracktemp = trackpos_local.X()*mod_Pxu[module] + trackpos_local.Y()*mod_Pyu[module];
		                        double vtracktemp = trackpos_local.X()*mod_Pxv[module] + trackpos_local.Y()*mod_Pyv[module];
		                        
		                        chi2 += pow( (uhittemp-utracktemp)/sigma_hitpos, 2 ) + pow( (vhittemp - vtracktemp)/sigma_hitpos, 2 );
                            }
                            tracksystem[itrack].chi2ndf = chi2/ndf;
                            
                            if (tracksystem[itrack].chi2ndf > TrackChi2Cut) tracksystem[itrack].status = false;
                        }
                        
                        //save the track into the track container of the event and set all its hits as used
                        for (unsigned int itrack = 0; itrack < tracksystem.size(); itrack++){
                            if (tracksystem[itrack].status)
                            allTracks.push_back(tracksystem[itrack]);
                        }
                        
                    }
                }
            }
        }
    }

    //finished track finding, now we save the tracks into output container
    std::sort(allTracks.begin(), allTracks.end(), SortTracks);
    
    map<int, vector<GridHit*>> goodHits;
    
    for (unsigned int i=0; i<allTracks.size(); i++){
        
        bool flag = false;
        for (unsigned int j=0; j<allTracks[i].hits.size(); j++){
            map< int, vector<GridHit*> >::iterator it = goodHits.find(allTracks[i].hits[j]->layer);
            if (it != goodHits.end()){
                for (UInt_t n = 0; n<(it->second).size(); n++){
                    if ((allTracks[i].hits[j]->layer  == ((it->second).at(n))->layer ) &&
                        (allTracks[i].hits[j]->module == ((it->second).at(n))->module) &&
                        (allTracks[i].hits[j]->index  == ((it->second).at(n))->index )) { flag = true; }
                }
            }
        }
        if (flag) continue;
        
        
        vector<int> modlisttemp,hitlisttemp;
	    vector<double> residxtemp,residytemp;

	    for (unsigned int j=0; j<allTracks[i].hits.size(); j++){
	        
	        goodHits[allTracks[i].hits[j]->layer].push_back(allTracks[i].hits[j]);
	        
	        modlisttemp.push_back(allTracks[i].hits[j]->module);
	        hitlisttemp.push_back(allTracks[i].hits[j]->index);
	        
	        //TODO: really shouldn't perform the same calculation twice here
            double zhittemp = mod_clusters[modlisttemp.back()].zglobal2D[hitlisttemp.back()];
            double uhittemp = mod_clusters[modlisttemp.back()].xclust2Dcorr[hitlisttemp.back()];
            double vhittemp = mod_clusters[modlisttemp.back()].yclust2Dcorr[hitlisttemp.back()];

            TVector3 trackpos_global( allTracks[i].trackPara[0]+allTracks[i].trackPara[1]*zhittemp, 
                                      allTracks[i].trackPara[2]+allTracks[i].trackPara[3]*zhittemp, 
                                      zhittemp );
            TVector3 modcenter_global( mod_x0[modlisttemp.back()], mod_y0[modlisttemp.back()], mod_z0[modlisttemp.back()] );
           
            
            TVector3 trackpos_local = mod_Rotinv[modlisttemp.back()]*(trackpos_global - modcenter_global);

            //Compute the local "u" and "v" coordinates of the track within the module, to allow for the
            //possibility of different strip orientations than X/Y:
            double utracktemp = trackpos_local.X()*mod_Pxu[modlisttemp.back()] + trackpos_local.Y()*mod_Pyu[modlisttemp.back()];
            double vtracktemp = trackpos_local.X()*mod_Pxv[modlisttemp.back()] + trackpos_local.Y()*mod_Pyv[modlisttemp.back()];
            
            residxtemp.push_back(uhittemp - utracktemp);
		    residytemp.push_back(vhittemp - vtracktemp);
	    }
	    
	    trackdata.nhitsontrack.push_back( allTracks[i].hits.size() );
	    trackdata.modlist_track.push_back( modlisttemp );
	    trackdata.hitlist_track.push_back( hitlisttemp );
	    trackdata.residx_hits.push_back( residxtemp );
	    trackdata.residy_hits.push_back( residytemp );
	    trackdata.eresidx_hits.push_back( residxtemp );
	    trackdata.eresidy_hits.push_back( residytemp );
	    
	    trackdata.Xtrack.push_back( allTracks[i].trackPara[0] );
	    trackdata.Xptrack.push_back( allTracks[i].trackPara[1] );
	    trackdata.Ytrack.push_back( allTracks[i].trackPara[2] );
	    trackdata.Yptrack.push_back( allTracks[i].trackPara[3] );
	    
	    trackdata.Chi2NDFtrack.push_back( allTracks[i].chi2ndf );
	    
	    trackdata.ntracks++;
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

    //prune_clusters( mod_clusters[module] ); this was moved to find_clusters
    
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

  if( layers_2Dmatch.size() >= TOTAL_REQUIRED_HIT ){
    //trackdata_t trackdatatemp;

    trackdata.ntracks = 0;
    
    bool foundtrack=true;
    while( foundtrack ){ //consider all possible combinations of one (2D matched) hit per layer:
      // This "brute force" approach will work reasonably well for cosmic data, but we'll need to be smarter when dealing with real data
      //under high-rate conditions
      //set<int> layers_2Dmatch; //list of layers with at least one 2D-matched hit:
      foundtrack = false;

      int nhitsrequired=layers_2Dmatch.size(); //on first iteration, require nhits = total number of layers with unused hits:

      while( nhitsrequired >= TOTAL_REQUIRED_HIT ){ //first iteration: determine number of layers with (unused) hits, populate lists of unused hits, count number of combinations,
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

	  //onbesttrack.clear();
	  
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

	      //clear certain arrays so they don't screw up chi^2 calculation:
	      ontrack.clear();
	      hitcombo.clear();
	      xresid_layer.clear();
	      yresid_layer.clear();
	      
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
		
		//		for(set<int>::iterator ilayer=layerswithfreehits.begin(); ilayer!=layerswithfreehits.end(); ++ilayer){
		for( set<int>::iterator ilayer = layerstotest.begin(); ilayer != layerstotest.end(); ++ilayer ){
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

int get_nearest_module( trackdata_t trdata, int ilayer, int itrack=0 ){
  int bestmod = -1;
  double r2min = 0.0;
  
  if( itrack >= 0 && itrack < trdata.ntracks && ilayer >= 0 && ilayer < nlayers ){
    //for( int ilayer=0; ilayer<nlayers; ilayer++ ){

    //Check if this layer is on the fitted track. If so, then there is no ambiguity:
    for( int ihit=0; ihit<trdata.nhitsontrack[itrack]; ihit++ ){
      int module = trdata.modlist_track[itrack][ihit];
      int layer = mod_layer[module];
      if( layer == ilayer ) return module;
    }
    
    //If we make it to this point, then ilayer is NOT on the track: find closest module to the track in the layer in question:
    for( int imodule=0; imodule<nmodules; imodule++ ){
      if( mod_layer[imodule] == ilayer ){ //then module is in this tracking layer:
	double xtrtemp = trdata.Xtrack[itrack] + mod_z0[imodule]*trdata.Xptrack[itrack];
	double ytrtemp = trdata.Ytrack[itrack] + mod_z0[imodule]*trdata.Yptrack[itrack];
	
	double r2mod = pow(xtrtemp-mod_x0[imodule],2)+pow(ytrtemp-mod_y0[imodule],2);
	
	if( bestmod == -1 || r2mod < r2min ){
	  bestmod = imodule;
	  r2min = r2mod;
	}
      }
    }
  }
  
  return bestmod; 
}

void event_display(){

  const char *filename = "/home/sbs-onl/analysis/long_cosmic_runs/jan_12_2021/gem_hit_2811_12.root";
  const char *configfilename = "./configs/config_UVA_EEL_5layer_Jan2021.txt";
  const TString outfilename = "2811_event_display";
  int nevents_save = 100;

  auto program_start = high_resolution_clock::now(); //starting time of the program --WX
  //Initialize walk correction parameters:
  double walkcor_mean_params[3] = {walkcor_mean_const, walkcor_mean_ADC0, walkcor_mean_exp};
  double walkcor_sigma_params[3] = {walkcor_sigma_const, walkcor_sigma_ADC0, walkcor_mean_exp};
  walkcor_mean_func->SetParameters(walkcor_mean_params);
  walkcor_sigma_func->SetParameters(walkcor_sigma_params);
  
  gROOT->ProcessLine(".x ~/rootlogon.C");
  gStyle->SetPalette(kRainBow);
  
  gStyle->SetOptStat(0);
  

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

	  if( skey == "TrackingAlgorithm" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    TrackingAlgorithmFlag = stemp.Atoi();
	  }
	  
	  if( skey == "reusestripsflag" ){ //flag to allow reuse of strips in multiple 2D clusters:
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    reusestripsflag = stemp.Atoi();
	  }
	  
	  if( skey == "taupulseshape" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    tau_pulseshape = stemp.Atof();
	  }
	  
	  if( skey == "dataformat" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    DataFormat = stemp;
	  }
	  
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

	  if( skey == "nADCsamples" ){ //ordinarily, this will always be 6, but we don't want to hardcode.
	    TString snsamp =  ( (TObjString*) (*tokens)[1] )->GetString();
	    nADCsamples = snsamp.Atoi();
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

	  if( skey == "mod_RYX" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      mod_RYX[i-1] = stemp.Atof();
	    }
	  }
	  
	  if( skey == "eventdisplay" && ntokens >= 2 ){
	    TString sevdisplay = ( (TObjString*) (*tokens)[1] )->GetString();

	    //eventdisplaymode = sevdisplay.Atoi();
	    eventdisplaymode = 1;
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

	    
	    usehitmaps=true;
	  }
	      
	}
	tokens->Delete();
      }
    }
  }


   //initialize the grid hit container here -- WX
  gridContainer.clear();
  
  nmodules_layer.clear();
  modlist_layer.clear();
  
  zavg_layer.resize(nlayers);
  int nmod_layer[nlayers];
  for( int ilayer=0; ilayer<nlayers; ilayer++ ){
    nmod_layer[ilayer] = 0;
    zavg_layer[ilayer] = 0.0;
    gridContainer.push_back(GridHitContainer(ilayer));

    nmodules_layer[ilayer] = 0;
    modlist_layer[ilayer].clear();
  }
  for( int imod=0; imod<nmodules; imod++ ){
    int layer = mod_layer[imod];
    nmod_layer[layer]++;
    zavg_layer[layer] += mod_z0[imod];
    for (unsigned int j=0; j<gridContainer.size(); j++){
        if (gridContainer[j].layer == layer) gridContainer[j].AddModule(imod);
    }

    nmodules_layer[layer]++;
    modlist_layer[layer].insert(imod);
    
  }

  int maxnmodperlayer=-1;
  
  for( int ilayer=0; ilayer<nlayers; ilayer++ ){
    zavg_layer[ilayer] /= double(nmod_layer[ilayer]);
    cout << "ilayer, zavg = " << ilayer << ", " << zavg_layer[ilayer] << endl;

    if( nmodules_layer[ilayer] > maxnmodperlayer ) maxnmodperlayer = nmodules_layer[ilayer];
  }

  cout << "max number modules per layer = " << maxnmodperlayer << endl;

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

    if( nlayersoncombo >= TOTAL_REQUIRED_HIT ){
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
  double CALOsum;
  int NGOODSCINT;
  int EventID;
  //int NstripX[nmodules];
  //int NstripY[nmodules];
  

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

    //Set default value in case the user hasn't defined something:
    if( mod_RYX.find( imodule ) == mod_RYX.end() ) mod_RYX[imodule] = 1.0;
  }
  
  
  TCanvas *c1 = new TCanvas("c1","c1",1600,1500);
  c1->Divide( maxnmodperlayer,nlayers,.001,.001);
  
  TCanvas *c_proj = new TCanvas("c1_proj","c1_proj",1200,800);
  c_proj->Divide(2,1);
  
  TClonesArray *hframe_modules = new TClonesArray("TH2D",nmodules);
  
    
  
  TChain *C = new TChain("GEMHit");

  C->Add(filename);

  //GEMHit_tree_HallAtest *T = new GEMHit_tree_HallAtest(C);
  //Instead of using a skeleton class, let's define all needed variables and initialize using SetBranchAddress depending on configuration:

  //Variables common to all three formats (so far):
  Int_t evtID; //event ID 
  Int_t nch;   //number of channels:
  Int_t Strip[MAXNCH]; //strip index within module:
  Int_t moduleID[MAXNCH]; //Unique identifier of module
  Int_t planeID[MAXNCH];  //Layer ID
  Int_t axis[MAXNCH];     //strip orientation flag;
  //ADC samples: 
  Int_t adc0[MAXNCH];
  Int_t adc1[MAXNCH];
  Int_t adc2[MAXNCH];
  Int_t adc3[MAXNCH];
  Int_t adc4[MAXNCH];
  Int_t adc5[MAXNCH];
  //Calorimeter info: specific to data format "HALLA":
  Int_t nfadc;
  Int_t fCH[MAXNFADC];
  Int_t ftimting[MAXNFADC];
  Int_t fadc[MAXNFADC];
  Int_t ntdc;
  Int_t tCH[MAXNTDC];
  Double_t ttiming[MAXNTDC];

  //Set Branch addresses common to all three data formats:
  C->SetBranchAddress("evtID",&evtID);
  C->SetBranchAddress("nch",&nch);
  C->SetBranchAddress("strip",Strip);
  C->SetBranchAddress("adc0",adc0);
  C->SetBranchAddress("adc1",adc1);
  C->SetBranchAddress("adc2",adc2);
  C->SetBranchAddress("adc3",adc3);
  C->SetBranchAddress("adc4",adc4);
  C->SetBranchAddress("adc5",adc5);
  
  if( DataFormat == "INFN" || DataFormat == "HALLA" || DataFormat == "GEP_MC"){
    C->SetBranchAddress("detID",moduleID);
    C->SetBranchAddress("planeID",axis);
    //No layer ID for INFN or Hall A data
    if( DataFormat == "HALLA" ){
      C->SetBranchAddress("nfadc",&nfadc);
      C->SetBranchAddress("ntdc",&ntdc);
      C->SetBranchAddress("fCH",fCH);
      C->SetBranchAddress("ftimting",ftimting);
      C->SetBranchAddress("fadc",fadc);
      C->SetBranchAddress("tCH",tCH);
      C->SetBranchAddress("ttiming",ttiming);
    }
  }

  if( DataFormat == "UVA" ){
    C->SetBranchAddress("planeID",planeID);
    C->SetBranchAddress("moduleID",moduleID);
    C->SetBranchAddress("axis",axis);
  }
  
  
  cout << "Total events = " << C->GetEntries() << endl;
  
  long nevent=0;

  
  TClonesArray *hxyhit_layer = new TClonesArray("TH2D",nlayers);
  
  long ntotal = C->GetEntries();

  //guesstimate binning for efficiency histos so that we can have ~few hundred events/bin
  double nbins_eff = double( ntotal )/400.0;

  //3.75*n^2 = ntot

  
    
  for( int ilayer=0; ilayer<nlayers; ilayer++ ){

    double Lx_layer = xgmax_layer[ilayer] - xgmin_layer[ilayer];
    double Ly_layer = ygmax_layer[ilayer] - ygmin_layer[ilayer]; 

    double nbinsy_eff,nbinsx_eff;


    int nbinsaxis_hitmap = 250;
    int nbinsx_hitmap,nbinsy_hitmap;
   
    //These formulas give:
    // nbins_tot = nbinsx*nbinsy = sqrt(nbins_eff/ratio)*ratio*sqrt(nbins_eff/ratio) = nbins_eff
  
    if( Lx_layer > Ly_layer ){
      double ratio = Lx_layer/Ly_layer;
      nbinsy_eff = sqrt( nbins_eff/ratio );
      nbinsx_eff = ratio*nbinsy_eff;

      nbinsy_hitmap = nbinsaxis_hitmap;
      nbinsx_hitmap = TMath::Nint(ratio*nbinsaxis_hitmap);
    
    } else {
      double ratio = Ly_layer/Lx_layer;
      nbinsx_eff = sqrt( nbins_eff/ratio );
      nbinsy_eff = ratio*nbinsx_eff;

      nbinsx_hitmap = nbinsaxis_hitmap;
      nbinsy_hitmap = TMath::Nint(ratio*nbinsaxis_hitmap);
    
    }
    
   
  }
    map<int,bool> mod_coord_flag;
    
  for( int imod=0; imod<nmodules; imod++ ){
    TString hnametemp;
    int layer = mod_layer[imod];


    //For event display, the long dimension of the layer is plotted horizontally, this means that we want to check coordinate system:

    double Lx_layer = xgmax_layer[layer] - xgmin_layer[layer];
    double Ly_layer = ygmax_layer[layer] - ygmin_layer[layer]; 

    if( Ly_layer < Lx_layer ){ //Y coordinate is "short" dimension of layer, plot Y (X) coordinate on the Y (X) axis (UVA-style)
      new( (*hframe_modules)[imod] ) TH2D( hnametemp.Format("hframe_layer%d_module%d",layer,imod), hnametemp.Format("Layer %d, Module %d", layer, imod),
					   200, -mod_Lx[imod]/2.0-15.0,mod_Lx[imod]/2.0+15.0,
					   200, -mod_Ly[imod]/2.0-15.0,mod_Ly[imod]/2.0+15.0 );

      ( (TH2D*) (*hframe_modules)[imod] )->SetXTitle("X local (mm)");
      ( (TH2D*) (*hframe_modules)[imod] )->SetYTitle("Y local (mm)");
      mod_coord_flag[imod] = true;
    } else { //X coordinate is "short" dimension of layer, plot X (Y) coordinate on Y (X) axis (INFN-style)
      new( (*hframe_modules)[imod] ) TH2D( hnametemp.Format("hframe_layer%d_module%d",layer,imod), hnametemp.Format("Layer %d, Module %d", layer, imod),
					   200, -mod_Ly[imod]/2.0-15.0,mod_Ly[imod]/2.0+15.0,
					   200, -mod_Lx[imod]/2.0-15.0,mod_Lx[imod]/2.0+15.0 );

      ( (TH2D*) (*hframe_modules)[imod] )->SetXTitle("Y local (mm)");
      ( (TH2D*) (*hframe_modules)[imod] )->SetYTitle("X local (mm)");
      mod_coord_flag[imod] = false;
    }
  }
    

  double zgmin_global = 10000;
  double zgmax_global = -10000;
  double xgmin_global = 10000;
  double xgmax_global = -10000;
  double ygmin_global = 10000;
  double ygmax_global = -10000;
  int i_min_x = 0;
  int i_max_x = 0;
  int i_min_y = 0;
  int i_max_y = 0;

	
  for(int i=0; i<nmodules; i++){
    
    if(mod_z0[i] < zgmin_global) zgmin_global = mod_z0[i];
    if(mod_z0[i] > zgmax_global) zgmax_global = mod_z0[i];
    
    if(mod_x0[i] < xgmin_global){
      xgmin_global = mod_x0[i];
      i_min_x = i;
    }
    if(mod_x0[i] > xgmax_global){
      xgmax_global = mod_x0[i];
      i_max_x = i;
    }

    if(mod_y0[i] < ygmin_global){
      ygmin_global = mod_y0[i];
      i_min_y = i;
    }
    if(mod_y0[i] > ygmax_global){
      ygmax_global = mod_y0[i];
      i_max_y = i;
    }
  }
  
  TH2D *proj_xz = new TH2D("proj_xz","z-x Projection;x global (cm);z global (cm)",200, (xgmin_global-mod_Lx[i_min_x])/10, (xgmax_global+mod_Lx[i_max_x])/10, 200, (zgmin_global - 100)/10, (zgmax_global + 100)/10);
  proj_xz->GetYaxis()->SetTitleOffset(1.4);
  
  
  TH2D *proj_yz = new TH2D("proj_yz","z-y Projection;y global (cm);z global (cm)",200, (ygmin_global-mod_Ly[i_min_y])/10, (ygmax_global+mod_Ly[i_max_y])/10, 200, (zgmin_global - 100)/10, (zgmax_global + 100)/10);
  proj_yz->GetYaxis()->SetTitleOffset(1.4);




  int n_saved = 0;

  while( C->GetEntry(nevent++) && (n_saved < nevents_save) ){

    EventID = evtID;
    
    //cout << "Processing event " << nevent << endl;
    
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
    
    if( DataFormat == "HALLA" ){ //Fill CALO and SCINT information:
      double fadcsum = 0.0;
      double fshowersum = 0.0;
      double fpreshowersum=0.0;
      
      for( int icalhit=0; icalhit<nfadc; icalhit++ ){
	//int fCH = fCH[icalhit];
	int chan = fCH[icalhit] % 256;
	int slot = fCH[icalhit] >> 8;
	
	int globalchan = chan + (slot-17)*16;
	
	//cout << "chan, slot, globalchan = " << chan << ", " << slot << ", " << globalchan << endl;
	
	int maxtimesample = ftimting[icalhit];
	
	int tmin = 30, tmax=45;
	if( globalchan == 19 ) tmin = 24;
	
	if( globalchan < 21 || globalchan > 26 ){
	  if( maxtimesample >= tmin && maxtimesample <= tmax ){
	    fadcsum += fadc[icalhit];
	  }
	}
      }
      


      CALOsum = fadcsum;

      int ngoodscint = 0;
      for( int itdc=0; itdc<ntdc; itdc++ ){
	if( tCH[itdc] < 10 && 110.0< ttiming[itdc] && ttiming[itdc] < 135.0 ){
	  ngoodscint++;
	}

      }



      NGOODSCINT = ngoodscint;
      
    }
    
    map<int,int> mod_maxstripX;
    map<int,double> mod_ADCmaxX;

    map<int,int> mod_maxstripY;
    map<int,double> mod_ADCmaxY;
    
    //determine the search region on each plane, 
    //for GEP we have two contraint points, for general purposes we might not 
    //have such contraints. so just set the search center to the center of the chamber
    //and size to be something larger than the size of the chamber --WX
    search_mod.clear();
    search_x.clear();
    search_y.clear();
    
    for (unsigned int i=0; i<gridContainer.size(); i++) gridContainer[i].Clear();
    
    for (int k=0; k<nmodules; k++){
        
        search_mod.push_back(k);
        search_x.push_back(mod_x0[k]);
        search_y.push_back(mod_y0[k]);
        
        for (unsigned int i=0; i<gridContainer.size(); i++){
            if (gridContainer[i].HasModule(k)){
                //slightly larger than the actual serach region just in case some numerical issue -- WX
                gridContainer[i].SetGeoInfo(k, search_x.back(), search_y.back(), mod_Lx[k]*1.01, mod_Ly[k]*1.01);
            }
        }
    }
    
    for( int ich=0; ich<nch; ich++ ){
      int strip = Strip[ich];
      int plane = axis[ich];
      int module = moduleID[ich];

      if( DataFormat == "UVA" ) module = moduleID[ich] + 4*(planeID[ich]-1);
      
      int layer = mod_layer[module];

      modules_hit.insert( module );
      
      layers_hit.insert(layer);
      
      double ADCsamples[nADCsamples];
      ADCsamples[0] = adc0[ich];
      ADCsamples[1] = adc1[ich];
      ADCsamples[2] = adc2[ich];
      ADCsamples[3] = adc3[ich];
      ADCsamples[4] = adc4[ich];
      ADCsamples[5] = adc5[ich];

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
	ModData[module].ADCsamp_xstrips[strip].resize(nADCsamples);
	
	// ADCsum_xstrips[module][strip] = 0.0;
	// ADCsamp_xstrips[module][strip].resize(6);
	double tsamp[nADCsamples];
	double dsamp[nADCsamples];
	double dtsamp[nADCsamples];
	
	for( int isamp=0; isamp<nADCsamples; isamp++ ){

	  ADCsamples[isamp] *= mod_RYX[module]; //EXPERIMENTAL: multiply X ADC samples by ratio of Y gain to X gain: everything else proceeds as before:
	  
	  ModData[module].ADCsamp_xstrips[strip][isamp] = ADCsamples[isamp];
	  ModData[module].ADCsum_xstrips[strip] += ADCsamples[isamp];
	  

	  double tsample = 12.5+25.0*isamp;

	  tsamp[isamp] = tsample;
	  dsamp[isamp] = 20.0;
	  dtsamp[isamp] = 12.5;
	  
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

	layers_hitX.insert(layer);

	NstripX_layer[layer] += 1;

	if( mod_maxstripX.find( module ) == mod_maxstripX.end() || ModData[module].ADCsum_xstrips[strip] > mod_ADCmaxX[module] ){
	  mod_ADCmaxX[module] = ModData[module].ADCsum_xstrips[strip];
	  mod_maxstripX[module] = strip;
	}
	
      } else if( plane == mod_vplaneID[module] && keepstrip ){ //y (horizontal/short) axis
	ModData[module].ystrips.insert( strip );

	ModData[module].ADCsum_ystrips[strip] = 0.0;
	ModData[module].ADCsamp_ystrips[strip].resize(nADCsamples);

	//double tsamp[nADCsamples];
	double tsamp[nADCsamples];
	double dsamp[nADCsamples];
	double dtsamp[nADCsamples];
	
	for( int isamp=0; isamp<nADCsamples; isamp++ ){
	  ModData[module].ADCsamp_ystrips[strip][isamp] = ADCsamples[isamp];
	  ModData[module].ADCsum_ystrips[strip] += ADCsamples[isamp];


	  double tsample = 12.5+25.0*isamp;

	  tsamp[isamp] = tsample;
	  dsamp[isamp] = 20.0;
	  dtsamp[isamp] = 12.5;
	  
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
	
	
	layers_hitY.insert(layer);

	NstripY_layer[layer] += 1;

	if( mod_maxstripY.find( module ) == mod_maxstripY.end() || ModData[module].ADCsum_ystrips[strip] > mod_ADCmaxY[module] ){
	  mod_ADCmaxY[module] = ModData[module].ADCsum_ystrips[strip];
	  mod_maxstripY[module] = strip;
	}
      }

      
    }


    for( set<int>::iterator ilay=layers_hitX.begin(); ilay != layers_hitX.end(); ++ilay ){
      if( layers_hitY.find( *ilay ) != layers_hitY.end() ){
	layers_hitXY.insert( *ilay );
      }
    }
    



    
    if( layers_hitXY.size() >= TOTAL_REQUIRED_HIT ){ //enough layers hit to (possibly) form a track: only bother with clustering and attempted track finding if this is the case: 

      if( eventdisplaymode != 0 ){
	for( int imodule=0; imodule<nmodules; imodule++ ){
	  ( (TH2D*) (*hframe_modules)[imodule] )->Reset();
	}
	( (TH2D*) (proj_xz) )->Reset();
	( (TH2D*) (proj_yz) )->Reset();
      }

      //we should probably restructure the earlier part of the code to eliminate the overhead of copying all this info:
      
      map<int,clusterdata_t> mod_clusters; //mapping of found clusters by module;

      //map<int,int> N2D_match_layer;
      
      for(set<int>::iterator imod=modules_hit.begin(); imod != modules_hit.end(); ++imod ){
	//moduledata_t datatemp;
	clusterdata_t clusttemp;


	int module = *imod;
	int layer = mod_layer[module];
	
	ModData[module].modindex = module;
	ModData[module].layerindex = layer;
	int strip = mod_maxstripX[module];
	
	
	strip = mod_maxstripY[module];

	int clusterflag = find_clusters_by_module_new( ModData[module], clusttemp );

	
	
	if( clusterflag != 0 ){
	  cout << "event " << evtID << ", module " << module << " too noisy, gave up" << endl;
	}

	//hNclust_module->Fill(clusttemp.nclust2D, module );
	
	//	cout << "ending cluster finding, ncluster =  " << clusttemp.nclust2D << endl;
	
	mod_clusters[module] = clusttemp;

	NclustX_layer[mod_layer[module]]+=clusttemp.nclustx;
	NclustY_layer[mod_layer[module]]+=clusttemp.nclusty;

	Nclust2D_layer[mod_layer[module]]+=clusttemp.nclust2D;
	

      }

      int nlayers_with_2Dclust = 0;
      for( int ilay=0; ilay<nlayers; ilay++ ){
	if( Nclust2D_layer[ilay] > 0 ) nlayers_with_2Dclust++;
      }
      
      trackdata_t tracktemp;

      tracktemp.ntracks = 0;

      //      cout << "Finding tracks, event..." << evtID << endl;

      if( nlayers_with_2Dclust >= TOTAL_REQUIRED_HIT ){

	//	cout << "Finding tracks, event... " << evtID << endl;
	
	auto start = high_resolution_clock::now();
	//define the forward and backward constraint points
	TVector3 fcp(0, 0, zavg_layer[0] - 10);
	TVector3 bkp(0, 0, zavg_layer[nlayers-1] + 10);

	if( TrackingAlgorithmFlag != 0 ){
	  new_find_tracks( mod_clusters, tracktemp, fcp, bkp);
	} else {
	  find_tracks( mod_clusters, tracktemp );
	}
	auto end = high_resolution_clock::now();
	totalTime += duration_cast<nanoseconds>(end - start);
	
	


	Ntracks = tracktemp.ntracks;

	if( Ntracks > 0 && DataFormat == "HALLA" ){ //Fill CALO information:
	  double fadcsum = 0.0;
	  double fshowersum = 0.0;
	  double fpreshowersum=0.0;
	  
	  for( int icalhit=0; icalhit<nfadc; icalhit++ ){
	    //int fCH = fCH[icalhit];
	    int chan = fCH[icalhit] % 256;
	    int slot = fCH[icalhit] >> 8;
	    
	    int globalchan = chan + (slot-17)*16;

	    //cout << "chan, slot, globalchan = " << chan << ", " << slot << ", " << globalchan << endl;
	    
	    int maxtimesample = ftimting[icalhit];
	    
	  	    
	    int tmin = 30, tmax=45;
	    if( globalchan == 19 ) tmin = 24;
	    
	    if( globalchan < 21 || globalchan > 26 ){
	      if( maxtimesample >= tmin && maxtimesample <= tmax ){
		fadcsum += fadc[icalhit];
	      }
	    }
	  }

	 
	  CALOsum = fadcsum;

	  int ngoodscint = 0;
	  for( int itdc=0; itdc<ntdc; itdc++ ){
	    if( tCH[itdc] < 10 && 110.0< ttiming[itdc] && ttiming[itdc] < 135.0 ){
	      ngoodscint++;
	    }
	 
	  }
	  
	  
	  NGOODSCINT = ngoodscint;
	  
	}
	
	//for( int itrack=0; itrack<tracktemp.ntracks; itrack++ ){
	for( int itrack=0; itrack<TMath::Min(tracktemp.ntracks,1); itrack++ ){
	  

	  //Only fill residuals if we have all four layers firing (may help clarify the situation):
	  //if( tracktemp.nhitsontrack[itrack] == 4 ){
	  //&& tracktemp.Chi2NDFtrack[itrack] < 1000.0 ){
	  for( int ihit=0; ihit<tracktemp.nhitsontrack[itrack]; ihit++ ){
	    int module= tracktemp.modlist_track[itrack][ihit];
	    int layer = mod_layer[module];

	    int iclust2D = tracktemp.hitlist_track[itrack][ihit];
	    

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


	    } 
	    
	   
	  }      


	}
      }

      
      
      if( eventdisplaymode != 0 ){

	TLine Ltemp;

	Ltemp.SetLineWidth(1);
	
	double stripADCmax=1.2e4;
	int ncolors = gStyle->GetNumberOfColors();

	map<int,int> module_ipad; //pad where module info is being plotted:

	// Plots for side view projections
	c_proj->cd(1);
	( (TH2D*) (proj_xz) )->Draw();

	c_proj->cd(2);
	( (TH2D*) (proj_yz) )->Draw();
	
	for( int ilayer=0; ilayer<nlayers; ilayer++ ){
	//   c1->cd( ilayer+1 );


	  int dummy=0;
	  for( set<int>::iterator imod=modlist_layer[ilayer].begin(); imod != modlist_layer[ilayer].end(); ++imod ){
	    int imodule = *imod;

	    //int ipad = ilayer + nlayers*dummy + 1;

	    int ipad = dummy + 1 + (nlayers-ilayer-1)*maxnmodperlayer;
	    
	    module_ipad[imodule] = ipad;

	    c1->cd(ipad);

	    ( (TH2D*) (*hframe_modules)[imodule] )->Draw();
	    
	    dummy++;


	    c_proj->cd(1);
	    Ltemp.SetLineColor(kBlack);
	    
	    Ltemp.DrawLine( (mod_x0[imodule] - mod_Lx[imodule]/2)/10, (mod_z0[imodule])/10, (mod_x0[imodule] + mod_Lx[imodule]/2)/10, (mod_z0[imodule])/10 );
	    
	    c_proj->cd(2);
	    
	    Ltemp.DrawLine( (mod_y0[imodule] - mod_Ly[imodule]/2)/10, mod_z0[imodule]/10, (mod_y0[imodule] + mod_Ly[imodule]/2)/10, mod_z0[imodule]/10 );
	  
	    gPad->Modified();
	    gPad->Update();
	    c_proj->Update();
	  
	  }
	}
	
	

	for( set<int>::iterator imod=modules_hit.begin(); imod != modules_hit.end(); ++imod ){
	  
	  int layer = mod_layer[*imod];
	  int module = *imod;

	  set<int> xstrips = ModData[module].xstrips;
	  set<int> ystrips = ModData[module].ystrips;
	  
	  c1->cd(module_ipad[module]);
	  
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

	      if( fabs(uhat.X()) >= 1.e-5 ){ //strip is not along Y, compute upper and lower limits in X:
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
		TVector3 point2local(xmax_strip, ymax_strip, 0.0 );
		// transformation to global coordinates no longer applicable since we now plot "module local" coordinates:
		// TVector3 modcenter(mod_x0[module],mod_y0[module],mod_z0[module]);
		// TVector3 point1global = mod_Rot[module] * point1local + modcenter;
		// TVector3 point2global = mod_Rot[module] * point2local + modcenter;

		//How to implement a logarithmic color scale: ADCstrip over ADCmax varies from threshold/max up to 1:
		double logADCmin = log(thresh_stripsum/stripADCmax);
		double logADCmax = 0.0;

		int logADCbin = int( (log(ModData[module].ADCsum_xstrips[strip]/stripADCmax)-logADCmin)/(logADCmax-logADCmin)*double(ncolors) );

		//cout << "logADCmin, ADC, logADCbin = " << logADCmin << ", " << ADCsum_xstrips[module][strip] << ", " << logADCbin << ", color = " << gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,logADCbin))) << endl;
		int ADCbin = int( (ModData[module].ADCsum_xstrips[strip]-thresh_stripsum)/(stripADCmax-thresh_stripsum)*double(ncolors) );
		
		Ltemp.SetLineColor( gStyle->GetColorPalette( TMath::Max(0,TMath::Min(ncolors-1,ADCbin))));

		//		Ltemp.DrawLine( point1global.Y(), point1global.X(), point2global.Y(), point2global.X() );
		if( mod_coord_flag[module] ){ //X = X, Y = Y:
		  Ltemp.DrawLine( point1local.X(), point1local.Y(), point2local.X(), point2local.Y() );
		} else { //Y = X, X = Y
		  Ltemp.DrawLine( point1local.Y(), point1local.X(), point2local.Y(), point2local.X() );
		}
	      } else { //strip IS along Y: in this case, we need to make some changes: 
		// double xmax_strip = strip_center_pos.X();
		// double xmin_strip = strip_center_pos.X();
		// double ymax_strip = mod_Ly[module]/2.0;
		// double ymin_strip = -mod_Ly[module]/2.0;

		double xmax_strip = mod_Lx[module]/2.0;
		double xmin_strip = -mod_Lx[module]/2.0;
		double ymax_strip = strip_center_pos.Y();
		double ymin_strip = strip_center_pos.Y();
		
		//These are local coordinates: need to convert to global coordinates; it should be sufficient to convert the two points:
		TVector3 point1local(xmin_strip, ymin_strip, 0.0);
		TVector3 point2local(xmax_strip, ymax_strip, 0.0 );
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
		
		//		Ltemp.DrawLine( point1global.Y(), point1global.X(), point2global.Y(), point2global.X() );
		if( mod_coord_flag[module] ){ //X = X, Y = Y:
		  Ltemp.DrawLine( point1local.X(), point1local.Y(), point2local.X(), point2local.Y() );
		} else { //Y = X, X = Y
		  Ltemp.DrawLine( point1local.Y(), point1local.X(), point2local.Y(), point2local.X() );
		}

		//Ltemp.DrawLine( point1global.Y(), point1global.X(), point2global.Y(), point2global.X() );
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
	    TVector3 vhat( mod_Pxv[module], mod_Pyv[module], 0.0); //(0, 1, 0)
	    TVector3 zaxis(0,0,1);

	    //cout << "vhat = " << endl;
	    //vhat.Print();
	    
	    //Unit vector perp. to U; i.e., ALONG strip:
	    TVector3 vperphat = -zaxis.Cross(vhat).Unit();

	    // cout << "vperphat = " << endl;
	    // vperphat.Print();
	    
	    TVector3 strip_center_pos = vlocal*vhat;

	    // cout << "strip center pos = " << endl;
	    // strip_center_pos.Print();
	    
	    if( fabs( strip_center_pos.X() ) <= mod_Lx[module]/2.0 &&
		fabs( strip_center_pos.Y() ) <= mod_Ly[module]/2.0 ){
	      //Then we can draw this strip:

	      //I think this covers all possible scenarios, SHOULD avoid divide-by-zero errors:
	      if( fabs(vhat.X()) >= 1.e-5 ){ //strip is not along Y, compute upper and lower limits in X:
		double xmin_strip = strip_center_pos.X() + vperphat.X()/vperphat.Y() * (-mod_Ly[module]/2.0-strip_center_pos.Y());
		double xmax_strip = strip_center_pos.X() + vperphat.X()/vperphat.Y() * (mod_Ly[module]/2.-strip_center_pos.Y());

		double ymin_strip = -mod_Ly[module]/2.0;
		double ymax_strip = mod_Ly[module]/2.0;
		
		xmin_strip = (xmin_strip < -mod_Lx[module]/2.0) ? -mod_Lx[module]/2.0 : xmin_strip;
		xmax_strip = (xmax_strip > mod_Lx[module]/2.0) ? mod_Lx[module]/2.0 : xmax_strip;

		//These conditions shouldn't be able to be satisfied if vperphat.X == 0; thus avoiding divide-by-zero errors:
		if( xmin_strip == -mod_Lx[module]/2.0 ) ymin_strip = strip_center_pos.Y() + vperphat.Y()/vperphat.X() * ( xmin_strip - strip_center_pos.X() );
		if( xmax_strip == mod_Lx[module]/2.0 ) ymax_strip = strip_center_pos.Y() + vperphat.Y()/vperphat.X() * ( xmax_strip - strip_center_pos.X() );

		//cout << "ystrip (xmin,ymin,xmax,ymax)=(" << xmin_strip << ", " << ymin_strip << ", " << xmax_strip << ", " << ymax_strip << ")" << endl;
		
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

		//Ltemp.DrawLine( point1global.Y(), point1global.X(), point2global.Y(), point2global.X() );
		//		Ltemp.DrawLine( point1global.Y(), point1global.X(), point2global.Y(), point2global.X() );
		if( mod_coord_flag[module] ){ //X = X, Y = Y:
		  Ltemp.DrawLine( point1local.X(), point1local.Y(), point2local.X(), point2local.Y() );
		} else { //Y = X, X = Y
		  Ltemp.DrawLine( point1local.Y(), point1local.X(), point2local.Y(), point2local.X() );
		}
		
	      } else { //strip IS along Y:

		double xmax_strip = mod_Lx[module]/2.0;
		double xmin_strip = -mod_Lx[module]/2.0;
		double ymax_strip = strip_center_pos.Y();
		double ymin_strip = strip_center_pos.Y();


		//	cout << "strip along Y, (xmin,ymin,xmax,ymax) = (" << xmin_strip << ", " << ymin_strip << ", " << xmax_strip << ", " << ymax_strip << ")" << endl;
		
		//These are local coordinates: need to convert to global coordinates; it should be sufficient to convert the two points:
		TVector3 point1local(xmin_strip, ymin_strip, 0.0);
		TVector3 point2local(xmax_strip, ymax_strip, 0.0 );
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

		//Ltemp.DrawLine( point1global.Y(), point1global.X(), point2global.Y(), point2global.X() );
		//		Ltemp.DrawLine( point1global.Y(), point1global.X(), point2global.Y(), point2global.X() );
		if( mod_coord_flag[module] ){ //X = X, Y = Y:
		  Ltemp.DrawLine( point1local.X(), point1local.Y(), point2local.X(), point2local.Y() );
		} else { //Y = X, X = Y
		  Ltemp.DrawLine( point1local.Y(), point1local.X(), point2local.Y(), point2local.X() );
		}
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
	  
	  if(clusttemp.nclust2D > 20) continue; //Some bad events have too many clusters and slow down the script

	  for( int iclust=0; iclust<clusttemp.nclust2D; iclust++ ){
	    

	    //global hit position:	    
	    double xhit = clusttemp.xglobal2D[iclust];
	    double yhit = clusttemp.yglobal2D[iclust];
	    double zhit = clusttemp.zglobal2D[iclust];
	    
	    //local hit position:	    
	    double uhit = clusttemp.xclust2D[iclust];
	    double vhit = clusttemp.yclust2D[iclust];
	    
	    
	    double det = mod_Pxu[module]*mod_Pyv[module] - mod_Pyu[module]*mod_Pxv[module];
	    
	    double Xlocal =  (mod_Pyv[module]*uhit - mod_Pyu[module]*vhit)/det;
	    double Ylocal =  (mod_Pxu[module]*vhit - mod_Pxv[module]*uhit)/det;

	    
	    c_proj->cd(1);
	    
	    Mhit.DrawMarker( xhit/10, zhit/10 );

	    c_proj->cd(2);
	    
	    Mhit.DrawMarker( yhit/10, zhit/10 );
	    
	    gPad->Modified();
	    gPad->Update();
	    c_proj->Update();


	    c1->cd(module_ipad[module]);

	    if( mod_coord_flag[module] ){
	      Mhit.DrawMarker( Xlocal, Ylocal );
	    } else {
	      Mhit.DrawMarker( Ylocal, Xlocal );
	    }
	    
	    gPad->Modified();
	    gPad->Update();
	    //gSystem->ProcessEvents();
	    c1->Update();
	    
	  }
	}
      
	if( tracktemp.ntracks > 0 && tracktemp.ntracks<9){ //Draw track information for each layer:
	  for( int itrack=0; itrack<tracktemp.ntracks; itrack++ ){

	    c_proj->cd(1);
	    Ltemp.DrawLine((tracktemp.Xtrack[itrack]+tracktemp.Xptrack[itrack]*(-5))/10, -5, (tracktemp.Xtrack[itrack]+tracktemp.Xptrack[itrack]*(zgmax_global + 5))/10, zgmax_global/10 + 5);
	    
	    c_proj->cd(2);
	    Ltemp.DrawLine((tracktemp.Ytrack[itrack]+tracktemp.Yptrack[itrack]*(-5))/10, -5, (tracktemp.Ytrack[itrack]+tracktemp.Yptrack[itrack]*(zgmax_global + 5))/10, zgmax_global/10 + 5);

	    
	    gPad->Modified();
	    gPad->Update();
	    c_proj->Update();
	    
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

	      TVector3 modcenter_global( mod_x0[module], mod_y0[module], mod_z0[module] );

	      //compute local X and Y of the track within the module:
	      TVector3 trackpos_local = mod_Rotinv[module]*(trackpos_global-modcenter_global); 

	      xtrack = trackpos_local.X();
	      ytrack = trackpos_local.Y();
	      
	      //now how about computing 
	      
	      c1->cd(module_ipad[module]);

	      // Mhit.DrawMarker( yhit_local_stripcoord, xhit_local_stripcoord+nstripsx*(module%3) );

	      Mtrack.SetMarkerStyle( mstyle_track[itrack] );
	      Mtrack.SetMarkerColor( mcolor_track[itrack] );

	      if( mod_coord_flag[module] ){
		Mtrack.DrawMarker( xtrack, ytrack );
	      } else {
		Mtrack.DrawMarker( ytrack, xtrack );
	      }
	      //Mhit.Draw("SAME");
	      //Mtrack.Draw("SAME");
		  
	      gPad->Modified();
	      gPad->Update();
	      //gSystem->ProcessEvents();
	      c1->Update();


	      c_proj->cd(1);
	      Mtrack.DrawMarker( trackpos_global.X()/10, trackpos_global.Z()/10);
	      
	      c_proj->cd(2);
	      Mtrack.DrawMarker( trackpos_global.Y()/10, trackpos_global.Z()/10);
	      
	      gPad->Modified();
	      gPad->Update();
	      c_proj->Update();
	      
	    }
	  }
	  
	}
	gSystem->ProcessEvents();

	n_saved++;

	if(n_saved == 1) {
	  c1->Print("./output/" + outfilename + ".pdf(");
	  c_proj->Print("./output/" + outfilename + ".pdf");
	}
	else if(n_saved == nevents_save) {
	  c1->Print("./output/" + outfilename + ".pdf");
	  c_proj->Print("./output/" + outfilename + ".pdf)");
	}
	else {
	  c1->Print("./output/" + outfilename + ".pdf");
	  c_proj->Print("./output/" + outfilename + ".pdf");
	  
	}	


      }
	  
	  
      //gSystem->Sleep(2500);

      //	  c1->Update();
	  
	  
    } else {
      
    }

    
  
    
  }

   
  
  auto program_end = high_resolution_clock::now();
  
  cout<<"total time spent in tracking: "<<(double)totalTime.count() / 1e9<<" seconds "<<endl;
  
  auto program_time = duration_cast<nanoseconds>(program_end - program_start); 
  
  cout<<"total program execution time: "<<(double)program_time.count() / 1e9<<" seconds" <<endl;
  
  cout<<"fraction of time spent in track-finding: "<<((double)totalTime.count()/(double)program_time.count())<<endl;
}

