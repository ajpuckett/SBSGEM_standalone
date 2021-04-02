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
#include <algorithm>
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

long NSKIPPED;

int TrackingAlgorithmFlag = 0; //0 = OLD, "brute force" algorithm, 1 = NEW, "grid container" algorithm"

int nstripsx = 1280;
int nstripsy = 1024;
int nmodules = 12;
int nlayers  = 4;
int nADCsamples = 6;

double strip_pitch = 0.4; //mm
double sigma_hitshape = strip_pitch; //mm

const double PI = TMath::Pi();

//These default values are guesses:
//int maxnhitspercluster=10;
int maxnstripspercluster=7;
int maxnstripXpercluster=9;
int maxnstripYpercluster=7;

int maxneighborsX=4; //+/-4 around local maximum:
int maxneighborsY=4; //+/-3 around local maximum:

//Defaults are to keep all strips even if the max occurs in the "wrong" time sample:
int sampmin_accept = -1; //If max ADC occurs below this bin, reject strip
int sampmax_accept = 6;  //If max ADC occurs above this bin, reject strip

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
double thresh_2ndmax_fraction = 0.2; //threshold to flag 2nd maximum as a fraction of first maximum
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
  map<int,double> ADCsum_xstrips_goodsamp;
  map<int,double> ADCsum_ystrips_goodsamp;
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

  void Clear(){ //maybe managing these things a bit better will help with the random segfaults?
    modindex=-1;
    layerindex=-1;
    xstrips.clear();
    ystrips.clear();
    ADCsum_xstrips.clear();
    ADCsum_ystrips.clear();
    ADCsum_xstrips_goodsamp.clear();
    ADCsum_ystrips_goodsamp.clear();
    ADCsamp_xstrips.clear();
    ADCsamp_ystrips.clear();
    Bestmatch_xstrips.clear();
    Bestmatch_ystrips.clear();
    BestCor_xstrips.clear();
    BestCor_ystrips.clear();
    xstrips_filtered.clear();
    ystrips_filtered.clear();
    ADCmax_xstrips.clear();
    ADCmax_ystrips.clear();

    isampmax_xstrips.clear();
    isampmax_ystrips.clear();

    Tfit_xstrips.clear();
    Tfit_ystrips.clear();
    dTfit_xstrips.clear();
    dTfit_ystrips.clear();
    Afit_xstrips.clear();
    dAfit_xstrips.clear();
    taufit_xstrips.clear();
    dtaufit_xstrips.clear();

    Tfit_Chi2NDF_xstrips.clear();

    Afit_ystrips.clear();
    dAfit_ystrips.clear();
    taufit_ystrips.clear();
    dtaufit_ystrips.clear();

    Tfit_Chi2NDF_ystrips.clear();

    Tmean_xstrips.clear();
    Tsigma_xstrips.clear();

    Tmean_ystrips.clear();
    Tsigma_ystrips.clear();

    Tmean_xstrips_walkcor.clear();
    Tsigma_xstrips_walkcor.clear();

    Tmean_ystrips_walkcor.clear();
    Tsigma_ystrips_walkcor.clear();
  }
  
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
  //Not clear we need these:
  //  vector<vector<double> > ADCsamp_xclust2D;
  //  vector<vector<double> > ADCsamp_yclust2D;
  
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
  vector<vector<double> > xstripADCsum; //individual strip ADCs, including effect of splitting fraction, as applicable
  
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
  vector<vector<double> > ystripADCsum; //individual strip ADCs, including effect of splitting fraction, as applicable

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
  void Clear(){

    nkeep2D = 0;
    nclust2D = 0;
    
    itrack_clust2D.clear();
    ixclust2D.clear();
    iyclust2D.clear();
    nstripx2D.clear();
    nstripy2D.clear();
    xclust2D.clear();
    yclust2D.clear();
    xclust2Dcorr.clear();
    yclust2Dcorr.clear();
    Eclust2D.clear();
    tclust2D.clear();
    tclust2Dwalkcorr.clear();
    dEclust2D.clear();
    dtclust2D.clear();
    dtclust2Dwalkcorr.clear();
    CorrCoeff2D.clear();
    CorrCoeffMaxStrips.clear();
    xglobal2D.clear();
    yglobal2D.clear();
    zglobal2D.clear();
    xmom2D.clear();
    ymom2D.clear();
    keepclust2D.clear();

    ADCsamp_xclust.clear();
    ADCsamp_yclust.clear();

    nclustx=0;
    nstripx.clear();
    ixstriplo.clear();
    ixstriphi.clear();
    ixstripmax.clear();
    xmean.clear();
    xsigma.clear();
    totalchargex.clear();
    txmean.clear();
    txsigma.clear();
    xstripADCsum.clear();

    nclusty=0;
    nstripy.clear();
    iystriplo.clear();
    iystriphi.clear();
    iystripmax.clear();
    ymean.clear();
    ysigma.clear();
    totalchargey.clear();
    tymean.clear();
    tysigma.clear();
    ystripADCsum.clear();
    
  }
  
  
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

  void Clear(){
    ntracks = 0;
    nhitsontrack.clear();
    modlist_track.clear();
    hitlist_track.clear();
    residx_hits.clear();
    residy_hits.clear();
    eresidx_hits.clear();
    eresidy_hits.clear();
    Xtrack.clear();
    Ytrack.clear();
    Xptrack.clear();
    Yptrack.clear();
    Chi2NDFtrack.clear();

  }
  
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

map<int,int> mod_nAPVX;
map<int,int> mod_nAPVY;
map<int,vector<double> > mod_Xgain;
map<int,vector<double> > mod_Ygain;

//Also useful to define number of modules per layer and list of modules by layer:
map<int,int> nmodules_layer;
map<int,set<int> > modlist_layer;

// Tracking layer combinatorics: populate these arrays once so we don't do it every event:
//For each possible number of layers from 3 up to the total number of layers, we list all possible combinations of n layers
map<int,vector<vector<int> > > layercombos;


vector<double> zavg_layer;
map<int,double> xgmin_layer, xgmax_layer, ygmin_layer, ygmax_layer;
//map<int,double> xbinwidth_layer,ybinwidth_layer;
//map<int,int> nbinsx_layer, nbinsy_layer;
double gridxbinwidth=10.0,gridybinwidth=10.0; //10 mm, eventually user can override
map<int,int> gridnbinsx_layer, gridnbinsy_layer;
map<int,double> gridxmin_layer, gridxmax_layer, gridymin_layer, gridymax_layer;

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
  
int find_clusters_by_module_newnew( moduledata_t mod_data, clusterdata_t &clust_data ){
  //In this approach, let us do 1D cluster-finding in X and Y, and then combine 1D clusters into 2D hit candidates:
  //This might affect some of the other routines that currently expect identical numbers of X and Y clusters:

  //

  clust_data.Clear();
  
  int module = mod_data.modindex;
  int layer = mod_data.layerindex;
  
  int nlocalmaxX=0, nlocalmaxY=0;

  clust_data.nclustx = 0;
  clust_data.nclusty = 0;
  clust_data.nclust2D = 0;
  clust_data.modindex = module;
  clust_data.layerindex = layer;
  clust_data.nkeep2D = 0;

  //create  temporary structures to hold local maxima:
  map<int,bool> islocalmax; //key = strip, value = true or false;
  set<int> localmaxima;     //list of strips that are local maxima
  //map<int,double> ADCmaxima; //key = strip, value = ADC sum of local max strip
  //vector<int> sortedmaxima;  //array of local max strips sorted by ADC amplitude: Do we really need this? perhaps not:
  
  //Should we do one loop to flag all local maxima, then a second loop to check for multiple
  //local maxima in a contiguous grouping of strips within maxneighbors of each other?

  islocalmax.clear();
  localmaxima.clear();
  
  
  for( set<int>::iterator ix=mod_data.xstrips.begin(); ix != mod_data.xstrips.end(); ++ix ){
    int strip = *ix;

    //    islocalmax[strip] = false;
    islocalmax.insert(std::pair<int,bool>(strip,false) );
    
    double sumstrip = mod_data.ADCsum_xstrips_goodsamp[strip];
    //check if strip is local max based on sums in immediately adjacent strips:
    double sumleft = 0.0;
    double sumright = 0.0;
    if( mod_data.xstrips.find( strip-1 ) != mod_data.xstrips.end() ){
      sumleft = mod_data.ADCsum_xstrips_goodsamp[strip-1];
    }
    if( mod_data.xstrips.find( strip+1 ) != mod_data.xstrips.end() ){
      sumright = mod_data.ADCsum_xstrips_goodsamp[strip+1]; 
    }

    if( sumstrip >= sumleft && sumstrip >= sumright ){ //new local maximum: add it to the list:
      islocalmax[strip] = true;
      localmaxima.insert(strip);
      //ADCmaxima[strip] = sumstrip;

      
    }
  }

  //Better idea: calculate the "prominence" of each peak other than the highest one, and apply a threshold
  //on the prominence to keep the peak as a seed for an independent cluster:

  //let's avoid erasing things from maps and sets until we get through the entire sorted list
  vector<int> peakstoerase;

  // for( int imax=0; imax<sortedmaxima.size(); imax++ ){ //these maxima are sorted in descending order of ADC
  //  int stripmax = sortedmaxima[imax]; 

  for( set<int>::iterator ix=localmaxima.begin(); ix != localmaxima.end(); ++ix ){
    int stripmax = *ix;

    double ADCmax = mod_data.ADCsum_xstrips_goodsamp[stripmax];
    
    double prominence = ADCmax;
    
    int striplo=stripmax,striphi=stripmax;

    double ADCminright = ADCmax;
    bool higherpeakright=false;
    int peakright=-1;
    
    while( mod_data.xstrips.find(striphi+1) != mod_data.xstrips.end() ){
      striphi++;
      if( mod_data.ADCsum_xstrips_goodsamp[striphi] < ADCminright && !higherpeakright ){ //as long as we haven't yet found a higher peak to the right, this is the lowest point between the current maximum and the next higher peak to the right:
	ADCminright = mod_data.ADCsum_xstrips_goodsamp[striphi];
      }

      if( islocalmax[striphi] && mod_data.ADCsum_xstrips_goodsamp[striphi] > ADCmax ){ //then this peak is in a contiguous group with another higher peak to the right:
	higherpeakright=true;
	peakright = striphi;
      }
    }

    double ADCminleft = ADCmax;
    bool higherpeakleft = false;
    int peakleft=-1;
    
    while( mod_data.xstrips.find(striplo-1) != mod_data.xstrips.end() ){
      striplo--;
      if( mod_data.ADCsum_xstrips_goodsamp[striplo] < ADCminleft && !higherpeakleft ){ //as long as we haven't yet found a higher peak to the right, this is the lowest point between the current maximum and the next higher peak to the right:
	ADCminleft = mod_data.ADCsum_xstrips_goodsamp[striplo];
      }

      if( islocalmax[striplo] && mod_data.ADCsum_xstrips_goodsamp[striplo] > ADCmax ){ //then this peak is in a contiguous groupd of strips with another higher peak to the left
	higherpeakleft=true;
	peakleft = striplo;
      }
    }

    // if this peak is contiguous with another higher peak to the left or the right, calculate prominence,
    // and reject if the prominence of this peak is below the theshold for significance, or below
    // a certain minimum percentage of the peak height:
    
    int nsamplesinsum = std::min(nADCsamples,sampmax_accept-sampmin_accept+1);

    double sigma_sum = sqrt(double(nsamplesinsum))*sigma_sample; // about 50 ADC channels

    //cout << "sigma_sum = " << sigma_sum << endl;
    
    if( !higherpeakleft ) ADCminleft = 0.0;
    if( !higherpeakright ) ADCminright = 0.0;
    
    if( higherpeakright || higherpeakleft ){ //contiguous with higher peaks on either the left or the right or both:
      prominence = ADCmax - std::max( ADCminleft, ADCminright );

      //cout << "ADC max, peak prominence = " << ADCmaxima[stripmax] << ", " << prominence << endl;
      
      bool peak_close=false;
      if( higherpeakleft && abs( peakleft - stripmax ) <= 2*maxneighborsX ) peak_close = true;
      if( higherpeakright && abs( peakright - stripmax ) <= 2*maxneighborsX ) peak_close = true;

      //peak_close = true;
      
      if( (prominence < thresh_2ndmax_nsigma * sigma_sum ||
	   prominence/ADCmax < thresh_2ndmax_fraction) && peak_close ){
	//	localmaxima.erase( stripmax );
	//islocalmax[stripmax] = false;
	//ADCmaxima.erase( stripmax );
	peakstoerase.push_back( stripmax );
      }
    }
  }
  //I think this is safer...

  for( int ipeak=0; ipeak<peakstoerase.size(); ipeak++ ){
    localmaxima.erase( peakstoerase[ipeak] );
    if( islocalmax.find(peakstoerase[ipeak]) != islocalmax.end() ){
      islocalmax[peakstoerase[ipeak]] = false;
    }
  }

  //cout << "After erasing insignificant peaks, nmaxima = " << localmaxima.size() << endl;
  
  
  //how do we deal with overlap? One idea is to collect all strips within +/- maxneighbors of each local max,
  //and then come up with a formula for assigning a fractional contribution of each nearby local maximum to each strip
  //How would we do this?
  //loop over all the strips
  //Count and collect a list of all local maxima within +/- maxneighbors of each strip
  //If we treat each local maximum as a separate hit, then
  // For each local maximum within +/- maxneighbors of any given strip,
  // which happens to be in a contiguous grouping with that strip, estimate a fractional contribution to that strip signal:
  //the contribution should be proportional to the ADC value of the nearby maximum, and should fall off with distance
  //according to a Lorentzian function with a width of approximately one strip pitch:
  // f(x-xstrip) = 1/(1 + (x-xstrip)^2/sigma^2)
  // For each contributing local maximum, we would express the ADC sum as:
  // fraction (hit i) = Ai/(1+(xi-xstrip)^2/sigma^2) / sum_j Aj /(1+(xj-xstrip)^2/sigma^2)
  // And then the ADC contribution of strip k to hit i would be f_ik * ADCk

  //The method we are using now finds too many clusters:
  //We should iterate over the sorted list of local maxima, and for each one after the first one, check
  //for overlap: if a local maximum is within +/- maxneighbors of another local maximum with larger amplitude,
  //then we have to check its "significance"
  
  for( set<int>::iterator ix = localmaxima.begin(); ix != localmaxima.end(); ++ix ){
    int stripmax = *ix;
    int striplo = stripmax;
    int striphi = stripmax;

    double ADCmax = mod_data.ADCsum_xstrips_goodsamp[stripmax];
    
    while( mod_data.xstrips.find( striplo-1 ) != mod_data.xstrips.end() &&
	   stripmax-striplo < maxneighborsX ) striplo--;
    while( mod_data.xstrips.find( striphi+1 ) != mod_data.xstrips.end() &&
	   striphi - stripmax < maxneighborsX ) striphi++;

    // ixlo_maxima[stripmax] = striplo;
    // ixhi_maxima[stripmax] = striphi;

    int nstrips_maxima = striphi-striplo+1; 
    
    //for each maximum, loop over all strips 

    

    double sumx = 0.0, sumx2 = 0.0, sumADC = 0.0, sumt = 0.0, sumt2 = 0.0;

    map<int,double> splitfraction;

    splitfraction.clear();
    
    vector<double> stripADCsum(nstrips_maxima);
    
    for( int istrip=striplo; istrip<=striphi; istrip++ ){
      int nmax_strip = 1;
      double sumweight = ADCmax/(1.0 + pow( (stripmax-istrip)*mod_ustrip_pitch[module]/sigma_hitshape, 2 ) );
      //double sweight = ADCmaxima[stripmax]/(1.0 + pow( (stripmax-istrip)*mod_ustrip_pitch[module]/sigma_hitshape, 2 ) );
      double maxweight = sumweight;

      //for every strip, search for any other (contiguous) local maxima within +/- maxneighbors and calculate the "weight" of any other local maxima in this strip
      for( int jstrip=(istrip-maxneighborsX); jstrip<=(istrip+maxneighborsX); jstrip++ ){
	if( localmaxima.find( jstrip ) != localmaxima.end() && jstrip != stripmax &&
	    mod_data.xstrips.find(jstrip) != mod_data.xstrips.end() ){ //then this strip has "overlap" with another local max: check if every strip in between fired:	  
	  
	  sumweight += mod_data.ADCsum_xstrips_goodsamp[jstrip]/(1.0 + pow( (jstrip-istrip)*mod_ustrip_pitch[module]/sigma_hitshape, 2 ) );
	}
      }

      //calculate "split fraction" for this strip as the ratio of the weight of the maximum CURRENTLY under consideration in this strip
      //to the sum of all weights of any other local maxima within +/- maxneighbors of this strip: 
      double ADCfrac_strip = maxweight/sumweight;

      splitfraction[istrip] = ADCfrac_strip;
      
      //Now we can start populating the 1D cluster info:
      double ulocal = ( istrip + 0.5 - 0.5*mod_nstripsu[module] ) * mod_ustrip_pitch[module];
      double ADCstrip = mod_data.ADCsum_xstrips_goodsamp[istrip]*ADCfrac_strip;
      double tstrip = mod_data.Tmean_xstrips[istrip];

      stripADCsum[istrip-striplo] = ADCstrip;
      
      sumx += ulocal * ADCstrip;
      sumADC += ADCstrip;
      sumx2 += pow(ulocal,2)*ADCstrip;
      sumt += tstrip * ADCstrip;
      sumt2 += pow(tstrip,2)*ADCstrip;
   
    }

    // if( sumADC >= thresh_clustersum ){ //don't use if the cluster ADC sum is below threshold:
    if( true ){

      clust_data.xstripADCsum.push_back( stripADCsum );
      
      clust_data.nstripx.push_back( nstrips_maxima );
      clust_data.ixstriplo.push_back( striplo );
      clust_data.ixstriphi.push_back( striphi );
      clust_data.ixstripmax.push_back( stripmax );
      
      clust_data.xmean.push_back( sumx / sumADC );
      clust_data.xsigma.push_back( sqrt( sumx2/sumADC - pow(sumx/sumADC,2) ) );
      clust_data.totalchargex.push_back( sumADC );
      clust_data.txmean.push_back( sumt / sumADC );
      clust_data.txsigma.push_back( sqrt( sumt2/sumADC - pow( sumt/sumADC, 2 ) ) );

      vector<double> xADCsamples(nADCsamples);
      for( int isamp=0; isamp<nADCsamples; isamp++ ){
	xADCsamples[isamp] = 0.0;
	for( int istrip=striplo; istrip<=striphi; istrip++ ){
	  xADCsamples[isamp] += mod_data.ADCsamp_xstrips[istrip][isamp]*splitfraction[istrip];
	}
      }

      clust_data.ADCsamp_xclust.push_back( xADCsamples );

      
      
      clust_data.nclustx++;
    }
  }

  //  cout << "finished X strips" << endl;

  //Now we need to repeat for Y: 

  //clear out containers for local maxima:
  islocalmax.clear();
  localmaxima.clear(); 
  // ADCmaxima.clear();
  // sortedmaxima.clear();

  //  cout << "starting Y strips" << endl;
  
  for( set<int>::iterator iy=mod_data.ystrips.begin(); iy != mod_data.ystrips.end(); ++iy ){
    int strip = *iy;

    //islocalmax[strip] = false;
    islocalmax.insert( std::pair<int,bool>(strip,false) );
    
    double sumstrip = mod_data.ADCsum_ystrips_goodsamp[strip];
    //check if strip is local max based on sums in immediately adjacent strips:
    double sumleft = 0.0;
    double sumright = 0.0;
    if( mod_data.ystrips.find( strip-1 ) != mod_data.ystrips.end() ){
      sumleft = mod_data.ADCsum_ystrips_goodsamp[strip-1];
    }
    if( mod_data.ystrips.find( strip+1 ) != mod_data.ystrips.end() ){
      sumright = mod_data.ADCsum_ystrips_goodsamp[strip+1]; 
    }

    if( sumstrip >= sumleft && sumstrip >= sumright ){ //new local maximum: add it to the list:
      islocalmax[strip] = true;
      localmaxima.insert(strip);
      
    }
  }

  //Better idea: calculate the "prominence" of each peak other than the highest one, and apply a threshold
  //on the prominence to keep the peak as a seed for an independent cluster:

  peakstoerase.clear();
  
  //for( int imax=0; imax<sortedmaxima.size(); imax++ ){ //these maxima 
  for( set<int>::iterator iy = localmaxima.begin(); iy != localmaxima.end(); ++iy ){
    int stripmax = *iy;

    double ADCmax = mod_data.ADCsum_ystrips_goodsamp[stripmax];
    double prominence = ADCmax;
    
    int striplo=stripmax,striphi=stripmax;

    double ADCminright = ADCmax;
    bool higherpeakright=false;
    int peakright = -1; 
    while( mod_data.ystrips.find(striphi+1) != mod_data.ystrips.end() ){
      striphi++;
      if( mod_data.ADCsum_ystrips_goodsamp[striphi] < ADCminright && !higherpeakright ){
	ADCminright = mod_data.ADCsum_ystrips_goodsamp[striphi];
      }

      if( islocalmax[striphi] && mod_data.ADCsum_ystrips_goodsamp[striphi] > ADCmax ){ //then this peak is in a contiguous group with another higher peak to the right:
	higherpeakright=true;
	peakright = striphi;
      }
    }

    double ADCminleft = ADCmax;
    bool higherpeakleft = false;
    int peakleft = -1;
    while( mod_data.ystrips.find(striplo-1) != mod_data.ystrips.end() ){
      striplo--;
      if( mod_data.ADCsum_ystrips_goodsamp[striplo] < ADCminleft && !higherpeakleft ){
	ADCminleft = mod_data.ADCsum_ystrips_goodsamp[striplo];
      }

      if( islocalmax[striplo] && mod_data.ADCsum_ystrips_goodsamp[striplo] > ADCmax ){ //then this peak is in a contiguous groupd of strips with another higher peak to the left
	higherpeakleft=true;
	peakleft = striplo;
      }
    }

    // if this peak is contiguous with another higher peak to the left or the right, calculate prominence,
    // and reject if the prominence of this peak is below the theshold for significance, or below
    // a certain minimum percentage of the peak height:
    
    int nsamplesinsum = std::min(nADCsamples,sampmax_accept-sampmin_accept+1);

    double sigma_sum = sqrt(double(nsamplesinsum))*sigma_sample; // about 50 ADC channels

    if( !higherpeakleft ) ADCminleft = 0.0;
    if( !higherpeakright ) ADCminright = 0.0;
    
    if( higherpeakright || higherpeakleft ){ //contiguous with higher peaks on either the left or the right or both:
      prominence = ADCmax - std::max( ADCminleft, ADCminright );

      bool peak_close=false;
      if( higherpeakleft && abs( peakleft - stripmax ) <= 2*maxneighborsY ) peak_close = true;
      if( higherpeakright && abs( peakright - stripmax ) <= 2*maxneighborsY ) peak_close = true;
      
      if( (prominence < thresh_2ndmax_nsigma * sigma_sum ||
	   prominence/ADCmax < thresh_2ndmax_fraction) && peak_close ){
	// localmaxima.erase( stripmax );
	// islocalmax[stripmax] = false;
	// ADCmaxima.erase( stripmax );

	peakstoerase.push_back( stripmax );
      }
    } 
  }

   for( int ipeak=0; ipeak<peakstoerase.size(); ipeak++ ){
    localmaxima.erase( peakstoerase[ipeak] );

    if( islocalmax.find(peakstoerase[ipeak]) != islocalmax.end() ){
      islocalmax[peakstoerase[ipeak]] = false;
    }
  }
  
  
  for( set<int>::iterator iy = localmaxima.begin(); iy != localmaxima.end(); ++iy ){
    int stripmax = *iy;
    int striplo = stripmax;
    int striphi = stripmax;

    double ADCmax = mod_data.ADCsum_ystrips_goodsamp[stripmax];
    
    while( mod_data.ystrips.find( striplo-1 ) != mod_data.ystrips.end() &&
	   stripmax-striplo < maxneighborsY ) striplo--;
    while( mod_data.ystrips.find( striphi+1 ) != mod_data.ystrips.end() &&
	   striphi - stripmax < maxneighborsY ) striphi++;

    // iylo_maxima[stripmax] = striplo;
    // iyhi_maxima[stripmax] = striphi;

    int nstrips_maxima = striphi-striplo+1; 
    
    //for each maximum, loop over all strips 

   

    double sumy = 0.0, sumy2 = 0.0, sumADC = 0.0, sumt = 0.0, sumt2 = 0.0;

    map<int,double> splitfraction;

    splitfraction.clear();
    
    vector<double> stripADCsum(nstrips_maxima);
    
    for( int istrip=striplo; istrip<=striphi; istrip++ ){
      int nmax_strip = 1;
      double sumweight = ADCmax/(1.0 + pow( (stripmax-istrip)*mod_vstrip_pitch[module]/sigma_hitshape, 2 ) );
      //double sweight = ADCmaxima[stripmax]/(1.0 + pow( (stripmax-istrip)*mod_ustrip_pitch[module]/sigma_hitshape, 2 ) );
      double maxweight = sumweight;

      //for every strip, search for any other local maxima within +/- maxneighbors and calculate the "weight" of any other local maxima in this strip
      for( int jstrip=istrip-maxneighborsY; jstrip<=istrip+maxneighborsY; jstrip++ ){
	if( localmaxima.find( jstrip ) != localmaxima.end() && jstrip != stripmax &&
	    mod_data.ystrips.find(jstrip) != mod_data.ystrips.end() ){ //then this strip has "overlap" with another local max
	  sumweight += mod_data.ADCsum_ystrips_goodsamp[jstrip]/(1.0 + pow( (jstrip-istrip)*mod_vstrip_pitch[module]/sigma_hitshape, 2 ) );
	}
      }

      //calculate "split fraction" for this strip as the ratio of the weight of the maximum CURRENTLY under consideration in this strip
      //to the sum of all weights of any other local maxima within +/- maxneighbors of this strip: 
      double ADCfrac_strip = maxweight/sumweight;

      splitfraction[istrip] = ADCfrac_strip;
      
      //Now we can start populating the 1D cluster info:
      double vlocal = ( istrip + 0.5 - 0.5*mod_nstripsv[module] ) * mod_vstrip_pitch[module];
      double ADCstrip = mod_data.ADCsum_ystrips_goodsamp[istrip]*ADCfrac_strip;
      double tstrip = mod_data.Tmean_ystrips[istrip];

      stripADCsum[istrip-striplo] = ADCstrip;
      
      sumy += vlocal * ADCstrip;
      sumADC += ADCstrip;
      sumy2 += pow(vlocal,2)*ADCstrip;
      sumt += tstrip * ADCstrip;
      sumt2 += pow(tstrip,2)*ADCstrip;
   
    }

    //if( sumADC >= thresh_clustersum ){
    if( true ){
      clust_data.ystripADCsum.push_back( stripADCsum );
      
      clust_data.nstripy.push_back( nstrips_maxima );
      clust_data.iystriplo.push_back( striplo );
      clust_data.iystriphi.push_back( striphi );
      clust_data.iystripmax.push_back( stripmax );
      
      clust_data.ymean.push_back( sumy / sumADC );
      clust_data.ysigma.push_back( sqrt( sumy2/sumADC - pow(sumy/sumADC,2) ) );
      clust_data.totalchargey.push_back( sumADC );
      clust_data.tymean.push_back( sumt / sumADC );
      clust_data.tysigma.push_back( sqrt( sumt2/sumADC - pow( sumt/sumADC, 2 ) ) );

      vector<double> yADCsamples(nADCsamples);
      for( int isamp=0; isamp<nADCsamples; isamp++ ){
	yADCsamples[isamp] = 0.0;
	for( int istrip=striplo; istrip<=striphi; istrip++ ){
	  yADCsamples[isamp] += mod_data.ADCsamp_ystrips[istrip][isamp]*splitfraction[istrip];
	}
      }

      clust_data.ADCsamp_yclust.push_back( yADCsamples );
      
      clust_data.nclusty++;
    }
  }

  islocalmax.clear();
  localmaxima.clear();
  peakstoerase.clear();

  //Now we are done with 1D clustering. Form 2D clusters: for now we form all possible X/Y (or U/V) associations:

  //long n2Dclustercandidates_raw = clust_data.nclustx * clust_data.nclusty;

  for( int iclustx=0; iclustx<clust_data.nclustx; iclustx++ ){
    for( int iclusty=0; iclusty<clust_data.nclusty; iclusty++ ){
      //to calculate cluster-level correlation coefficients, we need to use cluster-summed ADC samples calculated above
      
      double sumx=0.0, sumy=0.0, sumx2=0.0, sumy2=0.0, sumxy=0.0;

      double sumxmax=0.0, sumymax=0.0,sumx2max=0.0,sumy2max=0.0, sumxymax=0.0; 
      
      for( int isamp=0; isamp<nADCsamples; isamp++ ){
	sumx += clust_data.ADCsamp_xclust[iclustx][isamp];
	sumy += clust_data.ADCsamp_yclust[iclusty][isamp];
	sumx2 += pow( clust_data.ADCsamp_xclust[iclustx][isamp], 2 );
	sumy2 += pow( clust_data.ADCsamp_yclust[iclusty][isamp], 2 );
	sumxy += clust_data.ADCsamp_xclust[iclustx][isamp] * clust_data.ADCsamp_yclust[iclusty][isamp];

	sumxmax += mod_data.ADCsamp_xstrips[clust_data.ixstripmax[iclustx]][isamp];
	sumymax += mod_data.ADCsamp_ystrips[clust_data.iystripmax[iclusty]][isamp];
	sumx2max += pow( mod_data.ADCsamp_xstrips[clust_data.ixstripmax[iclustx]][isamp], 2 );
	sumy2max += pow( mod_data.ADCsamp_ystrips[clust_data.iystripmax[iclusty]][isamp], 2 );
	sumxymax += mod_data.ADCsamp_xstrips[clust_data.ixstripmax[iclustx]][isamp]*
	  mod_data.ADCsamp_ystrips[clust_data.iystripmax[iclusty]][isamp];
      }

      double nSAMP = double(nADCsamples);

      double mux = sumx/nSAMP;
      double muy = sumy/nSAMP;
      double varx = sumx2/nSAMP - pow(mux,2);
      double vary = sumy2/nSAMP - pow(muy,2);
      double sigx = sqrt(varx);
      double sigy = sqrt(vary);
      
      double ccor_clust = (sumxy - nSAMP*mux*muy)/(nSAMP*sigx*sigy);

      double muxmax = sumxmax/nSAMP;
      double muymax = sumymax/nSAMP;
      double varxmax = sumx2max/nSAMP - pow(muxmax,2);
      double varymax = sumy2max/nSAMP - pow(muymax,2);
      double sigxmax = sqrt(varxmax);
      double sigymax = sqrt(varymax);

      double ccor_maxstrips = (sumxymax - nSAMP*muxmax*muymax)/(nSAMP*sigxmax*sigymax);

      //Do we want to apply timing and/or threshold and/or ADC asymmetry cuts here?
      //At this stage we keep everything. Later we may wish to revisit that:
      clust_data.itrack_clust2D.push_back( -1 );
      clust_data.ixclust2D.push_back( iclustx );
      clust_data.iyclust2D.push_back( iclusty );
      //these are redundant, but ok:
      clust_data.nstripx2D.push_back( clust_data.nstripx[iclustx] );
      clust_data.nstripy2D.push_back( clust_data.nstripy[iclusty] );
      clust_data.xclust2D.push_back( clust_data.xmean[iclustx] );
      clust_data.yclust2D.push_back( clust_data.ymean[iclusty] );
      clust_data.xclust2Dcorr.push_back( clust_data.xmean[iclustx] );
      clust_data.yclust2Dcorr.push_back( clust_data.ymean[iclusty] );
      clust_data.Eclust2D.push_back( 0.5*(clust_data.totalchargex[iclustx]+clust_data.totalchargey[iclusty]) );
      clust_data.dEclust2D.push_back( clust_data.totalchargex[iclustx]-clust_data.totalchargey[iclusty] );
      clust_data.tclust2D.push_back( 0.5*(clust_data.txmean[iclustx]+clust_data.tymean[iclusty]) );
      clust_data.tclust2Dwalkcorr.push_back( clust_data.tclust2D[clust_data.nclust2D] );
      clust_data.dtclust2D.push_back( clust_data.txmean[iclustx]-clust_data.tymean[iclusty] );
      clust_data.dtclust2Dwalkcorr.push_back( clust_data.dtclust2D[clust_data.nclust2D] );
      clust_data.CorrCoeff2D.push_back( ccor_clust );
      clust_data.CorrCoeffMaxStrips.push_back( ccor_maxstrips );
      clust_data.keepclust2D.push_back( true );

      //now calculate global hit coordinates:
      
      double Utemp = clust_data.xmean[iclustx];
      double Vtemp = clust_data.ymean[iclusty];

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

      double xstripmax = (clust_data.ixstripmax[iclustx] + 0.5 - 0.5*mod_nstripsu[module] )*mod_ustrip_pitch[module];
      double ystripmax = (clust_data.iystripmax[iclusty] + 0.5 - 0.5*mod_nstripsv[module] )*mod_vstrip_pitch[module];

      clust_data.xmom2D.push_back( (clust_data.xmean[iclustx]-xstripmax)/mod_ustrip_pitch[module] );
      clust_data.ymom2D.push_back( (clust_data.ymean[iclusty]-ystripmax)/mod_vstrip_pitch[module] );
      
      clust_data.nclust2D++;
    } 
  }

  if( clust_data.nclust2D > 0 ){
    prune_clusters( clust_data );
  }
  return 0;
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

void find_tracks_old( map<int,clusterdata_t> mod_clusters, trackdata_t &trackdata ){
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
	  if( trackdata.ntracks == 0 && nhitsrequired == layers_2Dmatch.size() ){

	    //cout << "too many hit combos, skipping tracking for this event..." << endl;
	    NSKIPPED++;
	  }
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
		
	      } //if (nhits >= nhitsrequired)
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

	    //Experimental: eliminate all other hits containing either the same X maximum or same Y maximum as a hit already used in a track.
	    for( int ihittemp=0; ihittemp<hitlisttemp.size(); ihittemp++ ){
	      int modtemp = modlisttemp[ihittemp];
	      int iclusttemp = hitlisttemp[ihittemp]; //position in 2D cluster array:
	      int ixtemp = mod_clusters[modtemp].ixclust2D[iclusttemp];
	      int iytemp = mod_clusters[modtemp].iyclust2D[iclusttemp];
	      int layertemp = mod_layer[modtemp];

	      for( int jhittemp=0; jhittemp<N2Dhits_layer[layertemp]; jhittemp++ ){
	    	int modj = modindexhit2D[layertemp][jhittemp];
	    	int clustj = clustindexhit2D[layertemp][jhittemp];
	    	if( modj == modtemp ){
	    	  int ixj = mod_clusters[modj].ixclust2D[clustj];
	    	  int iyj = mod_clusters[modj].iyclust2D[clustj];
	    	  if( ixj == ixtemp || iyj == iytemp ){
	    	    hitused2D[layertemp][jhittemp] = true;
	    	  }
	    	}
	      }
	      
	    }
	    
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

    bool skipped = true;
    
    trackdata.ntracks = 0;
    
    bool foundtrack=true;
    while( foundtrack ){

      //as long as we found a track on previous track-finding iteration, look for more tracks:
      foundtrack = false;

      int nhitsrequired=layers_2Dmatch.size(); //initially, we require nhits = total number of layers with unused hits:

      //If we don't find a good track at nhitsrequired hits, we can decrement the required number of hits as long as it is greater than the minimum
      // of TOTAL_REQUIRED_HIT:

      bool foundgoodpair = false;

     
      
      while( nhitsrequired >= TOTAL_REQUIRED_HIT ){ //first iteration: determine number of layers with (unused) hits, populate lists of unused hits, count number of combinations,
	//etc:
	
	foundtrack = false;
	foundgoodpair = false;
	
	//here we populate the lists of free hits by layer: if we found a track on the previous iteration, some hits will have been marked as used
	//if we didn't find a track, then no hits will have been marked as used, but the number of hits required to
	//make a track will have been decremented:
	long ncombosremaining=1;

	map<int,int> Nfreehits_layer; //free hit count mapped by layer
	set<int> layerswithfreehits;  //list of layers with free hits
	map<int,vector<int> > freehitlist_layer; //list of free hits mapped by layer: index in the unchanging arrays defined above:
	map<int,int> freehitcounter; //counter for looping over combos:

	//GRID:
	map<int,vector<int> > Nfreehits_binxy_layer;
	map<int,vector<vector<int> > > freehitlist_binxy_layer;
	
	//Populate the "free hit list" information:
	for(set<int>::iterator ilay=layers_2Dmatch.begin(); ilay!=layers_2Dmatch.end(); ++ilay ){
	  int layer = *ilay;
	  Nfreehits_layer[layer] = 0;

	  int nbins_gridxy = gridnbinsx_layer[layer]*gridnbinsy_layer[layer];
	  
	  Nfreehits_binxy_layer[layer].resize( nbins_gridxy );
	  freehitlist_binxy_layer[layer].resize( nbins_gridxy );
	  for( int bin=0; bin<nbins_gridxy; bin++ ){
	    Nfreehits_binxy_layer[layer][bin] = 0;
	    freehitlist_binxy_layer[layer][bin].clear();
	  }
	  
	  for( int ihit=0; ihit<N2Dhits_layer[layer]; ihit++ ){
	    if( !hitused2D[layer][ihit] ) {
	      Nfreehits_layer[layer]++;
	      freehitlist_layer[layer].push_back( ihit );

	      double xgtemp = mod_clusters[modindexhit2D[layer][ihit]].xglobal2D[clustindexhit2D[layer][ihit]];
	      double ygtemp = mod_clusters[modindexhit2D[layer][ihit]].yglobal2D[clustindexhit2D[layer][ihit]];

	      int binxtemp = int( (xgtemp - gridxmin_layer[layer])/gridxbinwidth );
	      int binytemp = int( (ygtemp - gridymin_layer[layer])/gridybinwidth );

	      if( binxtemp >= 0 && binxtemp < gridnbinsx_layer[layer] &&
		  binytemp >= 0 && binytemp < gridnbinsy_layer[layer] ){
		int binxytemp = binxtemp + gridnbinsx_layer[layer]*binytemp;

		Nfreehits_binxy_layer[layer][binxytemp]++;
		freehitlist_binxy_layer[layer][binxytemp].push_back( ihit );
		
	      }
	    }
	  }
 
	  if( Nfreehits_layer[layer] > 0 ){
	    ncombosremaining *= Nfreehits_layer[layer];
	    layerswithfreehits.insert(layer);
	    freehitcounter[layer] = 0;
	    // minlayer = (layer < minlayer ) ? layer : minlayer;
	    // maxlayer = (layer > maxlayer ) ? layer : maxlayer;
	  }
	}

	//this will get us stuck in an infinite loop if we don't find a track: comment out
	//nhitsrequired = layerswithfreehits.size();

	//if( ncombosremaining > maxnhitcombinations || ncombosremaining < 0 ){
	  //   cout << "too many hit combos, skipping tracking for this event..." << endl;

	  //   if( trackdata.ntracks == 0 && nhitsrequired == layers_2Dmatch.size() ) NSKIPPED++;
	  
	//   return;
	// }

	//	cout << "Starting tracking, ncombosremaining = " << ncombosremaining << endl;
	
	if( layerswithfreehits.size() >= nhitsrequired ){ //If the number of layers with free hits exceeds the number of hits required to make a track, proceed:
	  

	  //number of layers to in total:
	  int nlayerstot = layerswithfreehits.size();

	  map<int,int> hitcombo; //hits mapped by layer:
	  map<int,bool> ontrack; //flag to indicate whether this layer actually ended up on the track; to handle cases when the number of layers with free hits exceeds the number of layers that end up on the track
	  map<int,double> xresid_layer; //"U" residual (generalized "X") for test combos, mapped by layer
	  map<int,double> yresid_layer; //"V" residual (generalized "Y") for test combos, mapped by layer
	  //bool first = true; //flag to indicate first hit combo to be considered

	  map<int,int> besthitcombo; //map to store "BEST" combo according to Track chi2, mapped by layer
	  map<int,bool> onbesttrack; //flag to indicate whether this layer actually ended up on the track; to handle cases when the number of layers with free hits exceeds the number of layers that end up on the track; i.e., when nhitsrequired < n layers with free hits
	  map<int,double> bestxresid; //"U" residual (generalized "X") for best hit combo, mapped by LAYER
	  map<int,double> bestyresid; //"V" residual (generalized "Y") for best hit combo, mapped by LAYER
	  
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

	  foundgoodpair = false;
	  
	  for( int icombo=0; icombo<layercombos[nhitsrequired].size(); icombo++ ){ //loop over all possible combinations of nhitsrequired LAYERS with available hits:

	    //nextcomboexists = true;
	    
	    set<int> layerstotest; //populate the list of layers to be tested:

	    int minlayer = nlayers+1;
	    int maxlayer = -1;
	    
	    //int testlayer1=minlayer,testlayer2=maxlayer;

	    long ncombosremaining = 1;
	    
	    for( int ihit=0; ihit<nhitsrequired; ihit++ ){

	      int layeri = layercombos[nhitsrequired][icombo][ihit];

	      if( layerswithfreehits.find(layeri) != layerswithfreehits.end() ){
	      
		layerstotest.insert( layeri );
	      //at the beginning of each layer combination, reset the free hit counter to zero for each layer being tested:
		//freehitcounter[layeri] = 0; we don't actually need this for new algorithm:
		//figure out lowest and highest layers in this combo:
		minlayer = (layeri < minlayer ) ? layeri : minlayer;
		maxlayer = (layeri > maxlayer ) ? layeri : maxlayer;

		ncombosremaining *= Nfreehits_layer[layeri];
		
	      }
	    }

	    if( layerstotest.size() < nhitsrequired ) { //skip this combination of layers if there aren't enough fired layers to make a track
	      continue;
	    }	 
	    
	    long ncombos_minmax = Nfreehits_layer[minlayer]*Nfreehits_layer[maxlayer];

	    // long nhits_otherlayers = 0;
	    // for( set<int>::iterator ilay=layerstotest.begin(); ilay != layerstotest.end(); ++ilay ){
	    //   if( *ilay != minlayer && *ilay != maxlayer ){ 
	    // 	nhits_otherlayers += Nfreehits_layer[*ilay];
	    //   }
	    // }
	    
	    long ncombostocheck = ncombos_minmax;
	    
	    foundgoodpair = ncombostocheck >= 1 && ncombostocheck <= maxnhitcombinations;   
	    
	    if( !foundgoodpair ){
	      //if the outermost layers have too many hit combinations, find the combination of layers
	      // with the largest lever arm in z such that the number of combinations is less than the
	      //maximum:
	      
	      int maxdiff=0;
	      for( set<int>::iterator ilay=layerstotest.begin(); ilay!=layerstotest.end(); ++ilay ){
		int layeri = *ilay;
		for( set<int>::iterator jlay=ilay; jlay!=layerstotest.end(); ++jlay ){
		  int layerj = *jlay;
		  if( layerj > layeri ){
		    long ncombostest = Nfreehits_layer[layeri]*Nfreehits_layer[layerj];

		    // nhits_otherlayers = 0;
		    // for( set<int>::iterator klay=layerstotest.begin(); klay != layerstotest.end(); ++klay ){
		    //   if( *klay != layeri && *klay != layerj ){
		    // 	nhits_otherlayers += Nfreehits_layer[*klay];
		    //   }
		    // }

		    //ncombostest *= nhits_otherlayers;
		    
		    if( ncombostest <= maxnhitcombinations && ncombostest >= 1 &&
			layerj - layeri > maxdiff ){ //we found a good pair
		      foundgoodpair = true;
		      minlayer = layeri;
		      maxlayer = layerj;
		      maxdiff = layerj - layeri;
		    }
		  }
		}
	      }
	    }
	    //if we don't find any pair with an acceptably low number of hit combinations to test for the current layer combination
	    //we skip this layer combination:
	    if( !foundgoodpair ){  
	      continue;
	    } else { //If we consider ANY combination of layers, we did not skip attempted track finding for this event:
	      skipped = false;
	    }
	    // if( ncombos_minmax > maxnhitcombinations || ncombos_minmax < 0 ) {
	    //   cout << "warning, skipping layer combination due to too many hit combinations, icombo, minlayer, maxlayer, nhitsrequired, ncombos, ntracks = "
	    // 	   << icombo << ", " << minlayer << ", " << maxlayer << ", " << nhitsrequired << ", " << ncombos_minmax << ", " << trackdata.ntracks << endl ;
	    //   continue;
	    // }
	    
	    double sumxhits=0.0,sumyhits=0.0,sumzhits=0.0,sumxzhits=0.0,sumyzhits=0.0,sumz2hits=0.0;
	    
	    double xtrtemp,ytrtemp,ztrtemp,xptrtemp,yptrtemp; //temporary storage for track parameters: 
	    
	    double varx,vary,varxp,varyp,covxxp,covyyp;
	    
	     
	    //Now loop over all combinations of one unused hit from minlayer and one unused hit from maxlayer:
	    for( int ihit=0; ihit<freehitlist_layer[minlayer].size(); ihit++ ){
	      for( int jhit=0; jhit<freehitlist_layer[maxlayer].size(); jhit++ ){

		//clear certain arrays before populating each combo:
		ontrack.clear();
		hitcombo.clear();
		xresid_layer.clear();
		yresid_layer.clear(); //at the beginning of each combination from inner most and outermost layers
		
		int hitmin=freehitlist_layer[minlayer][ihit];
		int hitmax=freehitlist_layer[maxlayer][jhit];

		//Next step is to calculate the straight line in 3D space connecting these two hits; then search for the hit in each intermediate layer closest to
		//the straight-line projection of the min and max layers:
		int modmin = modindexhit2D[minlayer][hitmin];
		int modmax = modindexhit2D[maxlayer][hitmax];
		int iclustmin = clustindexhit2D[minlayer][hitmin];
		int iclustmax = clustindexhit2D[maxlayer][hitmax];

		double xmin = mod_clusters[modmin].xglobal2D[iclustmin];		
		double ymin = mod_clusters[modmin].yglobal2D[iclustmin];		
		double zmin = mod_clusters[modmin].zglobal2D[iclustmin];
		
		double ymax = mod_clusters[modmax].yglobal2D[iclustmax];
		double xmax = mod_clusters[modmax].xglobal2D[iclustmax];
		double zmax = mod_clusters[modmax].zglobal2D[iclustmax];

		xptrtemp = (xmax-xmin)/(zmax-zmin);
		yptrtemp = (ymax-ymin)/(zmax-zmin);
		xtrtemp = 0.5*( xmin - zmin * xptrtemp + xmax - zmax * xptrtemp );
		ytrtemp = 0.5*( ymin - zmin * yptrtemp + ymax - zmax * yptrtemp );

		sumxhits = xmin + xmax;
		sumyhits = ymin + ymax;
		sumzhits = zmin + zmax;
		sumxzhits = xmin*zmin + xmax*zmax;
		sumyzhits = ymin*zmin + ymax*zmax;
		sumz2hits = pow(zmin,2)+pow(zmax,2);
		  
		nhits = 2;

		ontrack[minlayer] = true;
		ontrack[maxlayer] = true;

		hitcombo[minlayer] = hitmin;
		hitcombo[maxlayer] = hitmax;
		  
		//Now loop over all intermediate layers and find the hit which gives the smallest chi^2 when combined with the current layers on the track in a linear fit:

		//only test other layers for hit combinations that pass the track slope cuts:
		if( abs( xptrtemp ) <= TrackMaxSlopeX && abs(yptrtemp) <= TrackMaxSlopeY ){
		
		  //for( int layer=minlayer+1; layer<maxlayer; layer++ ){
		  for( set<int>::iterator ilay=layerstotest.begin(); ilay != layerstotest.end(); ++ilay ){

		    int layer = *ilay;
		
		    //if( layerstotest.find(layer)!=layerstotest.end() ){ //if this layer is in the list of layers we are currently considering:
		    if( layer != maxlayer && layer != minlayer ){
		      int khit_best=-1;
		      double minresid2 = 0.0;

		      ontrack[layer] = true; //initialize ontrack for this layer to false in order to avoid undefined behavior
		    
		      int modbest, hitbest=-1;
		      double xbest,ybest,zbest;

		      //now project the temporary track to this layer, and find the appropriate bin(s) in the grid to look for hits:
		      int nearest_module = -1;
		      double minr2=0.;
		      for( set<int>::iterator imod=modlist_layer[layer].begin(); imod!=modlist_layer[layer].end(); ++imod ){
			double trackxproj = xtrtemp + xptrtemp * mod_z0[*imod];
			double trackyproj = ytrtemp + yptrtemp * mod_z0[*imod];
			double r2 = pow( trackxproj - mod_x0[*imod], 2 ) + pow( trackyproj - mod_y0[*imod], 2 );

			if( nearest_module < 0 || r2 < minr2 ){
			  nearest_module = *imod;
			  minr2 = r2;
			}
		      }

		      //this exercise is only for the purpose of picking the appropriate grid bin within the layer:
		      double ztest = mod_z0[nearest_module];
		      double xtest = xtrtemp + xptrtemp*ztest;
		      double ytest = ytrtemp + yptrtemp*ztest;

		      int binxtest = int( (xtest-gridxmin_layer[layer])/gridxbinwidth );
		      int binytest = int( (ytest-gridymin_layer[layer])/gridybinwidth );

		      if( binxtest >= 0 && binxtest < gridnbinsx_layer[layer] && 
			  binytest >= 0 && binytest < gridnbinsy_layer[layer] ){ //track projection is within limits at this layer:
		      
			//if x and/or y are close to the upper or lower edge of a bin, also check neighboring bins:
			
			double binxdiff = (xtest - (gridxmin_layer[layer]+binxtest*gridxbinwidth) );
			double binydiff = (ytest - (gridymin_layer[layer]+binytest*gridybinwidth) );
			
			int binxtest_lo = binxtest, binxtest_hi = binxtest;
			int binytest_lo = binytest, binytest_hi = binytest;

			//If we are within 3 mm of the edge of the bin, test neighboring bins as well:
			if( binxdiff < 3.0 && binxtest > 0 ) binxtest_lo = binxtest-1;
			if( gridxbinwidth - binxdiff < 3.0 && binxtest + 1 < gridnbinsx_layer[layer] ) binxtest_hi = binxtest+1;
			if( binydiff < 3.0 && binytest > 0 ) binytest_lo = binytest-1;
			if( gridybinwidth - binydiff < 3.0 && binytest + 1 < gridnbinsy_layer[layer] ) binytest_hi = binytest+1;


			for( int binx=binxtest_lo; binx<=binxtest_hi; binx++ ){
			  for( int biny=binytest_lo; biny<=binytest_hi; biny++ ){
			    int binxy_test = binx + gridnbinsx_layer[layer]*biny;

			    for( int khit=0; khit<freehitlist_binxy_layer[layer][binxy_test].size(); khit++ ){
			      int hitk = freehitlist_binxy_layer[layer][binxy_test][khit];
			      int modk = modindexhit2D[layer][hitk];
			      int clustk = clustindexhit2D[layer][hitk];

			      double xhitk = mod_clusters[modk].xglobal2D[clustk];
			      double yhitk = mod_clusters[modk].yglobal2D[clustk];
			      double zhitk = mod_clusters[modk].zglobal2D[clustk];
			      
			      double uhitk = mod_clusters[modk].xclust2Dcorr[clustk];
			      double vhitk = mod_clusters[modk].yclust2Dcorr[clustk];

			      double trackxproj = xtrtemp + xptrtemp * zhitk;
			      double trackyproj = ytrtemp + yptrtemp * zhitk;
			      TVector3 trackpos_global( trackxproj, trackyproj, zhitk );
			      TVector3 modcenter_global( mod_x0[modk], mod_y0[modk], mod_z0[modk] );

			      TVector3 trackpos_local = mod_Rotinv[modk] * (trackpos_global - modcenter_global);

			      double utrack_k = trackpos_local.X()*mod_Pxu[modk] + trackpos_local.Y()*mod_Pyu[modk];
			      double vtrack_k = trackpos_local.X()*mod_Pxv[modk] + trackpos_local.Y()*mod_Pyv[modk];
			
			      double resid2hit = pow( uhitk - utrack_k, 2 ) + pow( vhitk - vtrack_k, 2 ); //

			      if( hitbest < 0 || resid2hit < minresid2 ){
				hitbest = hitk;
				modbest = modk;
				xbest = xhitk;
				ybest = yhitk;
				zbest = zhitk;
			      }
			      
			    } //end loop over unused hits in grid XY bin
			    
			  } //end loop over (up to two) ybins in the grid
			} //end loop over (up to two) xbins in the grid
		      } //end check whether nearest x and y bins are within layer limits
		      
		      // for( int khit=0; khit<freehitlist_layer[layer].size(); khit++ ){
		      // 	int hitk = freehitlist_layer[layer][khit];
		      // 	int modk = modindexhit2D[layer][hitk];
		      // 	int clustk = clustindexhit2D[layer][hitk];

		      // 	double xhitk = mod_clusters[modk].xglobal2D[clustk];
		      // 	double yhitk = mod_clusters[modk].yglobal2D[clustk];
		      // 	double zhitk = mod_clusters[modk].zglobal2D[clustk];

		      // 	double uhitk = mod_clusters[modk].xclust2Dcorr[clustk];
		      // 	double vhitk = mod_clusters[modk].yclust2Dcorr[clustk];

		      // 	double trackxproj = xtrtemp + xptrtemp * zhitk;
		      // 	double trackyproj = ytrtemp + yptrtemp * zhitk;
		      // 	TVector3 trackpos_global( trackxproj, trackyproj, zhitk );
		      // 	TVector3 modcenter_global( mod_x0[modk], mod_y0[modk], mod_z0[modk] );

		      // 	TVector3 trackpos_local = mod_Rotinv[modk] * (trackpos_global - modcenter_global);

		      // 	double utrack_k = trackpos_local.X()*mod_Pxu[modk] + trackpos_local.Y()*mod_Pyu[modk];
		      // 	double vtrack_k = trackpos_local.X()*mod_Pxv[modk] + trackpos_local.Y()*mod_Pyv[modk];
			
		      // 	double resid2hit = pow( uhitk - utrack_k, 2 ) + pow( vhitk - vtrack_k, 2 ); //

		      // 	if( khit_best < 0 || resid2hit < minresid2 ){ //either this is the first hit we are considering in this layer or it has the smallest
		      // 	  //squared residual with the current track:
		      // 	  khit_best = hitk; //use index in the free hit list for convenience. We may want to revisit this choice
		      // 	  minresid2 = resid2hit;

		      // 	  modbest = modk;
		      // 	  hitbest = hitk;
		      // 	  xbest = xhitk;
		      // 	  ybest = yhitk;
		      // 	  zbest = zhitk;
			
		      // 	}
		      // } //end loop over unused hits in each intermediate layer
		    
		      //if the closest hit in this layer to the current straight-line projection of the track is on a good track
		      double ndf_temp = double(2*(nhits+1)-4);
		      //chi2 contribution of this hit wrt the straight-line fit to all previous hits is resid2/sigma^2
		      // chi2/ndf < TrackChi2Cut 
		      if( hitbest >= 0 ){
			//this might be a good track: add this hit to the combo and update the linear fit:
			ontrack[layer] = true;
			hitcombo[layer] = hitbest;

			//update sums:
			sumxhits += xbest;
			sumyhits += ybest;
			sumzhits += zbest;
			sumxzhits += xbest*zbest;
			sumyzhits += ybest*zbest;
			sumz2hits += pow(zbest,2);

			nhits++;

			//update sums,
			// but don't update track parameters until chi2 calculation
			//update best track based on fit to all hits in current combo:
			// double denom = (sumz2hits*nhits - pow(sumzhits,2));
			// xptrtemp = (nhits*sumxzhits - sumxhits*sumzhits)/denom;
			// yptrtemp = (nhits*sumyzhits - sumyhits*sumzhits)/denom;
			// xtrtemp = (sumz2hits*sumxhits - sumzhits*sumxzhits)/denom;
			// ytrtemp = (sumz2hits*sumyhits - sumzhits*sumyzhits)/denom;
			
		      } //check if we found an unused hit sufficiently close to the track at this layer
		    } //test if current intermediate layer is in the list of layers to test
		  } //end loop over intermediate layers
		} //end check if the slope of the track formed from the min and max hits is within the allowed limits
		
		
		
		//
		if( nhits >= nhitsrequired ){ //calculate chi2, and update best hit combo as applicable:
		  //update best track based on fit to all hits in current combo:
		  double denom = (sumz2hits*nhits - pow(sumzhits,2));
		  xptrtemp = (nhits*sumxzhits - sumxhits*sumzhits)/denom;
		  yptrtemp = (nhits*sumyzhits - sumyhits*sumzhits)/denom;
		  xtrtemp = (sumz2hits*sumxhits - sumzhits*sumxzhits)/denom;
		  ytrtemp = (sumz2hits*sumyhits - sumzhits*sumyzhits)/denom;

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
		      
		    } //end check if layer on track
		  } //end loop over layers for chi2 and residual calculation
		
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

		    ngoodcombos++;
		  }
		} //check number of hits on track candidate at least nhitsrequired
	      } //end loop over free hits in max layer (topmost layer)
	    } //end loop over free hits in min layer (bottom-most layer)
	  } //end loop over all possible combinations of nhitsrequired layers
	  

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

	    //Experimental: eliminate all other hits containing either the same X maximum or same Y maximum as a hit already used in a track.
	    for( int ihittemp=0; ihittemp<hitlisttemp.size(); ihittemp++ ){
	      int modtemp = modlisttemp[ihittemp];
	      int iclusttemp = hitlisttemp[ihittemp]; //position in 2D cluster array:
	      int ixtemp = mod_clusters[modtemp].ixclust2D[iclusttemp];
	      int iytemp = mod_clusters[modtemp].iyclust2D[iclusttemp];
	      int layertemp = mod_layer[modtemp];

	      for( int jhittemp=0; jhittemp<N2Dhits_layer[layertemp]; jhittemp++ ){
	    	int modj = modindexhit2D[layertemp][jhittemp];
	    	int clustj = clustindexhit2D[layertemp][jhittemp];
	    	if( modj == modtemp ){
	    	  int ixj = mod_clusters[modj].ixclust2D[clustj];
	    	  int iyj = mod_clusters[modj].iyclust2D[clustj];
	    	  if( ixj == ixtemp || iyj == iytemp ){
	    	    hitused2D[layertemp][jhittemp] = true;
	    	  }
	    	}
	      }
	      
	    }
	    
	    trackdata.ntracks++;
	  } //end check on number of good combinations > 0 and track with smallest chi^2 passes cut
	  
	} //end check on "layers with free hits >= nhits required"
	  //if( !foundtrack ) break; //prevents getting stuck in infinite loop if we don't find any tracks on first iteration:

	//if we fail to find a track at the current hit requirement; we reduce the number of hits required:
	//If the number of hits required falls below 3, we exit the loop:
	//If we DO find a track at the current hit requirement, we look for more tracks at the current
	//hit requirement before dropping the hit requirement.
	
	if( !foundtrack ) {
	  //if( nhitsrequired == TOTAL_REQUIRED_HIT && !foundgoodpair && trackdata.ntracks == 0 ) NSKIPPED++; 
	  nhitsrequired--;
	}
	
      } //while( nhitsrequired >= TOTAL_REQUIRED_HIT )
	    
	//cout << "nhitsrequired = " << nhitsrequired << endl;
    } //while( foundtrack ) //If we found a track on this try, look for more tracks!
    if( skipped ) NSKIPPED++;
  } //total number of layers with at least one 2D cluster >= TOTAL_REQUIRED_HIT
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

void GEM_reconstruct_standalone_consolidated( const char *filename, const char *configfilename, const char *outfilename="temp.root" ){
  auto program_start = high_resolution_clock::now(); //starting time of the program --WX
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

  long firstevent=0;
  
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

	  if( skey == "mod_Xgain" && ntokens >= 3 ){
	    TString smod = ( (TObjString*) (*tokens)[1] )->GetString();

	    int module = smod.Atoi();

	    TString snAPVs = ( (TObjString*) (*tokens)[2] )->GetString();

	    int nAPVs = snAPVs.Atoi();

	    if( ntokens - 3 >= nAPVs ){ //read in gain coefficients:
	      for( int iAPV=0; iAPV<nAPVs; iAPV++ ){
		TString sgain = ( (TObjString*) (*tokens)[iAPV+3] )->GetString();

		mod_Xgain[module].push_back(sgain.Atof());
	      }
	    }
	  }

	  if( skey == "mod_Ygain" && ntokens >= 3 ){
	    TString smod = ( (TObjString*) (*tokens)[1] )->GetString();

	    int module = smod.Atoi();

	    TString snAPVs = ( (TObjString*) (*tokens)[2] )->GetString();

	    int nAPVs = snAPVs.Atoi();

	    if( ntokens - 3 >= nAPVs ){ //read in gain coefficients:
	      for( int iAPV=0; iAPV<nAPVs; iAPV++ ){
		TString sgain = ( (TObjString*) (*tokens)[iAPV+3] )->GetString();

		mod_Ygain[module].push_back(sgain.Atof());
	      }
	    }
	  }
	  
	  if( skey == "gridxbinwidth" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    gridxbinwidth = skey.Atof();
	  }

	  if( skey == "gridybinwidth" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    gridybinwidth = skey.Atof();
	  }
	  
	  if( skey == "firstevent" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    firstevent = stemp.Atoi();
	  }
	  
	  if( skey == "maxneighborsx" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    //maxnhitcombinations = stemp.Atoi();
	    maxneighborsX = stemp.Atoi();
	  }

	   if( skey == "maxneighborsy" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    //maxnhitcombinations = stemp.Atoi();
	    maxneighborsY = stemp.Atoi();
	   }
	  
	  if( skey == "maxhitcombos" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    maxnhitcombinations = stemp.Atoi();
	  }
	  
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
	tokens->Delete();
      }
    }
  }

  // //PulseShape->FixParameter(2,tau_pulseshape);
  // PulseShape->SetParameter(2,tau_pulseshape);
  // //PulseShape->SetParLimits(1,-150,150);
  // PulseShape->SetParLimits(0,0.0,1.e5);
  // // PulseShape->SetParLimits(2,0.0, 300.0);

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


  //Check if module gain has been initialized: 
  
  for( int imod=0; imod<nmodules; imod++ ){
    int layer = mod_layer[imod];
    nmod_layer[layer]++;
    zavg_layer[layer] += mod_z0[imod];
    for (unsigned int j=0; j<gridContainer.size(); j++){
        if (gridContainer[j].layer == layer) gridContainer[j].AddModule(imod);
    }

    nmodules_layer[layer]++;
    modlist_layer[layer].insert(imod);

    //The following lines insure that mod_Xgain and mod_Ygain have been properly initialized
    //if the user has not provided their own values:
    
    mod_nAPVX[imod] = mod_nstripsu[imod]/128;
    if( mod_nstripsu[imod]%128 > 0 ) mod_nAPVX[imod]++;
   
    mod_nAPVY[imod] = mod_nstripsv[imod]/128;
    if( mod_nstripsv[imod]%128 > 0 ) mod_nAPVY[imod]++;

    if( mod_Xgain.find(imod) == mod_Xgain.end() ){ //X gain for this module has NOT been initialized: set all to 1
      for( int i=0; i<mod_nAPVX[imod]; i++ ){
	mod_Xgain[imod].push_back( 1.0 );
      }
    } else {
      if( mod_Xgain[imod].size() < mod_nAPVX[imod] ){
	for( int i=0; i<mod_nAPVX[imod] - mod_Xgain[imod].size(); i++ ){
	  mod_Xgain[imod].push_back( 1.0 );
	}
      }
    }

    if( mod_Ygain.find(imod) == mod_Ygain.end() ){ //X gain for this module has NOT been initialized: set all to 1
      for( int i=0; i<mod_nAPVY[imod]; i++ ){
	mod_Ygain[imod].push_back( 1.0 );
      }
    } else {
      if( mod_Ygain[imod].size() < mod_nAPVY[imod] ){
	for( int i=0; i<mod_nAPVY[imod] - mod_Ygain[imod].size(); i++ ){
	  mod_Ygain[imod].push_back( 1.0 );
	}
      }
    }
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
  int HitXstripMax[nlayers];
  int HitYstripMax[nlayers];
  int HitXstripLo[nlayers];
  int HitYstripLo[nlayers];
  int HitXstripHi[nlayers];
  int HitYstripHi[nlayers];

  Tout->Branch("EventID",&EventID,"EventID/I");
  Tout->Branch("Ntracks",&Ntracks,"Ntracks/I");
  Tout->Branch("CALOsum",&CALOsum,"CALOsum/D");
  Tout->Branch("NGOODSCINT", &NGOODSCINT, "NGOODSCINT/I");
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

  //Add strip index information for each hit:
  Tout->Branch("HitXstripMax",HitXstripMax,"HitXstripMax[TrackNhits]/I");
  Tout->Branch("HitYstripMax",HitYstripMax,"HitYstripMax[TrackNhits]/I");

  Tout->Branch("HitXstripLo",HitXstripLo,"HitXstripLo[TrackNhits]/I");
  Tout->Branch("HitYstripLo",HitYstripLo,"HitYstripLo[TrackNhits]/I");

  Tout->Branch("HitXstripHi",HitXstripHi,"HitXstripHi[TrackNhits]/I");
  Tout->Branch("HitYstripHi",HitYstripHi,"HitYstripHi[TrackNhits]/I");
  
  double xgmin_all=1e9, xgmax_all=-1e9, ygmin_all=1e9, ygmax_all=-1e9;
  //  map<int,double> xgmin_layer, xgmax_layer, ygmin_layer, ygmax_layer;

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

  //calculate also grid for tracking speed improvement:
  for( int layer=0; layer<nlayers; layer++ ){
    double xstart = xgmin_layer[layer]-gridxbinwidth;
    int nbinsx=0;
    //double xstop = xstart + nbinsx*gridxbinwidth;
    while( xstart + nbinsx*gridxbinwidth <= xgmax_layer[layer]+gridxbinwidth ){ 
      nbinsx++;
    }

    double xstop = xstart + nbinsx*gridxbinwidth;

    gridnbinsx_layer[layer] = nbinsx;
    gridxmin_layer[layer] = xstart;
    gridxmax_layer[layer] = xstop;
    
    double ystart = ygmin_layer[layer]-gridybinwidth;
    int nbinsy=0;
    while( ystart + nbinsy*gridybinwidth <= ygmax_layer[layer]+gridybinwidth ){
      nbinsy++;
    }

    double ystop = ystart + nbinsy*gridybinwidth;

    gridnbinsy_layer[layer] = nbinsy;
    gridymin_layer[layer] = ystart;
    gridymax_layer[layer] = ystop;
    
  }
  
  //Is the size of this canvas biased against people with low-resolution monitors? Probably yes.
  TCanvas *c1 = new TCanvas("c1","c1",1600,1500);
  // Let's retool this event display; instead, let's divide the canvas by layers and modules instead of by layer alone:  
  c1->Divide( maxnmodperlayer,nlayers,.001,.001);

  //Now let's create a histogram to display the strip ADC histos event by event:
  TCanvas *c2x = new TCanvas("c2x","c2x",2400,1500); 
  c2x->Divide( maxnmodperlayer, nlayers, .001,.001);

  //Now let's create a histogram to display the strip ADC histos event by event:
  TCanvas *c2y = new TCanvas("c2y","c2y",2400,1500); 
  c2y->Divide( maxnmodperlayer, nlayers, .001,.001);
  
  TCanvas *c_proj = new TCanvas("c1_proj","c1_proj",1200,800);
  c_proj->Divide(2,1);
  
  //  TCanvas *c2 = new TCanvas("c2","c2",1600,1200);
  
  //TClonesArray *hframe_layers = new TClonesArray("TH2D",nlayers);

  TClonesArray *hframe_modules = new TClonesArray("TH2D",nmodules);
  
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

  //Now suppose we want to plot the strip signals (U and V) by module to debug clustering and cluster-splitting:


  //These are for event display only. We need different ones for software gain-match:
  //Histograms for ALL strips:
  TClonesArray *hmod_ADC_vs_Xstrip = new TClonesArray("TH1D",nmodules);
  TClonesArray *hmod_ADC_vs_Ystrip = new TClonesArray("TH1D",nmodules);
  //Histograms for strips included in good 2D clusters:
  TClonesArray *hmod_ADC_vs_Xstrip_good = new TClonesArray("TH1D",nmodules);
  TClonesArray *hmod_ADC_vs_Ystrip_good = new TClonesArray("TH1D",nmodules);

  //Here we can also make ADC distribution plots for "gain matching" analysis:

  TClonesArray *hmod_ADCmax_vs_APV_Xstrips = new TClonesArray("TH2D",nmodules);
  TClonesArray *hmod_ADCmax_vs_APV_Ystrips = new TClonesArray("TH2D",nmodules);

  TClonesArray *hmod_ADCsum_vs_APV_Xstrips = new TClonesArray("TH2D",nmodules);
  TClonesArray *hmod_ADCsum_vs_APV_Ystrips = new TClonesArray("TH2D",nmodules);

  TClonesArray *hmod_ADCmax_vs_APV_Xstrips_all = new TClonesArray("TH2D",nmodules);
  TClonesArray *hmod_ADCmax_vs_APV_Ystrips_all = new TClonesArray("TH2D",nmodules);

  TClonesArray *hmod_ADCsum_vs_APV_Xstrips_all = new TClonesArray("TH2D",nmodules);
  TClonesArray *hmod_ADCsum_vs_APV_Ystrips_all = new TClonesArray("TH2D",nmodules);

  //Now let's also collect the cluster charges by APV:
  TClonesArray *hmod_clustersumX_vs_APV = new TClonesArray("TH2D",nmodules);
  TClonesArray *hmod_clustersumY_vs_APV = new TClonesArray("TH2D",nmodules);
  
  for( int imod=0; imod<nmodules; imod++ ){
    TString histname;
    histname.Form("hmod_ADC_vs_Xstrip_module%d",imod);
    TString histtitle; 
    histtitle.Form("X strips, Module %d", imod );

    int nbins = mod_nstripsu[imod];
    double min=-0.5, max=nbins-0.5;
    
    new( (*hmod_ADC_vs_Xstrip)[imod] ) TH1D( histname.Data(), histtitle.Data(), nbins, min, max );

    histname.Form("hmod_ADC_vs_Ystrip_module%d",imod);
    histtitle.Form("Y strips, Module %d", imod );

    nbins = mod_nstripsv[imod];
    min = -0.5;
    max = nbins-0.5;
    
    new( (*hmod_ADC_vs_Ystrip)[imod] ) TH1D( histname.Data(), histtitle.Data(), nbins, min, max );

    //now repeat for the "good" clusters:
    histname.Form("hmod_ADC_vs_Xstrip_good_module%d", imod );
    histtitle.Form("X strips, Module %d, good clusters", imod );

    nbins = mod_nstripsu[imod];
    min=-0.5;
    max=nbins-0.5;
    new( (*hmod_ADC_vs_Xstrip_good)[imod] ) TH1D( histname.Data(), histtitle.Data(), nbins, min, max );

    //now repeat for the "good" clusters:
    histname.Form("hmod_ADC_vs_Ystrip_good_module%d", imod );
    histtitle.Form("Y strips, Module %d, good clusters", imod );

    nbins = mod_nstripsv[imod];
    min = -0.5;
    max = nbins-0.5;

    new( (*hmod_ADC_vs_Ystrip_good)[imod] ) TH1D( histname.Data(), histtitle.Data(), nbins, min, max );

    histname.Form( "hmod_ADCmax_vs_APV_Xstrips_module%d", imod );
    nbins = mod_nstripsu[imod]/128;

    new( (*hmod_ADCmax_vs_APV_Xstrips)[imod] ) TH2D( histname.Data(), "Max ADC in max X strip on cluster in good track", nbins, -0.5, nbins+0.5, 256, 0.0, 4096.0 );

    histname.Form( "hmod_ADCmax_vs_APV_Ystrips_module%d", imod );
    nbins = mod_nstripsv[imod]/128;

    new( (*hmod_ADCmax_vs_APV_Ystrips)[imod] ) TH2D( histname.Data(), "Max ADC in max Y strip on cluster in good track", nbins, -0.5, nbins+0.5, 256, 0.0, 4096.0 );



    histname.Form( "hmod_ADCsum_vs_APV_Xstrips_module%d", imod );
    nbins = mod_nstripsu[imod]/128;

    new( (*hmod_ADCsum_vs_APV_Xstrips)[imod] ) TH2D( histname.Data(), "ADC sum in max X strip on cluster in good track", nbins, -0.5, nbins+0.5, 256, 0.0, 1.5e4 );

    histname.Form( "hmod_ADCsum_vs_APV_Ystrips_module%d", imod );
    nbins = mod_nstripsv[imod]/128;

    new( (*hmod_ADCsum_vs_APV_Ystrips)[imod] ) TH2D( histname.Data(), "ADC sum in max Y strip on cluster in goood track", nbins, -0.5, nbins+0.5, 256, 0.0, 1.5e4 );


    histname.Form( "hmod_ADCmax_vs_APV_Xstrips_all_module%d", imod );
    nbins = mod_nstripsu[imod]/128;

    new( (*hmod_ADCmax_vs_APV_Xstrips_all)[imod] ) TH2D( histname.Data(), "Max ADC in any X strip on cluster in good track", nbins, -0.5, nbins+0.5, 256, 0.0, 4096.0 );

    histname.Form( "hmod_ADCmax_vs_APV_Ystrips_all_module%d", imod );
    nbins = mod_nstripsv[imod]/128;

    new( (*hmod_ADCmax_vs_APV_Ystrips_all)[imod] ) TH2D( histname.Data(), "Max ADC in any Y strip on cluster in good track", nbins, -0.5, nbins+0.5, 256, 0.0, 4096.0 );



    histname.Form( "hmod_ADCsum_vs_APV_Xstrips_all_module%d", imod );
    nbins = mod_nstripsu[imod]/128;

    new( (*hmod_ADCsum_vs_APV_Xstrips_all)[imod] ) TH2D( histname.Data(), "ADC sum in any X strip on cluster in good track", nbins, -0.5, nbins+0.5, 256, 0.0, 1.5e4 );

    histname.Form( "hmod_ADCsum_vs_APV_Ystrips_all_module%d", imod );
    nbins = mod_nstripsv[imod]/128;

    new( (*hmod_ADCsum_vs_APV_Ystrips_all)[imod] ) TH2D( histname.Data(), "ADC sum in any Y strip on cluster in goood track", nbins, -0.5, nbins+0.5, 256, 0.0, 1.5e4 );

    histname.Form( "hmod_clustersumX_vs_APV_module%d", imod );

    nbins = mod_nstripsu[imod]/128;
    new( (*hmod_clustersumX_vs_APV)[imod] ) TH2D( histname.Data(), "Cluster sum in X strips, cluster on good track", nbins, -0.5, nbins+0.5, 256, 0.0, 5e4 );

    histname.Form( "hmod_clustersumY_vs_APV_module%d", imod );
    nbins = mod_nstripsv[imod]/128;

    new( (*hmod_clustersumY_vs_APV)[imod] ) TH2D( histname.Data(), "Cluster sum in Y strips, cluster on good track", nbins, -0.5, nbins+0.5, 256, 0.0, 5e4 );
    
  }
  
  TChain *C = new TChain("GEMHit");

  C->Add(filename);

  // int nch;
  // int *strip;
  // int *adc0;

  // C->SetBranchAddress("nch",&nch);
  // C->SetBranchAddress("strip",strip);
  // C->SetBranchAddress("adc0",adc0);

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
  
  long nevent=firstevent;

  fout->cd();

  
  TH2D *hCALfadc_vs_ch = new TH2D("hCALfadc_vs_ch","",32,-0.5,31.5,1000,0.0,25000.0);
  TH2D *hCALftime_vs_ch = new TH2D("hCALftime_vs_ch","",32,-0.5,31.5,51,-0.5,50.5);
  TH1D *hCALfadc_sum = new TH1D("hCALfadc_sum","",2000,0.0,50000.0);

  //TH1D *hCALfadc_sum_nocuts = new TH1D("hCALfadc_sum_nocuts","",2000,0.0,50000.0);
  // TH1D *hCALfadc_shsum = new TH1D("hCALfadc_shsum","",2000,0.0,50000.0);
  // TH1D *hCALfadc_pssum = new TH1D("hCALfadc_pssum","",2000,0.0,50000.0);
  //TH2D *hCALfadc_sh_vs_ps = new TH2D("hCALfadc_sh_vs_ps","",250,0.0,50000.0,250,0.0,50000.0);

  //Scintillator info:
  TH2D *hSCINT_tdc_vs_ch = new TH2D("hSCINT_tdc_vs_ch", "", 25,-0.5,24.5, 600,0.0,150.0 );
  TH1D *hNgoodScint = new TH1D("hNgoodScint", "N scint. TDC hits", 11, -0.5, 10.5 ); 


  TH2D *hCALfadc_vs_ch_goodtrack = new TH2D("hCALfadc_vs_ch_goodtrack","",32,-0.5,31.5,1000,0.0,25000.0);
  TH2D *hCALftime_vs_ch_goodtrack = new TH2D("hCALftime_vs_ch_goodtrack","",32,-0.5,31.5,51,-0.5,50.5);
  TH1D *hCALfadc_sum_goodtrack = new TH1D("hCALfadc_sum_goodtrack","",2000,0.0,50000.0);

  // TH1D *hCALfadc_sum_nocuts = new TH1D("hCALfadc_sum_nocuts","",2000,0.0,50000.0);
  // TH1D *hCALfadc_shsum = new TH1D("hCALfadc_shsum","",2000,0.0,50000.0);
  // TH1D *hCALfadc_pssum = new TH1D("hCALfadc_pssum","",2000,0.0,50000.0);
  //TH2D *hCALfadc_sh_vs_ps = new TH2D("hCALfadc_sh_vs_ps","",250,0.0,50000.0,250,0.0,50000.0);

  //Scintillator info:
  TH2D *hSCINT_tdc_vs_ch_goodtrack = new TH2D("hSCINT_tdc_vs_ch_goodtrack", "", 25,-0.5,24.5, 600, 0.0, 150.0);
  TH1D *hNgoodScint_goodtrack = new TH1D("hNgoodScint_goodtrack", "N scint. TDC hits", 11, -0.5, 10.5 ); 

  TH1D *hNtracks_nocuts = new TH1D("hNtracks_nocuts", "Number of tracks/trigger", 11, -0.5, 10.5 );
  TH1D *hNtracks_GoodScintAndCALO = new TH1D("hNtracks_GoodScintAndCALO", "Number of tracks/trigger with good scint. and CALO hits", 11, -0.5, 10.5);
  
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
  
  TH1D *hclustwidth_x = new TH1D("hclustwidth_x","Width of X clusters in strips",100,0.5,100.5);
  TH1D *hclustwidth_y = new TH1D("hclustwidth_y","Width of Y clusters in strips",100,0.5,100.5);
  TH1D *hclustADCsum_x = new TH1D("hclustADCsum_x","ADC sum in X strips",200,0.0,100000.0);
  TH1D *hclustADCsum_y = new TH1D("hclustADCsum_y","ADC sum in Y strips",200,0.0,100000.0);

  TH2D *hclustwidth_x_module = new TH2D("hclustwidth_x_module","Width of clusters in X strips",maxnstripXpercluster,0.5,maxnstripXpercluster+0.5,nmodules,-0.5,nmodules-0.5);
  TH2D *hclustwidth_y_module = new TH2D("hclustwidth_y_module","Width of clusters in Y strips",maxnstripYpercluster,0.5,maxnstripYpercluster+0.5,nmodules,-0.5,nmodules-0.5);

  TH2D *hmaxtimebin_xstrip_module = new TH2D("hmaxtimebin_xstrip_module","Max time bin for max X strip on cluster in track", nADCsamples,-0.5,nADCsamples-0.5,nmodules,-0.5,nmodules-0.5);
  TH2D *hmaxtimebin_ystrip_module = new TH2D("hmaxtimebin_ystrip_module","Max time bin for max Y strip on cluster in track", nADCsamples,-0.5,nADCsamples-0.5,nmodules,-0.5,nmodules-0.5);
  //TH2D *hmaxtimebin_xclust_module = new TH2D("hmaxtimebin_xclust_module","Max time bin for cluster-summed X samples on cluster in track", nADCsamples,-0.5,nADCsamples-0.5,nmodules,-0.5,nmodules-0.5);
  //TH2D *hmaxtimebin_yclust_module = new TH2D("hmaxtimebin_yclust_module","Max time bin for cluster-summed Y samples on cluster in track", nADCsamples,-0.5,nADCsamples-0.5);

  TH2D *hADCasym_vs_module = new TH2D("hADCasym_vs_module","(ADCX-ADCY)/(ADCX+ADCY) asymmetry by module, clusters on track",nmodules,-0.5,nmodules-0.5,200,-1.05,1.05);
  TH2D *hdT_vs_module = new TH2D("hdT_vs_module","t_{x} - t_{y} (ns), clusters on track",nmodules,-0.5,nmodules-0.5,200,-50.0,50.0);

  TH2D *hADCsum_Xstrip_max_module = new TH2D("hADCsum_Xstrip_max_module", "ADC sum on max X strip, cluster in track",500,0.0,1.5e4,nmodules,-0.5,nmodules-0.5);
  TH2D *hADCsum_Ystrip_max_module = new TH2D("hADCsum_Ystrip_max_module", "ADC sum on max Y strip, cluster in track",500,0.0,1.5e4,nmodules,-0.5,nmodules-0.5);

  TH2D *hADCsampmax_Xstrip_module = new TH2D("hADCsampmax_Xstrip_module", "Max ADC sample on max X strip, cluster in track",512,0.0,4096.0,nmodules,-0.5,nmodules-0.5);
  TH2D *hADCsampmax_Ystrip_module = new TH2D("hADCsampmax_Ystrip_module", "Max ADC sample on max Y strip, cluster in track",512,0.0,4096.0,nmodules,-0.5,nmodules-0.5);

  TH2D *hADCprodXYstrip_max_module = new TH2D("hADCprodXYstrip_max_module", "#sqrt{ADCX*ADCY} sums for max x,y strips",500,0.0,1.5e4,nmodules,-0.5,nmodules-0.5);
  TH2D *hADCprodXYsamp_max_module = new TH2D("hADCprodXYsamp_max_module", "#sqrt{ADCX*ADCY} max samples, max xy strips",512,0.0,4096.0,nmodules,-0.5,nmodules-0.5);
  TH2D *hADCprodXYclust_module = new TH2D("hADCprodXYclust_module", "#sqrt{ADCX*ADCY} cluster sums",500,0.0,5e4,nmodules,-0.5,nmodules-0.5);
  TH2D *hADCsumXclust_module = new TH2D("hADCsumXclust_module","ADCX cluster sum",500,0.0,5e4,nmodules,-0.5,nmodules-0.5);
  TH2D *hADCsumYclust_module = new TH2D("hADCsumYclust_module","ADCY cluster sum",500,0.0,5e4,nmodules,-0.5,nmodules-0.5);
  TH2D *hStrip_maxcor_module = new TH2D("hStrip_maxcor_module","Correlation coefficient for cluster seed pixel",1000,-1.1,1.1,nmodules,-0.5,nmodules-0.5);
  TH2D *hClust_corr_module = new TH2D("hClust_corr_module","Cluster corr. coeff.",500,-1.1,1.1,nmodules,-0.5,nmodules-0.5);
  
  
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

  TH2D *hClust2D_NstripXvsNstripY = new TH2D("hClust2D_NstripX_vs_NstripY","Nstrips X vs Y",maxnstripXpercluster,0.5,maxnstripXpercluster+0.5,maxnstripYpercluster,0.5,maxnstripYpercluster+0.5);

  TH2D *hClust2D_Xmom_vs_NstripX = new TH2D("hClust2D_Xmom_vs_NstripX","x - xstripmax vs nstripx",maxnstripXpercluster,0.5,maxnstripXpercluster+0.5,600,-3,3);
  TH2D *hClust2D_Ymom_vs_NstripY = new TH2D("hClust2D_Ymom_vs_NstripY","y - ystripmax vs nstripy",maxnstripYpercluster,0.5,maxnstripYpercluster+0.5,600,-3,3);
  
  TH1D *hADCsumXstrip_max = new TH1D("hADCsumXstrip_max","ADC sum for max X strip, cluster on track only",500,0.0,1.5e4);
  TH1D *hADCsumYstrip_max = new TH1D("hADCsumYstrip_max","ADC sum for max Y strip, cluster on track only",500,0.0,1.5e4);

  TH2D *hADCsampXstrip_max = new TH2D("hADCsampXstrip_max","ADC samples for max strip, clusters on track only",6,-0.5,5.5,500,0.0,4000.0);
  TH2D *hADCsampYstrip_max = new TH2D("hADCsampYstrip_max","ADC samples for max strip, clusters on track only",6,-0.5,5.5,500,0.0,4000.0);

  TH2D *hADCsampXclust = new TH2D("hADCsampXclust","Cluster-summed X ADC samples, cluster on track",6,-0.5,5.5,500,0.0,10000.0);
  TH2D *hADCsampYclust = new TH2D("hADCsampYclust","Cluster-summed Y ADC samples, cluster on track",6,-0.5,5.5,500,0.0,10000.0);
  
  TH1D *hADCsampmax_Xstrip = new TH1D("hADCsampmax_Xstrip","max ADC sample for max strip, clusters on track",512,0.0,4096.0);
  TH1D *hADCsampmax_Ystrip = new TH1D("hADCsampmax_Ystrip","max ADC sample for max strip, clusters on track",512,0.0,4096.0);

  TH1D *hADCprodXYstrip_max = new TH1D("hADCprodXYstrip_max","sqrt(ADCx*ADCy) sums for max x,y strips",500,0.0,1.5e4);
  TH1D *hADCprodXYsamp_max = new TH1D("hADCprodXYsamp_max","sqrt(ADCx*ADCy) max samples for max x,y strips",500,0.0,3000.0);
  TH1D *hADCprodXYcluster = new TH1D("hADCprodXYcluster","sqrt(ADCX*ADCY) cluster sums",500.0,0.0,5.0e4);
  
  TH1D *hADCsumXclust_max = new TH1D("hADCsumXclust_max","Cluster ADC sum X, clusters on tracks",500,0.0,5e4);
  TH1D *hADCsumYclust_max = new TH1D("hADCsumYclust_max","Cluster ADC sum Y, clusters on tracks",500,0.0,5e4);
  
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

  // TH1D *hStrip_Afit = new TH1D("hStrip_Afit", "A strip fit", 1000,0.0,300.0);
  // TH1D *hStrip_taufit = new TH1D("hStrip_taufit","tau strip fit", 1000.0,0.0,300.0);
  
  // TH1D *hStrip_Tfit = new TH1D("hStrip_Tfit","t_{strip} fit (ns), strips in clusters on tracks",1000,-150.0,150.0);
  // TH1D *hStrip_TfitX = new TH1D("hStrip_TfitX","t_{strip} fit (ns), X strips in clusters on tracks",1000,-150.0,150.0);
  // TH1D *hStrip_TfitY = new TH1D("hStrip_TfitY","t_{strip} fit (ns), Y strips in clusters on tracks",1000,-150.0,150.0);

  // TH2D *hStrip_TfitXvsY = new TH2D("hStrip_TfitXvsY","Strip tfit_{x} vs. tfit_{y}, strips in cluster on track",300,-150,150,300,-150,150);
  
  TClonesArray *hdidhit_layer = new TClonesArray("TH2D",nlayers);
  TClonesArray *hshouldhit_layer = new TClonesArray("TH2D",nlayers);

  TClonesArray *hxyhit_layer = new TClonesArray("TH2D",nlayers);
  //TClonesArray *hxytrack_layer = new TClonesArray("TH2D",nlayers);

  TClonesArray *hdidhit_module = new TClonesArray("TH2D",nmodules);
  TClonesArray *hshouldhit_module = new TClonesArray("TH2D",nmodules);
  TClonesArray *hxyhit_module = new TClonesArray("TH2D",nmodules);
  
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
    
    TString hnametemp;
    
    new( (*hdidhit_layer)[ilayer] ) TH2D( hnametemp.Format("hdidhit_layer%d",ilayer), "Hit on track in layer", TMath::Nint(nbinsy_eff), ygmin_layer[ilayer]-10.0, ygmax_layer[ilayer]+10.0, TMath::Nint(nbinsx_eff), xgmin_layer[ilayer]-25.0,xgmax_layer[ilayer]+25.0 );
    new( (*hshouldhit_layer)[ilayer] ) TH2D( hnametemp.Format("hshouldhit_layer%d",ilayer), "Track passed through", TMath::Nint(nbinsy_eff), ygmin_layer[ilayer]-10.0, ygmax_layer[ilayer]+10.0, TMath::Nint(nbinsx_eff), xgmin_layer[ilayer]-25.0,xgmax_layer[ilayer]+25.0);
    
    //add a few-mm "buffer zone" at the edges of this histogram:
    new( (*hxyhit_layer)[ilayer] ) TH2D( hnametemp.Format("hxyhit_layer%d",ilayer), "X vs. Y of hits on tracks",
					 nbinsy_hitmap, ygmin_layer[ilayer]-10.0, ygmax_layer[ilayer]+10.0,
					 nbinsx_hitmap, xgmin_layer[ilayer]-25.0, xgmax_layer[ilayer]+25.0 );
    
    
    // new( (*hframe_layers)[ilayer] ) TH2D( hnametemp.Format("hframe_layer%d", ilayer), hnametemp.Format("Layer %d", ilayer),
    // 					  200, ygmin_layer[ilayer]-10.0, ygmax_layer[ilayer]+10.0,
    // 					  200, xgmin_layer[ilayer]-25.0, xgmax_layer[ilayer]+25.0 );
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



  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  double nbins_eff_module = double(ntotal)/400.0*double(nlayers)/double(nmodules);
  double nbinsx_eff_module,nbinsy_eff_module;
  
  for( int imodule=0; imodule<nmodules; imodule++ ){ //Add creation of module-specific efficiency histograms and hit maps:
    TRotation Rtemp;
    Rtemp.RotateX(mod_ax[imodule]);
    Rtemp.RotateY(mod_ay[imodule]);
    Rtemp.RotateZ(mod_az[imodule]);

    mod_Rot[imodule] = Rtemp;
    mod_Rotinv[imodule] = Rtemp.Inverse();

    int nbinsaxis_hitmap = 250;
    int nbinsy_hitmap,nbinsx_hitmap;
    
    if( mod_Lx[imodule] > mod_Ly[imodule] ){
      double ratio = mod_Lx[imodule]/mod_Ly[imodule];
      nbinsy_eff_module = sqrt(nbins_eff_module/ratio);
      nbinsx_eff_module = ratio*nbinsy_eff_module;

      nbinsy_hitmap = nbinsaxis_hitmap;
      nbinsx_hitmap = TMath::Nint(ratio*nbinsaxis_hitmap);
    } else {
      double ratio = mod_Ly[imodule]/mod_Lx[imodule];
      nbinsx_eff_module = sqrt(nbins_eff_module/ratio);
      nbinsy_eff_module = ratio * nbinsx_eff_module;

      nbinsx_hitmap = nbinsaxis_hitmap;
      nbinsy_hitmap = TMath::Nint(ratio*nbinsaxis_hitmap);
    }

    TString hnametemp;
    new( (*hdidhit_module)[imodule] ) TH2D( hnametemp.Format("hdidhit_module%d",imodule),"Hit on track in module",
					    TMath::Nint(nbinsy_eff_module), -mod_Ly[imodule]/2.0-10.0,mod_Ly[imodule]/2.0+10.0,
					    TMath::Nint(nbinsx_eff_module), -mod_Lx[imodule]/2.0-10.0,mod_Lx[imodule]/2.0+10.0 );

    new( (*hshouldhit_module)[imodule] ) TH2D( hnametemp.Format("hshouldhit_module%d",imodule), "Track passed through module",
					       TMath::Nint(nbinsy_eff_module), -mod_Ly[imodule]/2.0-10.0,mod_Ly[imodule]/2.0+10.0,
					       TMath::Nint(nbinsx_eff_module), -mod_Lx[imodule]/2.0-10.0,mod_Lx[imodule]/2.0+10.0 );


    new( (*hxyhit_module)[imodule] ) TH2D( hnametemp.Format("hxyhit_module%d",imodule), "X vs Y of hits on tracks",
					   nbinsy_hitmap, -mod_Ly[imodule]/2.0-10.0,mod_Ly[imodule]/2.0+10.0,
					   nbinsx_hitmap, -mod_Lx[imodule]/2.0-10.0,mod_Lx[imodule]/2.0+10.0 );
    
  }

  NSKIPPED = 0;
  
  while( C->GetEntry(nevent++) && (NMAX < 0 || nevent-firstevent < NMAX ) ){
    /*if( nevent % 1000 == 0 )*/ //cout << nevent << endl;

    EventID = evtID;
    
    if( (nevent-firstevent)%1000 == 0 ){

      cout << "Processing event " << EventID << ", total event count = " << nevent << endl;

    }
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

    for( int imodule=0; imodule<nmodules; imodule++ ){
      moduledata_t modtemp;
      modtemp.Clear();
      ModData.insert( std::pair<int,moduledata_t>(imodule, modtemp) ); 
    }
    
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
	
	hCALfadc_vs_ch->Fill( globalchan, fadc[icalhit] );
	hCALftime_vs_ch->Fill( globalchan, ftimting[icalhit] );
	
	int tmin = 30, tmax=45;
	if( globalchan == 19 ) tmin = 24;
	
	if( globalchan < 21 || globalchan > 26 ){
	  if( maxtimesample >= tmin && maxtimesample <= tmax ){
	    fadcsum += fadc[icalhit];
	  }
	}
      }
      
      hCALfadc_sum->Fill(fadcsum);

      CALOsum = fadcsum;

      int ngoodscint = 0;
      for( int itdc=0; itdc<ntdc; itdc++ ){
	if( tCH[itdc] < 10 && 110.0< ttiming[itdc] && ttiming[itdc] < 135.0 ){
	  ngoodscint++;
	}
	hSCINT_tdc_vs_ch->Fill( tCH[itdc], ttiming[itdc] );
      }

      hNgoodScint->Fill( ngoodscint );

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
      double maxsamp=0.0,sumsamp=0.0,sumgoodsamp=0.0;

      int isamp_max = -1;

      //should we only sum the samples between sampmin and sampmax? not yet. That would probably break too many other parts of the code:
      
      
      for( int isamp=0; isamp<nADCsamples; isamp++ ){
	sumsamp += ADCsamples[isamp];
	//maxsamp = (ADCsamples[isamp] > maxsamp ) ? ADCsamples[isamp] : maxsamp;
	if( isamp_max < 0 || ADCsamples[isamp] > maxsamp ){
	  maxsamp = ADCsamples[isamp];
	  isamp_max = isamp;
	}
	if( isamp >= sampmin_accept && isamp <= sampmax_accept ){
	  sumgoodsamp += ADCsamples[isamp];
	}
      }

      if( maxsamp >= thresh_maxsample && sumsamp >= thresh_stripsum &&
	  isamp_max >= sampmin_accept && isamp_max <= sampmax_accept ) keepstrip = true;
      
      double tsum = 0.0;
      double tsum2 = 0.0;

      double tmean = 0.0;
      double tsigma = 0.0;

      double tcorr;
      
      if( plane == mod_uplaneID[module] && keepstrip ){ //x (vertical/long) axis
	//mod_xstrips_hit[module].insert( strip );

	ModData[module].xstrips.insert( strip );
	ModData[module].ADCsum_xstrips[strip] = 0.0;
	ModData[module].ADCsum_xstrips_goodsamp[strip] = 0.0;
	ModData[module].ADCsamp_xstrips[strip].resize(nADCsamples);
	
	// ADCsum_xstrips[module][strip] = 0.0;
	// ADCsamp_xstrips[module][strip].resize(6);
	double tsamp[nADCsamples];
	double dsamp[nADCsamples];
	double dtsamp[nADCsamples];

	int iAPV = strip/128;
	
	for( int isamp=0; isamp<nADCsamples; isamp++ ){

	  ADCsamples[isamp] *= mod_RYX[module]*mod_Xgain[module][iAPV]; //EXPERIMENTAL: multiply X ADC samples by ratio of Y gain to X gain: everything else proceeds as before:
	  
	  ModData[module].ADCsamp_xstrips[strip][isamp] = ADCsamples[isamp];
	  ModData[module].ADCsum_xstrips[strip] += ADCsamples[isamp];
	  if( isamp >= sampmin_accept && isamp <= sampmax_accept ){
	    ModData[module].ADCsum_xstrips_goodsamp[strip] += ADCsamples[isamp];
	  }
	  // ADCsamp_xstrips[module][strip][isamp] = ADCsamples[isamp];
	  // ADCsum_xstrips[module][strip] += ADCsamples[isamp];

	  hADCvsSampleAllStrips->Fill( isamp, ADCsamples[isamp] );
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

	// TGraphErrors *gtemp = new TGraphErrors(nADCsamples, tsamp, ADCsamples, dtsamp, dsamp);
	// // max of pulse shape occurs at t-t0= tau:
	// // therefore ADCmax = p0*tau*exp(-1) --> p0 = ADCmax/tau*exp(1)
	// PulseShape->SetParameter(0,ModData[module].ADCmax_xstrips[strip]/tau_pulseshape*exp(1.0));
	// double tguess = 12.5+25.0*ModData[module].isampmax_xstrips[strip];
	
	// PulseShape->SetParameter(1,tguess-tau_pulseshape);
	// //PulseShape->FixParameter(2,tau_pulseshape);
	// PulseShape->SetParError(0,20.0);
	// PulseShape->SetParError(1,10.0);
	// PulseShape->SetParError(2,10.0);
	
	// // c2->cd();
	// // gtemp->SetMarkerStyle(20);
	// // gtemp->Draw("AP");
	// gtemp->Fit(PulseShape,"SQ0");

	// gPad->Modified();
	// c2->Update();

	// gSystem->Sleep(10);
	
	// ModData[module].Tfit_xstrips[strip] = PulseShape->GetParameter(1);
	// ModData[module].dTfit_xstrips[strip] = PulseShape->GetParError(1);
	// ModData[module].Afit_xstrips[strip] = PulseShape->GetParameter(0);
	// ModData[module].dAfit_xstrips[strip] = PulseShape->GetParError(0);
	// ModData[module].taufit_xstrips[strip] = PulseShape->GetParameter(2);
	// ModData[module].dtaufit_xstrips[strip] = PulseShape->GetParError(2);
	
	//delete gtemp;
	
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
	ModData[module].ADCsum_ystrips_goodsamp[strip] = 0.0;
	ModData[module].ADCsamp_ystrips[strip].resize(nADCsamples);

	//double tsamp[nADCsamples];
	double tsamp[nADCsamples];
	double dsamp[nADCsamples];
	double dtsamp[nADCsamples];

	int iAPV = strip/128;
	
	for( int isamp=0; isamp<nADCsamples; isamp++ ){

	  ADCsamples[isamp] *= mod_Ygain[module][iAPV];
	  
	  ModData[module].ADCsamp_ystrips[strip][isamp] = ADCsamples[isamp];
	  ModData[module].ADCsum_ystrips[strip] += ADCsamples[isamp];

	  if( isamp >= sampmin_accept && isamp <= sampmax_accept ){
	    ModData[module].ADCsum_ystrips_goodsamp[strip] += ADCsamples[isamp];
	  }
	  
	  hADCvsSampleAllStrips->Fill( isamp, ADCsamples[isamp] );

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

	// double tguess = 12.5+25.0*ModData[module].isampmax_ystrips[strip];
	
	// TGraphErrors *gtemp = new TGraphErrors(nADCsamples, tsamp, ADCsamples, dtsamp, dsamp);
	// // max of pulse shape occurs at t-t0= tau:
	// // therefore ADCmax = p0*tau*exp(-1) --> p0 = ADCmax/tau*exp(1)
	// PulseShape->SetParameter(0,ModData[module].ADCmax_ystrips[strip]/tau_pulseshape*exp(1.0));
	// PulseShape->SetParameter(1,tguess-tau_pulseshape);
	// //	PulseShape->FixParameter(2,tau_pulseshape);
	// PulseShape->SetParError(0,20.0);
	// PulseShape->SetParError(1,10.0);
	// PulseShape->SetParError(2,10.0);

	// // c2->cd();
	// // gtemp->SetMarkerStyle(20);
	// // gtemp->Draw("AP");
	// gtemp->Fit(PulseShape,"SQ0");

	// // gPad->Modified();
	// // c2->Update();
	// // gSystem->Sleep(10);
	
	// ModData[module].Tfit_ystrips[strip] = PulseShape->GetParameter(1);
	// ModData[module].dTfit_ystrips[strip] = PulseShape->GetParError(1);
	// ModData[module].Afit_ystrips[strip] = PulseShape->GetParameter(0);
	// ModData[module].dAfit_ystrips[strip] = PulseShape->GetParError(0);
	// ModData[module].taufit_ystrips[strip] = PulseShape->GetParameter(2);
	// ModData[module].dtaufit_ystrips[strip] = PulseShape->GetParError(2);

	// delete gtemp;
	
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

    if( layers_hitXY.size() >= TOTAL_REQUIRED_HIT ){ //enough layers hit to (possibly) form a track: only bother with clustering and attempted track finding if this is the case: 

      if( eventdisplaymode != 0 ){
	
	for( int imodule=0; imodule<nmodules; imodule++ ){
	  ( (TH2D*) (*hframe_modules)[imodule] )->Reset();

	  ( (TH1D*) (*hmod_ADC_vs_Xstrip)[imodule] )->Reset();
	  ( (TH1D*) (*hmod_ADC_vs_Ystrip)[imodule] )->Reset();

	  ( (TH1D*) (*hmod_ADC_vs_Xstrip_good)[imodule] )->Reset();
	  ( (TH1D*) (*hmod_ADC_vs_Ystrip_good)[imodule] )->Reset();
	  
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

	clusttemp.Clear();

	//mod_clusters.insert( std::pair<int,
	
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

	// hNstripsX_module->Fill( ModData[module].xstrips.size(),  module );
	// hNstripsY_module->Fill( ModData[module].ystrips.size(),  module );
	
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
	
	
	// cout << "starting cluster finding, event " << evtID << ", module " << module << ". (nstrip X, nstrip Y) = (" << ModData[module].xstrips.size()
	//      << ", " << ModData[module].ystrips.size() << ")...";
	//int clusterflag = find_clusters_by_module_new( ModData[module], clusttemp );

	
	
	int clusterflag = find_clusters_by_module_newnew( ModData[module], clusttemp );

	//cout << "   ...done" << endl;
	
	//NEW NEW clustering works, we need to track down and eliminate places where we assume
	// that nclustx == nclusty == nclust2D
	
	if( clusterflag != 0 ){
	  cout << "event " << evtID << ", module " << module << " too noisy, gave up" << endl;
	}

	//hNclust_module->Fill(clusttemp.nclust2D, module );
	
	//cout << "ending cluster finding, ncluster =  " << clusttemp.nclust2D << endl;
	
	//mod_clusters[module] = clusttemp;

	//maybe it's safer to use map::insert here:
	mod_clusters.insert( std::pair<int,clusterdata_t>(module, clusttemp) );
	
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

      }

      int nlayers_with_2Dclust = 0;
      for( int ilay=0; ilay<nlayers; ilay++ ){
	if( Nclust2D_layer[ilay] > 0 ) nlayers_with_2Dclust++;
      }

      hNlayers_2Dclust->Fill( nlayers_with_2Dclust );
      
      trackdata_t tracktemp;

      tracktemp.Clear();
      
      tracktemp.ntracks = 0;

      //      cout << "Finding tracks, event..." << evtID << endl;

      if( nlayers_with_2Dclust >= TOTAL_REQUIRED_HIT ){

	//cout << "Finding tracks, event... " << evtID << ", total event count = " << nevent << endl;
	
	auto start = high_resolution_clock::now();
	//define the forward and backward constraint points
	TVector3 fcp(0, 0, zavg_layer[0] - 10);
	TVector3 bkp(0, 0, zavg_layer[nlayers-1] + 10);

	if( TrackingAlgorithmFlag == 1 ){
	  new_find_tracks( mod_clusters, tracktemp, fcp, bkp);
	} else if (TrackingAlgorithmFlag == 0 ){
	  find_tracks_old( mod_clusters, tracktemp );
	} else {
	  find_tracks(mod_clusters, tracktemp );
	}
	auto end = high_resolution_clock::now();
	totalTime += duration_cast<nanoseconds>(end - start);
	
	//cout << "track finding successful, event..." << evtID << endl;
	
	hNtracks_found->Fill( tracktemp.ntracks );
	
	for( int ilay=0; ilay<nlayers; ilay++ ){
	  Nclustperlayer[ilay] = Nclust2D_layer[ilay];
	}

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
	    
	    hCALfadc_vs_ch_goodtrack->Fill( globalchan, fadc[icalhit] );
	    hCALftime_vs_ch_goodtrack->Fill( globalchan, ftimting[icalhit] );
	    
	    int tmin = 30, tmax=45;
	    if( globalchan == 19 ) tmin = 24;
	    
	    if( globalchan < 21 || globalchan > 26 ){
	      if( maxtimesample >= tmin && maxtimesample <= tmax ){
		fadcsum += fadc[icalhit];
	      }
	    }
	  }

	  hCALfadc_sum_goodtrack->Fill(fadcsum);

	  CALOsum = fadcsum;

	  int ngoodscint = 0;
	  for( int itdc=0; itdc<ntdc; itdc++ ){
	    if( tCH[itdc] < 10 && 110.0< ttiming[itdc] && ttiming[itdc] < 135.0 ){
	      ngoodscint++;
	    }
	    hSCINT_tdc_vs_ch_goodtrack->Fill( tCH[itdc], ttiming[itdc] );
	  }
	  
	  hNgoodScint_goodtrack->Fill( ngoodscint );
	  
	  NGOODSCINT = ngoodscint;
	  
	}
	
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

	    //In order to make the efficiency determination less biased, let's collect a list
	    //of all the layers on the track, and require that at least three OTHER layers
	    //than the layer in question fired.
	    //We will have to do the same when we fill the "should hit" histogram
	    set<int> layersontrack_temp;
	    for( int ihit=0; ihit<tracktemp.nhitsontrack[itrack]; ihit++ ){
	      int modtemp = tracktemp.modlist_track[itrack][ihit];
	      int layertemp = mod_layer[modtemp];

	      layersontrack_temp.insert( layertemp );
	    }
	      
	    for( int ilay=0; ilay<nlayers; ilay++ ){
	      double xtracktemp = TrackX + zavg_layer[ilay]*TrackXp;
	      double ytracktemp = TrackY + zavg_layer[ilay]*TrackYp;

	      // Check whether we have at least three layers OTHER than the layer in question
	      // on this track:

	      int minhits = 3;
	      if( layersontrack_temp.find( ilay ) != layersontrack_temp.end() ){
		//this layer IS on the track; we require at least three OTHER layers to have fired:
		minhits = 4;
		//If the layer is NOT on the track, then we only require three hits in any other layers to fill the "should hit" histogram for that layer for that event.
	      } 

	      if( tracktemp.nhitsontrack[itrack] >= minhits ){
		( (TH2D*) (*hshouldhit_layer)[ilay] )->Fill( ytracktemp, xtracktemp );

		int nearest_module = get_nearest_module(tracktemp, ilay, itrack);

		if( nearest_module >= 0 && nearest_module < nmodules ){
		  ( (TH2D*) (*hshouldhit_module)[nearest_module] )->Fill( ytracktemp - mod_y0[nearest_module],
									  xtracktemp - mod_x0[nearest_module] );

		  hNstripsX_module->Fill( ModData[nearest_module].xstrips.size(), nearest_module );
		  hNstripsY_module->Fill( ModData[nearest_module].ystrips.size(), nearest_module );

		  hNclust_module->Fill( mod_clusters[nearest_module].nclust2D, nearest_module );
		
		}
	      }

	      //Now HERE is where we should fill ADC type histograms:
	      
	      //To fill "should hit" histograms by module, we need to find the closest module in each layer
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

	      hTrackXeresid_vs_layer->Fill( layer, utracktemp - uhittemp ); //resid = utrack - uhit --> uhit = utrack - resid
	      hTrackYeresid_vs_layer->Fill( layer, vtracktemp - vhittemp );

	      hTrackXeresid_vs_module->Fill( module, utracktemp - uhittemp );
	      hTrackYeresid_vs_module->Fill( module, vtracktemp - vhittemp );
	      
	    }

	    int APVx = clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]] / 128 ;
	    int APVy = clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]] / 128 ;

	    ( (TH2D*) (*hmod_ADCmax_vs_APV_Xstrips)[module] )->Fill( APVx, ModData[module].ADCmax_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]] );
	    ( (TH2D*) (*hmod_ADCmax_vs_APV_Ystrips)[module] )->Fill( APVy, ModData[module].ADCmax_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]] );

	    //( (TH2D*) (*hmod_ADCsum_vs_APV_Xstrips)[module] )->Fill( APVx, ModData[module].ADCsum_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]] );
	    //( (TH2D*) (*hmod_ADCsum_vs_APV_Ystrips)[module] )->Fill( APVy, ModData[module].ADCsum_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]] );

	    //These ADC values account for any cluster-splitting which may or may not have occurred:
	    ( (TH2D*) (*hmod_ADCsum_vs_APV_Xstrips)[module] )->Fill( APVx, clusttemp.xstripADCsum[clusttemp.ixclust2D[iclust2D]][clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]-clusttemp.ixstriplo[clusttemp.ixclust2D[iclust2D]]]);
	    ( (TH2D*) (*hmod_ADCsum_vs_APV_Ystrips)[module] )->Fill( APVy, clusttemp.ystripADCsum[clusttemp.iyclust2D[iclust2D]][clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]-clusttemp.iystriplo[clusttemp.iyclust2D[iclust2D]]] );
	    
	    //now fill the "all strips" histograms:

	    int APVmin=APVx, APVmax=APVx;
	    for( int ixstrip=clusttemp.ixstriplo[clusttemp.ixclust2D[iclust2D]]; ixstrip<=clusttemp.ixstriphi[clusttemp.ixclust2D[iclust2D]]; ixstrip++ ){
	      APVx = ixstrip/128;

	      //Check whether this cluster straddles two APVs
	      if( APVx<APVmin ) APVmin = APVx; 
	      if( APVx>APVmax ) APVmax = APVx;
	      
	      double ADCmax_strip = ModData[module].ADCmax_xstrips[ixstrip];
	      double ADCsum_strip = clusttemp.xstripADCsum[clusttemp.ixclust2D[iclust2D]][ixstrip-clusttemp.ixstriplo[clusttemp.ixclust2D[iclust2D]]];
	      
	      ( (TH2D*) (*hmod_ADCmax_vs_APV_Xstrips_all)[module] )->Fill( APVx, ADCmax_strip );
	      ( (TH2D*) (*hmod_ADCsum_vs_APV_Xstrips_all)[module] )->Fill( APVx, ADCsum_strip );
	    }

	    //Fill cluster sums by APV if all strips on the cluster are in the same APV:
	    if( APVmin == APVmax && clusttemp.nstripx2D[iclust2D] >= 2 ){
	      ( (TH2D*) (*hmod_clustersumX_vs_APV)[module] )->Fill( APVmax, clusttemp.totalchargex[clusttemp.ixclust2D[iclust2D]] );
	    }

	    APVmin = APVy;
	    APVmax = APVy;
	    //now fill the "all strips" histograms:
	    for( int iystrip=clusttemp.iystriplo[clusttemp.iyclust2D[iclust2D]]; iystrip<=clusttemp.iystriphi[clusttemp.iyclust2D[iclust2D]]; iystrip++ ){
	      APVy = iystrip/128;

	      //Check whether this cluster straddles two APVs
	      if( APVy<APVmin ) APVmin = APVy; 
	      if( APVy>APVmax ) APVmax = APVy;
	      
	      double ADCmax_strip = ModData[module].ADCmax_ystrips[iystrip];
	      double ADCsum_strip = clusttemp.ystripADCsum[clusttemp.iyclust2D[iclust2D]][iystrip-clusttemp.iystriplo[clusttemp.iyclust2D[iclust2D]]];
	      
	      ( (TH2D*) (*hmod_ADCmax_vs_APV_Ystrips_all)[module] )->Fill( APVy, ADCmax_strip );
	      ( (TH2D*) (*hmod_ADCsum_vs_APV_Ystrips_all)[module] )->Fill( APVy, ADCsum_strip );
	    }

	    if( APVmin == APVmax && clusttemp.nstripy2D[iclust2D] >= 2 ){
	      ( (TH2D*) (*hmod_clustersumY_vs_APV)[module] )->Fill( APVmax, clusttemp.totalchargey[clusttemp.iyclust2D[iclust2D]]);
	    }

	    
	    hclustwidth_x_module->Fill( clusttemp.nstripx2D[iclust2D], module );
	    hclustwidth_y_module->Fill( clusttemp.nstripy2D[iclust2D], module );

	    hmaxtimebin_xstrip_module->Fill( ModData[module].isampmax_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]], module );
	    hmaxtimebin_ystrip_module->Fill( ModData[module].isampmax_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]], module );

	    hADCasym_vs_module->Fill( module, clusttemp.dEclust2D[iclust2D]/(2.0*clusttemp.Eclust2D[iclust2D]) );
	    hdT_vs_module->Fill( module, clusttemp.dtclust2D[iclust2D] );

	    hADCsum_Xstrip_max_module->Fill( ModData[module].ADCsum_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]], module );
	    hADCsum_Ystrip_max_module->Fill( ModData[module].ADCsum_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]], module );

	    hADCsampmax_Xstrip_module->Fill( ModData[module].ADCmax_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]], module );
	    hADCsampmax_Ystrip_module->Fill( ModData[module].ADCmax_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]], module );

	    if( clusttemp.nstripx2D[iclust2D] >= 2 )
	      hADCsumXclust_module->Fill( clusttemp.totalchargex[clusttemp.ixclust2D[iclust2D]], module );
	    if( clusttemp.nstripy2D[iclust2D] >= 2 )
	      hADCsumYclust_module->Fill( clusttemp.totalchargey[clusttemp.iyclust2D[iclust2D]], module );

	    hADCprodXYstrip_max_module->Fill( sqrt(ModData[module].ADCsum_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]]*
						   ModData[module].ADCsum_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]] ), module );
	    hADCprodXYsamp_max_module->Fill( sqrt(ModData[module].ADCmax_xstrips[clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]]]*
						  ModData[module].ADCmax_ystrips[clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]]] ), module );
	    hADCprodXYclust_module->Fill( sqrt( clusttemp.totalchargex[clusttemp.ixclust2D[iclust2D]]*clusttemp.totalchargey[clusttemp.iyclust2D[iclust2D]] ), module );

	    hStrip_maxcor_module->Fill( clusttemp.CorrCoeffMaxStrips[iclust2D], module );
	    hClust_corr_module->Fill( clusttemp.CorrCoeff2D[iclust2D], module );
	    
	  
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

	      hADCsampXclust->Fill( isamp, clusttemp.ADCsamp_xclust[clusttemp.ixclust2D[iclust2D]][isamp] );
	      hADCsampYclust->Fill( isamp, clusttemp.ADCsamp_yclust[clusttemp.iyclust2D[iclust2D]][isamp] );
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

		// hStrip_Tfit->Fill( ModData[module].Tfit_xstrips[ixstrip] );
		// hStrip_TfitX->Fill( ModData[module].Tfit_xstrips[ixstrip] );

		// hStrip_Afit->Fill( ModData[module].Afit_xstrips[ixstrip] );
		// hStrip_taufit->Fill( ModData[module].taufit_xstrips[ixstrip]);
		
	      }
	    }

	    //loop over strips in cluster and fill strip deltaT histograms:
	    for( int iystrip=clusttemp.iystriplo[clusttemp.iyclust2D[iclust2D]]; iystrip<=clusttemp.iystriphi[clusttemp.iyclust2D[iclust2D]]; iystrip++ ){
	      if( ModData[module].Tmean_ystrips.find( iystrip ) != ModData[module].Tmean_ystrips.end() ){
		hStrip_dT->Fill( ModData[module].Tmean_ystrips[iystrip] - clusttemp.tclust2D[iclust2D] );
		hStrip_dT_vs_ADCsum->Fill( ModData[module].ADCsum_ystrips[iystrip], ModData[module].Tmean_ystrips[iystrip] - clusttemp.tclust2D[iclust2D] );

		hStrip_dTwalkcor->Fill( ModData[module].Tmean_ystrips_walkcor[iystrip] - clusttemp.tclust2Dwalkcorr[iclust2D] );
		hStrip_dTwalkcor_vs_ADC->Fill( ModData[module].ADCsum_ystrips[iystrip], ModData[module].Tmean_ystrips_walkcor[iystrip] - clusttemp.tclust2Dwalkcorr[iclust2D] );

		// hStrip_Tfit->Fill( ModData[module].Tfit_ystrips[iystrip] );
		// hStrip_TfitY->Fill( ModData[module].Tfit_ystrips[iystrip] );

		// hStrip_Afit->Fill( ModData[module].Afit_ystrips[iystrip] );
		// hStrip_taufit->Fill( ModData[module].taufit_ystrips[iystrip]);
		
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
	      // IF we are here, then it is because a hit DID occur in this layer/module
	      // on this track;
	      // This means that we have to require that the total number of hits on the
	      // track be at least equal to four to fill the "did hit" histogram
	      
	      if( tracktemp.nhitsontrack[itrack] >= 4 ){
		( (TH2D*) (*hdidhit_layer)[layer] )->Fill( TrackY + zavg_layer[layer]*TrackYp,
							   TrackX + zavg_layer[layer]*TrackXp );
	      }
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

	      HitXstripMax[ihit] = clusttemp.ixstripmax[clusttemp.ixclust2D[iclust2D]];
	      HitYstripMax[ihit] = clusttemp.iystripmax[clusttemp.iyclust2D[iclust2D]];

	      HitXstripLo[ihit] = clusttemp.ixstriplo[clusttemp.ixclust2D[iclust2D]];
	      HitYstripLo[ihit] = clusttemp.iystriplo[clusttemp.iyclust2D[iclust2D]];

	      HitXstripHi[ihit] = clusttemp.ixstriphi[clusttemp.ixclust2D[iclust2D]];
	      HitYstripHi[ihit] = clusttemp.iystriphi[clusttemp.iyclust2D[iclust2D]];
	      
	      hClust2D_Xmom_vs_NstripX->Fill( HitNstripX[ihit], HitXmom[ihit] );
	      hClust2D_Ymom_vs_NstripY->Fill( HitNstripY[ihit], HitYmom[ihit] );

	      ( (TH2D*) (*hxyhit_module)[module] )->Fill( HitYlocal[ihit], HitXlocal[ihit] );
	      
	      //To get a proper efficiency numerator, we need this check:

	      //if( get_nearest_module( tracktemp, layer, itrack ) == module ){
	      //this check shouldn't be needed any more, since get_nearest_module(track,layer) always defaults to the module with a hit
	      // if the hit actually occurred.
	      if( tracktemp.nhitsontrack[itrack] >= 4 ){
		( (TH2D*) (*hdidhit_module)[module] )->Fill( TrackY + zavg_layer[layer]*TrackYp - mod_y0[module],
							     TrackX + zavg_layer[layer]*TrackXp - mod_x0[module] );
	      }
	      //}
	    }
	  }
      

	  hTrackXY->Fill( tracktemp.Ytrack[itrack],tracktemp.Xtrack[itrack] );
	  hTrackXp->Fill( tracktemp.Xptrack[itrack] );
	  hTrackYp->Fill( tracktemp.Yptrack[itrack] );
	}
      }

      hNtracks_nocuts->Fill( tracktemp.ntracks );
      if( CALOsum > 5500.0 && NGOODSCINT >= 2 ){
	hNtracks_GoodScintAndCALO->Fill( tracktemp.ntracks );
      }
      
      if( tracktemp.ntracks > 0 ){
	Tout->Fill();
      }
      
      if( eventdisplaymode != 0 ){

	cout << "starting event display, ievent = " << nevent << endl;
	
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

	    //Fill strip ADC histogram. We'll need a canvas for this later
	    ( (TH1D *) (*hmod_ADC_vs_Xstrip)[module] )->Fill( strip, ModData[module].ADCsum_xstrips_goodsamp[strip] );
	    
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

	    //Fill strip ADC histogram. We'll need a canvas for this later
	    ( (TH1D *) (*hmod_ADC_vs_Ystrip)[module] )->Fill( strip, ModData[module].ADCsum_ystrips_goodsamp[strip] );
	    
	    
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

	  //Now let's draw the strip ADC histograms:
	  c2x->cd(module_ipad[module]);
	  ( (TH1D*) (*hmod_ADC_vs_Xstrip)[module] )->Draw("hist");
	  c2y->cd(module_ipad[module]);
	  ( (TH1D*) (*hmod_ADC_vs_Ystrip)[module] )->Draw("hist");
	  
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

	  int mstyle_stripADC[] = {2,3,4,5,20,21,22,23,25,25};
	  int mcol_stripADC[] = {1,2,4,6,8,9,30,28,40};
	  
	  for( int iclust=0; iclust<clusttemp.nclust2D; iclust++ ){
	    
	    if( clusttemp.keepclust2D[iclust] ){
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

	    //Now also draw the strip signals for this hit: We may eventually wish to only draw clusters that passed 2D filtering criteria and/or ended up on good
	    //tracks, but we would need to modify the cluster data structure to make that work:
	    for( int ixclust = 0; ixclust<clusttemp.nclustx; ixclust++ ){
	      TMarker MstripADC;
	      MstripADC.SetMarkerStyle( mstyle_stripADC[ixclust%10] );
	      MstripADC.SetMarkerColor( mcol_stripADC[ixclust%9] );
	      
	      int stripcount=0;
	      for( int istripx=clusttemp.ixstriplo[ixclust]; istripx<=clusttemp.ixstriphi[ixclust]; istripx++ ){
		double xADC = clusttemp.xstripADCsum[ixclust][stripcount];
		
		c2x->cd( module_ipad[module] );
		MstripADC.DrawMarker( istripx,xADC );
		stripcount++;
		//( (*TH1D) (*hmod_ADC_vs_Xstrip_good)[module] )->Fill( istripx, xADC );
	      }
	    }

	    for( int iyclust = 0; iyclust<clusttemp.nclusty; iyclust++ ){
	      TMarker MstripADC;
	      MstripADC.SetMarkerStyle( mstyle_stripADC[iyclust%10] );
	      MstripADC.SetMarkerColor( mcol_stripADC[iyclust%9] );

	      int stripcount = 0;
	      for( int istripy=clusttemp.iystriplo[iyclust]; istripy<=clusttemp.iystriphi[iyclust]; istripy++ ){
		double yADC = clusttemp.ystripADCsum[iyclust][stripcount];

		c2y->cd( module_ipad[module] );
		MstripADC.DrawMarker( istripy,yADC );
		stripcount++;
	      }

	      gPad->Modified();
	      gPad->Update();
	      //gSystem->ProcessEvents();
	      c1->Update();
	      
	      //( (*TH1D) (*hmod_ADC_vs_Xstrip_good)[module] )->Fill( istripx, xADC );
	    }
	    
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
      
	cout << "press any key to continue (q to quit, p to print):" << endl;
	TString reply;
	reply.ReadLine(cin,kFALSE);

	if( reply.BeginsWith("q") ) break;
	if( reply.BeginsWith("p") ){
	  TString pdffilename;
	  c1->Print( pdffilename.Format( "plots/c1_event_%d.pdf", nevent ), "pdf");
	  c2x->Print( pdffilename.Format( "plots/c2x_event_%d.pdf", nevent ), "pdf");
	  c2y->Print( pdffilename.Format( "plots/c2y_event_%d.pdf", nevent ), "pdf");
	  c_proj->Print( pdffilename.Format( "plots/cproj_event_%d.pdf", nevent ), "pdf");
	}
      }
	  
	  
      //gSystem->Sleep(2500);

      //	  c1->Update();
	  
	  
    } else {
      hNtracks_nocuts->Fill( 0 );
      if( CALOsum >= 5500.0 && NGOODSCINT >= 2 ){
	hNtracks_GoodScintAndCALO->Fill( 0 );
      }
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

  TClonesArray *heff_module = new TClonesArray("TH2D",nmodules);

  for( int imod=0; imod<nmodules; imod++ ){
    TString histname;

    TH2D *htemp = ( (TH2D*) (*hdidhit_module)[imod] );
    new( (*heff_module)[imod] ) TH2D( *htemp );
    
    ( (TH2D*) (*heff_module)[imod] )->SetName( histname.Format( "heff_module%d", imod ) );
    ( (TH2D*) (*heff_module)[imod] )->SetTitle( histname.Format( "Track-based efficiency, module %d", imod ) );
    
    ( (TH2D*) (*heff_module)[imod] )->Divide( (TH2D*) (*hshouldhit_module)[imod] ); //
    //Fix ROOT's 
    
  }
  
  fout->Write();
  
  auto program_end = high_resolution_clock::now();
  
  cout<<"total time spent in tracking: "<<(double)totalTime.count() / 1e9<<" seconds "<<endl;
  
  auto program_time = duration_cast<nanoseconds>(program_end - program_start); 
  
  cout<<"total program execution time: "<<(double)program_time.count() / 1e9<<" seconds" <<endl;
  
  cout<<"fraction of time spent in track-finding: "<<((double)totalTime.count()/(double)program_time.count())<<endl;

  cout << "Number of events with tracking skipped due to too many hit combos = " << NSKIPPED << ", Total events processed = "
       << nevent - firstevent << endl;
  cout << "Analysis rate = " << double(nevent-firstevent)/((double) program_time.count()/1.e9) << endl;
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
