#include "GEM_cosmic_tracks.C"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TEventList.h"
#include "TCut.h"
#include <iostream>
#include <fstream>
#include "TMinuit.h"
#include "TMath.h"

double PI = TMath::Pi();

int nlayers=4;
int nmodules=12;

//Make these global for chi^2 function for numerical minimization:

map<int,int> mod_layer;
map<int,double> mod_x0; //key = module index, value = x center position
map<int,double> mod_y0; //key = module index, value = y center position
map<int,double> mod_z0; //key = module index, value = z position
map<int,double> mod_ax; //key = module index, value = X axis rotation (yaw)
map<int,double> mod_ay; //key = module index, value = Y axis rotation (pitch)
map<int,double> mod_az; //key = module index, value = Z axis rotation (roll)
map<int,double> mod_dx0;
map<int,double> mod_dy0;
map<int,double> mod_dz0;
map<int,double> mod_dax;
map<int,double> mod_day;
map<int,double> mod_daz;
//Need "U" and "V" strip angle definitions to generalize to the case of arbitrary strip orientation:
map<int,double> mod_uangle; //key = module index, value = angle between module x axis and "U" direction
map<int,double> mod_vangle; //key = module index, value = angle between modyle x axis and "V" direction
map<int,double> mod_Pxu;    //cos(uangle);
map<int,double> mod_Pyu;    //sin(uangle);
map<int,double> mod_Pxv;    //cos(vangle);
map<int,double> mod_Pyv;    //cos(vangle);

map<int,bool> fixmod; //allowing fixing the position and orientation of arbitrary combinations of modules:

long NMAX;

//Let's see if we can improve things by doing one iteration of linearized, then use TMinuit:
int NTRACKS;
vector<double> XTRACK,YTRACK,XPTRACK,YPTRACK;
vector<int> TRACKNHITS;
vector<vector<int> > HITMOD;
vector<vector<double> > HITX,HITY;

void CHI2_FCN( int &npar, double *gin, double &f, double *par, int flag){

  double chi2=0.0;

  //for( int pass=0; pass<2; pass++ ){ //First pass, re-fit tracks, Second pass: fit parameters:
  for( int itr=0; itr<NTRACKS; itr++ ){
    //First: re-fit track using latest parameters:
    
    for( int ihit=0; ihit<TRACKNHITS[itr]; ihit++ ){
      double ulocal = HITX[itr][ihit];
      double vlocal = HITY[itr][ihit];

      int module = HITMOD[itr][ihit];

      int ipar_x0 = 6*module;
      int ipar_y0 = 6*module+1;
      int ipar_z0 = 6*module+2;
      int ipar_ax = 6*module+3;
      int ipar_ay = 6*module+4;
      int ipar_az = 6*module+5;
      
      double det = mod_Pxu[module]*mod_Pyv[module] - mod_Pyu[module]*mod_Pxv[module];

      double xlocal = (mod_Pyv[module]*ulocal - mod_Pyu[module]*vlocal)/det; //(sin(alphav)*U - sin(alphau)*V)/det = U = X for alphau = 0, alphav = 90
      double ylocal = (mod_Pxu[module]*vlocal - mod_Pxv[module]*ulocal)/det; //(cos(alphau)*V - cos(alphav)*U)/det = V = Y for alphau = 0, alphav = 90

      TVector3 hitpos_local(xlocal,ylocal,0);
      TRotation R;
      R.RotateX( mod_ax[module] + par[ipar_ax] );
      R.RotateY( mod_ay[module] + par[ipar_ay] );
      R.RotateZ( mod_az[module] + par[ipar_az] );

      TVector3 modcenter_global( mod_x0[module]+par[ipar_x0],
				 mod_y0[module]+par[ipar_y0],
				 mod_z0[module]+par[ipar_z0] );

      TVector3 hitpos_global = modcenter_global + R*hitpos_local;

      TVector3 trackpos_global( XTRACK[itr]+XPTRACK[itr]*hitpos_global.Z(),
				YTRACK[itr]+YPTRACK[itr]*hitpos_global.Z(),
				hitpos_global.Z() );
	
	
	
      TVector3 diff = hitpos_global - trackpos_global;
	
      chi2 += pow(diff.X(),2)+pow(diff.Y(),2); //
	
      //update tracks? this is the tricky part:
    }

    // if( pass == 0 ){ //update tracks:
    // 	TVectorD Track = Atrack.Invert() * btrack;

    // 	XTRACK[itr] = Track(0);
    // 	XPTRACK[itr] = Track(1);
    // 	YTRACK[itr] = Track(2);
    // 	YPTRACK[itr] = Track(3);
    // } //loop over hits
  } //loop over tracks
    //loop over two passes

  f = chi2;
  
}



void GEM_align( const char *inputfilename, const char *configfilename, const char *outputfilename="newGEMalignment.txt" ){
  ifstream configfile(configfilename);
  
  int niter=1; //number of alignment iterations:

  //default to 100000
  NMAX=100000; //limit number of events so we can maximize number of alignment iterations.
  
  int refmod=4;
  
  TCut globalcut = "";

  int offsetsonlyflag = 0;
  int rotationsonlyflag = 0;

  int fixz = 0; //fix z coordinate of layer if our data don't give us enough sensitivity to determine the z coordinates:
  int fixax=0, fixay=0, fixaz=0;

  double sigma_hitpos=0.14; //mm
  
  //Copied from GEM_reconstruct: For this routine we are only interested in the number of layers and the number of modules, and the geometrical information:
  if( configfile ){
    TString currentline;
    
    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig")){
      if( !currentline.BeginsWith("#") ){
	TObjArray *tokens = currentline.Tokenize(" ");

	int ntokens = tokens->GetEntries();

	if( ntokens >= 2 ){
	  TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	  if( skey == "niter" ){
	    TString sniter =  ( (TObjString*) (*tokens)[1] )->GetString();
	    niter = sniter.Atoi();
	  }
	  
	  if( skey == "nlayers" ){
	    TString snlayers = ( (TObjString*) (*tokens)[1] )->GetString();
	    nlayers = snlayers.Atoi();
	  }

	  if( skey == "offsetsonly" ){
	    TString sflag = ( (TObjString*) (*tokens)[1] )->GetString();
	    offsetsonlyflag = sflag.Atoi();
	  }

	  if( skey == "rotationsonly" ){
	    TString sflag = ( (TObjString*) (*tokens)[1] )->GetString();
	    rotationsonlyflag = sflag.Atoi();
	  }
	  
	  if( skey == "nmodules" ){
	    TString snmodules = ( (TObjString*) (*tokens)[1] )->GetString();
	    nmodules = snmodules.Atoi();
	  }

	  if( skey == "refmod" ){
	    TString sflag = ( (TObjString*) (*tokens)[1] )->GetString();
	    refmod = sflag.Atoi();
	  }

	  if( skey == "fixmod" && ntokens >= nmodules + 1 ){
	    for( int i=1; i<ntokens; i++ ){
	      TString stemp = ( (TObjString*) (*tokens)[i] )->GetString();
	      int flagtemp = stemp.Atoi();
	      fixmod[i-1] = (flagtemp != 0 );
	    }
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

	  if( skey == "fixz" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    fixz = stemp.Atoi();

	    cout << "setting fixz = " << fixz << endl;
	  }
	      
	  if( skey == "fixax" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    fixax = stemp.Atoi();
	  }

	  if( skey == "fixay" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    fixay = stemp.Atoi();
	  }

	  if( skey == "fixaz" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    fixaz = stemp.Atoi();
	  }

	  if( skey == "sigma" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    sigma_hitpos = stemp.Atof();
	  }

	  if( skey == "NMAX" && ntokens >= 2 ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    NMAX = stemp.Atoi();
	  }
	      
	  // if( skey == "eventdisplay" && ntokens >= 2 ){
	  //   TString sevdisplay = ( (TObjString*) (*tokens)[1] )->GetString();

	  //   eventdisplaymode = sevdisplay.Atoi();
	  // }

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

	  // if( skey == "maxstripspercluster" && ntokens >= 2 ){ //largest number of hits to be considered as part of same cluster:
	  //   TString smaxstrips =  ( (TObjString*) (*tokens)[1] )->GetString();
	  //   maxnstripspercluster = smaxstrips.Atoi();
	  // }

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

	 
	  // if( skey == "maxcor_threshold" && ntokens >= 2 ){ //threshold on MAX correlation coefficient:
	  //   TString sthreshold = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   maxstripcorthreshold = sthreshold.Atof();
	  // }
	  
	  // if( skey == "stripcor_threshold" && ntokens >= 2 ){
	  //   TString sthreshold = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   stripcorthreshold = sthreshold.Atof();
	  // }

	  // if( skey == "clustcor_threshold" && ntokens >= 2 ){
	  //   TString sthreshold = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   clustcorthreshold = sthreshold.Atof();
	  // }

	  // if( skey == "clust2D_ADCasymcut" && ntokens >= 2 ){
	  //   TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   cluster2Dmatch_asymcut = scut.Atof();
	  // }

	  // if( skey == "clust2D_dTcut" && ntokens >= 2 ){
	  //   TString scut = ( (TObjString*) (*tokens)[1] )->GetString();

	  //   cluster2Dmatch_tcut = scut.Atof();
	  // }

	  // if( skey == "threshold_maxsample" && ntokens >= 2 ){
	  //   TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   thresh_maxsample = scut.Atof();
	  // }

	  // if( skey == "threshold_stripsum" && ntokens >= 2 ){
	  //   TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   thresh_stripsum = scut.Atof();
	  // }

	  // if( skey == "threshold_clustersum" && ntokens >= 2 ){
	  //   TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   thresh_clustersum = scut.Atof();
	  // }

	  // if( skey == "trackchi2cut" && ntokens >= 2 ){
	  //   TString scut = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   TrackChi2Cut = scut.Atof();
	  // }
	  
	}
      }
    }
    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut") ){
      if( !currentline.BeginsWith("#") ){
	globalcut += currentline;
      }
    }
  } else {
    return;
  }

  if( offsetsonlyflag != 0 && rotationsonlyflag != 0 ){
    cout << "nothing to align, quitting..." << endl;
    return;   
  }

  TChain *C = new TChain("Tout");
  C->Add(inputfilename);

  TEventList *elist = new TEventList("elist");

  C->Draw(">>elist",globalcut);
  
  GEM_cosmic_tracks *T = new GEM_cosmic_tracks(C);

  long nevent=0;

  
  cout << "Alignment starting solution: " << endl;
  for( map<int,double>::iterator imod=mod_x0.begin(); imod!=mod_x0.end(); ++imod ){
    int module = imod->first;
    cout << "Module " << module << ": (x0,y0,z0,ax,ay,az)=("
	 << mod_x0[module] << ", " << mod_y0[module] << ", " << mod_z0[module] << ", "
	 << mod_ax[module] << ", " << mod_ay[module] << ", " << mod_az[module] << ")" << endl;
  }

 
  //niter = 1;
  
  for( int iter=0; iter<=niter; iter++ ){
    
    nevent=0;
    
    //For each alignment iteration, let's define our chi^2 function of the global geometrical alignment parameters
    //and then try to linearize the problem:
    //    int nparam = nmodules*6; //order of parameters is x0,y0,z0,ax,ay,az
    int nparam = 0;

    int nfreemodules=0;
    vector<int> freemodlist;
    map<int,int> freemodindex;
    
    for( int imod=0; imod<nmodules; imod++ ){
      if( !fixmod[imod] ){
	nparam += 6;
	freemodlist.push_back( imod );
	freemodindex[imod] = nfreemodules;
	nfreemodules++;
      }
    }
    
    if( refmod >= 0 && refmod < nmodules ){
      nparam = (nmodules-1)*6;
    }

    if( nparam == 0 ){
      cout << "all modules fixed, nothing to align... quitting" << endl;
      break;
    }

    cout << "nparam = " << nparam << endl;
    
    TMatrixD M(nparam,nparam);
    TVectorD b(nparam);
    
    for( int ipar=0; ipar<nparam; ipar++ ){
      for( int jpar=0; jpar<nparam; jpar++ ){
	M(ipar,jpar) = 0.0;
      }
      b(ipar) = 0.0;
    }

    //To simplify the special cases for alignment only or rotation only, we define separate matrices for alignment only or position only fits:
    int nparam_rot = nmodules*3;
    int nparam_pos = nmodules*3;

    if( refmod >= 0 && refmod < nmodules ){
      nparam_rot = (nmodules-1)*3;
      nparam_pos = (nmodules-1)*3;
    }
    
    TMatrixD Mrot(nparam_rot,nparam_rot);
    TVectorD brot(nparam_rot);

    TMatrixD Mpos(nparam_pos,nparam_pos);
    TVectorD bpos(nparam_pos);

    for( int ipar=0; ipar<nparam_rot; ipar++ ){
      for(int jpar=0; jpar<nparam_rot; jpar++ ){
	Mrot(ipar,jpar) = 0.0;
	Mpos(ipar,jpar) = 0.0;
      }
      brot(ipar) = 0.0;
      bpos(ipar) = 0.0;
    }
      
    
    // We wish to minimize the sum of squared residuals between all hits and tracks by varying the x,y,z position offsets and ax,ay,az
    // rotation angles of all nmodules modules:
    //For example:
    // chi^2 = sum_{i=1}^Nevent sum_{j=1}^{Nhit} (xhit_ij - xtrack_ij)^2/sigxij^2 + (yhit_ij - ytrack_ij)^2/sigyij^2
    // poslocal = (xlocal,ylocal,0);
    // posglobal = R*poslocal + modcenter
    // R = Rz*Ry*Rx
    // Rx mixes y and z components:
    // Rx = | 1       0         0        |  x' = x    
    //      | 0       cos(ax)   -sin(ax) |  y' = cos(ax)*y -sin(ax)*z ~= y
    //      | 0       sin(ax)   cos(ax)  |  z' = sin(ax)*y +cos(ax)*z ~= +ax*y
    // Similarly for Ry, we have:
    // x' = cos(ay)*x + sin(ay)*z ~= x
    // y' = y
    // z' = -sin(ay)*x + cos(ay)*z ~= -ay*x
    // Similarly for Rz, we have:
    // x' = cos(az)*x - sin(az)*y ~= x - az*y
    // y' = sin(az)*x + cos(az)*y ~= az*x + y
    // Linearized global transformation:
    // xglobal = x - modaz*y + modx0
    // yglobal = y + modaz*x + mody0
    // zglobal = modax*y - moday*x + modz0
    // chi^2 = sum_i,j (xlocal - modaz*ylocal + modx0 - (xtrack + xptrack*(modz0 + modax*ylocal - moday*xlocal)))^2/sigx^2 +
    //                 (ylocal + modaz*xlocal + mody0 - (ytrack + yptrack*(modz0 + modax*ylocal - moday*xlocal)))^2/sigy^2
    // dchi2/dx0 = 2*(xlocal - modaz*ylocal + modx0 - (xtrack + xptrack*(modz0 + modax*ylocal - moday*xlocal))/sigx^2*1
    // dchi2/dy0 = 2*(ylocal + modaz*xlocal + mody0 - (ytrack + yptrack*(modz0 + modax*ylocal - moday*xlocal))/sigy^2*1
    // dchi2/dz0 = 2*(xlocal - modaz*ylocal + modx0 - (xtrack + xptrack*(modz0 + modax*ylocal - moday*xlocal))/sigx^2*-xptrack +
    //             2*(ylocal + modaz*xlocal + mody0 - (ytrack + yptrack*(modz0 + modax*ylocal - moday*xlocal))/sigy^2*-yptrack
    //dchi2/dax =

    // if( iter < niter ){
    //   offsetsonlyflag = true;
    // } else {
    //   offsetsonlyflag = false;
    //   rotationsonlyflag = false;
    // }

    cout << "Starting linearized alignment procedure, iteration = " << iter << endl;
    
    if( iter == niter ){ //only fill these arrays on the last iteration:
      NTRACKS = 0;
    }
    while( T->GetEntry(elist->GetEntry(nevent++))){
      
      if( nevent % 100000 == 0 ){
	cout << "Linearized alignment, nevent = " << nevent << endl;
      }
      
      double xptrack,yptrack,xtrack,ytrack;
      if( iter < 0 ){ //on first iteration use track from ROOT tree:
	xptrack = T->TrackXp;
	yptrack = T->TrackYp;
	xtrack = T->TrackX;
	ytrack = T->TrackY;
      } else { //on ALL iterations, we re-fit the track using updated alignment parameters (or the initial ones from the config file)
	double sumX=0.0, sumY = 0.0, sumZ = 0.0, sumXZ = 0.0, sumYZ = 0.0, sumZ2 = 0.0;

	for( int ihit=0; ihit<T->TrackNhits; ihit++ ){
	  int module = T->HitModule[ihit];

	  double ulocal = T->HitXlocal[ihit]; //"U" local: generalized "X"
	  double vlocal = T->HitYlocal[ihit]; //"V" local: generalized "Y"
	  
	  double det = mod_Pxu[module]*mod_Pyv[module] - mod_Pyu[module]*mod_Pxv[module]; //cos( alphau) * sin(alphav) - sin(alphau)*cos(alphav) = 1 for alphau = 0, alphav = 90
	  
	  double xlocal = (mod_Pyv[module]*ulocal - mod_Pyu[module]*vlocal)/det; //(sin(alphav)*U - sin(alphau)*V)/det = U = X for alphau = 0, alphav = 90
	  double ylocal = (mod_Pxu[module]*vlocal - mod_Pxv[module]*ulocal)/det; //(cos(alphau)*V - cos(alphav)*U)/det = V = Y for alphau = 0, alphav = 90

	  TVector3 hitpos_local(xlocal,ylocal,0);
	  TRotation R;
	  R.RotateX( mod_ax[module] );
	  R.RotateY( mod_ay[module] );
	  R.RotateZ( mod_az[module] );
	
	  TVector3 modcenter_global( mod_x0[module],mod_y0[module],mod_z0[module] );
	  TVector3 hitpos_global = modcenter_global + R*hitpos_local;
	
	  double sigma = 0.1;
	  double weight = pow(sigma_hitpos,-2);

	  weight = 1.0;
	  
	  sumX += hitpos_global.X()*weight;
	  sumY += hitpos_global.Y()*weight;
	  sumZ += hitpos_global.Z()*weight;
	  sumXZ += hitpos_global.X()*hitpos_global.Z()*weight;
	  sumYZ += hitpos_global.Y()*hitpos_global.Z()*weight;
	  sumZ2 += pow(hitpos_global.Z(),2)*weight;	  

	}

	double nhits = double(T->TrackNhits);
	
	double denom = (sumZ2*nhits - pow(sumZ,2));
	xptrack = (nhits*sumXZ - sumX*sumZ)/denom;
	yptrack = (nhits*sumYZ - sumY*sumZ)/denom;
	xtrack = (sumZ2*sumX - sumZ*sumXZ)/denom;
	ytrack = (sumZ2*sumY - sumZ*sumYZ)/denom;

	// cout << "Old track (xp,yp,x,y)=(" << T->TrackXp << ", " << T->TrackYp << ", " << T->TrackX << ", " << T->TrackY
	//      << ")" << endl;
	// cout << "New track (xp,yp,x,y)=(" << xptrack << ", " << yptrack << ", " << xtrack << ", " << ytrack << ")" << endl;
	
      }
	
      if( nevent < NMAX && iter == niter ){ //fill TRACK arrays 
	NTRACKS++;
	XTRACK.push_back( xtrack );
	XPTRACK.push_back( xptrack );
	YTRACK.push_back( ytrack );
	YPTRACK.push_back( yptrack );
	TRACKNHITS.push_back( T->TrackNhits );
      }
    
      //we want to modify this to compute the CHANGE in module parameters required to minimize chi^2;
      // so the starting parameters are taken as given.
      // x_0 --> x_0 + dx0
      // y_0 --> y_0 + dy0
      // z_0 --> z_0 + dz0
      // ax --> ax + dax
      // ay --> ay + day
      // az --> az + daz
      // The coefficients of the changes in the parameters should stay the same as those of the parameters themselves, but the RHS needs modified:
      vector<int> HITMODTEMP;
      vector<double> HITXTEMP,HITYTEMP;
      
      for( int ihit=0; ihit<T->TrackNhits; ihit++ ){
	int module = T->HitModule[ihit];
	double sigx = T->HitSigX[ihit];
	double sigy = T->HitSigY[ihit];
	double ulocal = T->HitXlocal[ihit]; //"U" local: generalized "X"
	double vlocal = T->HitYlocal[ihit]; //"V" local: generalized "Y"

	double det = mod_Pxu[module]*mod_Pyv[module] - mod_Pyu[module]*mod_Pxv[module]; //cos( alphau) * sin(alphav) - sin(alphau)*cos(alphav) = 1 for alphau = 0, alphav = 90
	
	double xlocal = (mod_Pyv[module]*ulocal - mod_Pyu[module]*vlocal)/det; //(sin(alphav)*U - sin(alphau)*V)/det = U = X for alphau = 0, alphav = 90
	double ylocal = (mod_Pxu[module]*vlocal - mod_Pxv[module]*ulocal)/det; //(cos(alphau)*V - cos(alphav)*U)/det = V = Y for alphau = 0, alphav = 90

	//On subsequent iterations after the first, we want to fit the changes in the parameters relative to the previous iteration. How can we do this properly?
	//We need to come up with a new definition for the "local" coordinates that properly accounts for the new coordinate system:
	//We already re-fit the track; this means that 

	if( nevent < NMAX && iter == niter ){
	  HITMODTEMP.push_back( module );
	  HITXTEMP.push_back( ulocal );
	  HITYTEMP.push_back( vlocal );
	}
	TVector3 hitpos_local(xlocal,ylocal,0);
	TRotation R;
	R.RotateX( mod_ax[module] );
	R.RotateY( mod_ay[module] );
	R.RotateZ( mod_az[module] );
	
	TVector3 modcenter_global( mod_x0[module],mod_y0[module],mod_z0[module] );
	TVector3 hitpos_global = modcenter_global + R*hitpos_local;
	
	double sigma = 0.1;
	double weight = pow(sigma_hitpos,-2);

	int ipar_fix[3] = {3*module,3*module+1,3*module+2};

	if( refmod >= 0 && refmod < nmodules ){
	  if( module > refmod ){
	  // 	ipar_x0 = 6*(module-1);
	  // 	ipar_y0 = 6*(module-1)+1;
	  // 	ipar_z0 = 6*(module-1)+2;
	  // 	ipar_ax = 6*(module-1)+3;
	  // 	ipar_ay = 6*(module-1)+4;
	  // 	ipar_az = 6*(module-1)+5;
	  
	    ipar_fix[0] = 3*(module-1);
	    ipar_fix[1] = 3*(module-1)+1;
	    ipar_fix[2] = 3*(module-1)+2;
	    //   }
	    // }
	  }
	}

	double xcoeff[6] = {1.0, 0.0, -xptrack, -xptrack*ylocal, xptrack*xlocal, -ylocal };
	double ycoeff[6] = {0.0, 1.0, -yptrack, -yptrack*ylocal, yptrack*xlocal, xlocal };
	
	if(freemodindex.find(module) != freemodindex.end()){
	  int modidx = freemodindex[module];
	  
	  int ipar_x0 = 6*modidx;
	  int ipar_y0 = 6*modidx+1;
	  int ipar_z0 = 6*modidx+2;
	  int ipar_ax = 6*modidx+3;
	  int ipar_ay = 6*modidx+4;
	  int ipar_az = 6*modidx+5;

	  int ipar[6] = {ipar_x0, ipar_y0, ipar_z0, ipar_ax, ipar_ay, ipar_az };
      
	  for( int i=0; i<6; i++ ){
	    for( int j=0; j<6; j++ ){
	      M(ipar[i], ipar[j]) += weight*(xcoeff[i]*xcoeff[j] + ycoeff[i]*ycoeff[j]);
	    }
	    b(ipar[i]) += weight*(xcoeff[i]*(xtrack - xlocal) + ycoeff[i]*(ytrack-ylocal));
	    //b(ipar[i]) += xcoeff[i]*xRHS + ycoeff[i]*yRHS;
	  }
	}
	
	for( int i=0; i<3; i++ ){
	  for( int j=0; j<3; j++ ){
	    Mpos( ipar_fix[i], ipar_fix[j] ) += weight*(xcoeff[i]*xcoeff[j]+ycoeff[i]*ycoeff[j]);
	    Mrot( ipar_fix[i], ipar_fix[j] ) += weight*(xcoeff[i+3]*xcoeff[j+3]+ycoeff[i+3]*ycoeff[j+3]);
	  }
	  //For the positional offsets, we need to subtract the sum of all alphax, alphay, alphaz dependent terms from the RHS:
	  // so this is like -xcoeff[i]*(xcoeff[3]*ax + xcoeff[4]*ay + xcoeff[5]*az)-ycoeff[i]*(ycoeff[3]*ax+ycoeff[4]*ay+ycoeff[5]*az)
	  //For the rotational offsets, the opposite is true
	  bpos( ipar_fix[i] ) += weight*( xcoeff[i]*(xtrack-xlocal - (xcoeff[3]*mod_ax[module]+xcoeff[4]*mod_ay[module]+xcoeff[5]*mod_az[module])) +
					  ycoeff[i]*(ytrack-ylocal - (ycoeff[3]*mod_ax[module]+ycoeff[4]*mod_ay[module]+ycoeff[5]*mod_az[module])) );
	  brot( ipar_fix[i] ) += weight*( xcoeff[i+3]*(xtrack-xlocal - (xcoeff[0]*mod_x0[module]+xcoeff[1]*mod_y0[module]+xcoeff[2]*mod_z0[module])) +
					  ycoeff[i+3]*(ytrack-ylocal - (ycoeff[0]*mod_x0[module]+ycoeff[1]*mod_y0[module]+ycoeff[2]*mod_z0[module])) );
	}
      }
    
    
      if( nevent < NMAX && iter == niter ){
	HITMOD.push_back( HITMODTEMP );
	HITX.push_back( HITXTEMP );
	HITY.push_back( HITYTEMP );
      }
    }
    
    
    // if( offsetsonlyflag  != 0 ){
    //   for( int ipar=0; ipar<nparam; ipar++ ){
    // 	bool isroti = (ipar%6 > 2);
    // 	for( int jpar=0; jpar<nparam; jpar++ ){
    // 	  bool isrotj = (jpar%6 > 2);
    // 	  if( jpar != ipar ){
    // 	    if( isroti || isrotj ){
    // 	      M(ipar,jpar)=0.0;
    // 	      M(jpar,ipar) = 0.0;
    // 	    }
    // 	  }
    // 	}
    // 	if( isroti ){
    // 	  //b(ipar) = 0.0;
    // 	  M(ipar,ipar) = 1.0;
    // 	}
    //   }
      
    // } else if( rotationsonlyflag != 0 ){
    //   for( int ipar=0; ipar<nparam; ipar++ ){
    // 	bool isroti = (ipar%6 > 2);
    // 	for( int jpar=0; jpar<nparam; jpar++ ){
    // 	  bool isrotj = (jpar%6 > 2);
    // 	  if( jpar != ipar ){
    // 	    if( !isroti || !isrotj ){
    // 	      M(ipar,jpar)=0.0;
    // 	      M(jpar,ipar) = 0.0;
    // 	    }
    // 	  }
    // 	}
    // 	if( !isroti ){
    // 	  //b(ipar) = 0.0;
    // 	  M(ipar,ipar) = 1.0;
    // 	}
    //   }
    // }

    // if( fixz > 0 ){
    //   for( int ipar=0; ipar<nparam; ipar++ ){
    // 	bool ismodzi = (ipar%6 == 2);
    // 	for( int jpar=0; jpar<nparam; jpar++ ){
    // 	  bool ismodzj = (jpar%6 == 2);
    // 	  if( jpar != ipar ){
    // 	    if( ismodzi || ismodzj ){
    // 	      M(ipar,jpar)=0.0;
    // 	      M(jpar,ipar)=0.0;
    // 	    }
    // 	  }
    // 	}
    // 	if( ismodzi ){
    // 	  b(ipar) = 1.0;
    // 	  M(ipar,ipar) = 1.0;
    // 	}
    //   } 
    // }

    // if( fixax > 0 ){
    //   for( int ipar=0; ipar<nparam; ipar++ ){
    // 	bool ismodzi = (ipar%6 == 3);
    // 	for( int jpar=0; jpar<nparam; jpar++ ){
    // 	  bool ismodzj = (jpar%6 == 3);
    // 	  if( jpar != ipar ){
    // 	    if( ismodzi || ismodzj ){
    // 	      M(ipar,jpar)=0.0;
    // 	      M(jpar,ipar)=0.0;
    // 	    }
    // 	  }
    // 	}
    // 	if( ismodzi ){
    // 	  b(ipar) = 1.0;
    // 	  M(ipar,ipar) = 1.0;
    // 	}
    //   }
    // }

    // if( fixay > 0 ){
    //   for( int ipar=0; ipar<nparam; ipar++ ){
    // 	bool ismodzi = (ipar%6 == 4);
    // 	for( int jpar=0; jpar<nparam; jpar++ ){
    // 	  bool ismodzj = (jpar%6 == 4);
    // 	  if( jpar != ipar ){
    // 	    if( ismodzi || ismodzj ){
    // 	      M(ipar,jpar)=0.0;
    // 	      M(jpar,ipar)=0.0;
    // 	    }
    // 	  }
    // 	}
    // 	if( ismodzi ){
    // 	  b(ipar) = 0.0;
    // 	  M(ipar,ipar) = 1.0;
    // 	}
    //   }
      
    // }

    // if( fixaz > 0 ){
    //   for( int ipar=0; ipar<nparam; ipar++ ){
    // 	bool ismodzi = (ipar%6 == 5);
    // 	for( int jpar=0; jpar<nparam; jpar++ ){
    // 	  bool ismodzj = (jpar%6 == 5);
    // 	  if( jpar != ipar ){
    // 	    if( ismodzi || ismodzj ){
    // 	      M(ipar,jpar)=0.0;
    // 	      M(jpar,ipar)=0.0;
    // 	    }
    // 	  }
    // 	}
    // 	if( ismodzi ){
    // 	  b(ipar) = 0.0;
    // 	  M(ipar,ipar) = 1.0;
    // 	}
    //   }
      
    // }
  
    
    //Okay, wish me luck:

    //M.Print();
    //b.Print();

    cout << "Matrix symmetric? = " << M.IsSymmetric() << endl;

    // Here is an idea: to fit the changes of the parameters instead of the parameters themselves, we subtract from the RHS another vector:
    TVectorD PreviousSolution(nparam);
    TVectorD PreviousSolution_posonly(nparam_pos);
    TVectorD PreviousSolution_rotonly(nparam_rot);
    for( int imodule=0; imodule<nmodules; imodule++ ){

      if( freemodindex.find(imodule) != freemodindex.end() ){
	int modidx = freemodindex[imodule];
      
	int ipar_x0 = modidx*6;
	int ipar_y0 = modidx*6+1;
	int ipar_z0 = modidx*6+2;
	int ipar_ax = modidx*6+3;
	int ipar_ay = modidx*6+4;
	int ipar_az = modidx*6+5;

	PreviousSolution(ipar_x0) = mod_x0[imodule];
	PreviousSolution(ipar_y0) = mod_y0[imodule];
	PreviousSolution(ipar_z0) = mod_z0[imodule];
	PreviousSolution(ipar_ax) = mod_ax[imodule];
	PreviousSolution(ipar_ay) = mod_ay[imodule];
	PreviousSolution(ipar_az) = mod_az[imodule];
	
      }
      int iparx = 3*imodule;
      int ipary = 3*imodule+1;
      int iparz = 3*imodule+2;
      
     
      if( refmod >= 0 && refmod < nmodules && imodule > refmod ){
	// ipar_x0 = 6*(imodule-1);
	// ipar_y0 = 6*(imodule-1)+1;
	// ipar_z0 = 6*(imodule-1)+2;
	// ipar_ax = 6*(imodule-1)+3;
	// ipar_ay = 6*(imodule-1)+4;
	// ipar_az = 6*(imodule-1)+5;

	iparx = 3*(imodule-1);
	ipary = 3*(imodule-1)+1;
	iparz = 3*(imodule-1)+2;
      }

      if( imodule != refmod ){


	PreviousSolution_posonly(iparx) = mod_x0[imodule];
	PreviousSolution_posonly(ipary) = mod_y0[imodule];
	PreviousSolution_posonly(iparz) = mod_z0[imodule];

	PreviousSolution_rotonly(iparx) = mod_ax[imodule];
	PreviousSolution_rotonly(ipary) = mod_ay[imodule];
	PreviousSolution_rotonly(iparz) = mod_az[imodule];
      }

     
      
    }

    TVectorD bshift = b - M*PreviousSolution;
    TVectorD bshiftpos = bpos - Mpos*PreviousSolution_posonly;
    TVectorD bshiftrot = brot - Mrot*PreviousSolution_rotonly;
    
    TVectorD Solution(nparam);
    TVectorD Solution_posonly(nparam_pos);
    TVectorD Solution_rotonly(nparam_rot);
    
    M.Invert();
    Mpos.Invert();
    Mrot.Invert();
    
    if( iter >= 0 ){
      Solution = M*bshift;
      Solution_posonly = Mpos*bshiftpos;
      Solution_rotonly = Mrot*bshiftrot;
    } else {
      Solution = M*b;
      Solution_posonly = Mpos * bpos;
      Solution_rotonly = Mrot * brot;
    }

    
    
    //M.Print();
    
    //Solution.Print();

    map<int,double> prev_x0 = mod_x0;
    map<int,double> prev_y0 = mod_y0;
    map<int,double> prev_z0 = mod_z0;
    map<int,double> prev_ax = mod_ax;
    map<int,double> prev_ay = mod_ay;
    map<int,double> prev_az = mod_az;

    double startpar[6*nmodules];
  
    for( int imodule=0; imodule<nmodules; imodule++ ){
      if( freemodindex.find(imodule) != freemodindex.end() ){
	int modidx = freemodindex[imodule];
	int ipar_x0 = modidx*6;
	int ipar_y0 = modidx*6+1;
	int ipar_z0 = modidx*6+2;
	int ipar_ax = modidx*6+3;
	int ipar_ay = modidx*6+4;
	int ipar_az = modidx*6+5;

	mod_x0[imodule] = Solution(ipar_x0);
	mod_y0[imodule] = Solution(ipar_y0);
	mod_z0[imodule] = Solution(ipar_z0);
	mod_ax[imodule] = Solution(ipar_ax);
	mod_ay[imodule] = Solution(ipar_ay);
	mod_az[imodule] = Solution(ipar_az);

	if( iter >= 0 ){
	  mod_x0[imodule] += prev_x0[imodule];
	  mod_y0[imodule] += prev_y0[imodule];
	  mod_z0[imodule] += prev_z0[imodule];
	  mod_ax[imodule] += prev_ax[imodule];
	  mod_ay[imodule] += prev_ay[imodule];
	  mod_az[imodule] += prev_az[imodule];
	}

	//Strictly speaking, I don't think this is necessary:
	PreviousSolution(ipar_x0) = mod_x0[imodule];
	PreviousSolution(ipar_y0) = mod_y0[imodule];
	PreviousSolution(ipar_z0) = mod_z0[imodule];
	PreviousSolution(ipar_ax) = mod_ax[imodule];
	PreviousSolution(ipar_ay) = mod_ay[imodule];
	PreviousSolution(ipar_az) = mod_az[imodule];
	
	//	}
	mod_dx0[imodule] = sqrt(fabs(M(ipar_x0,ipar_x0)));
	mod_dy0[imodule] = sqrt(fabs(M(ipar_y0,ipar_y0)));
	mod_dz0[imodule] = sqrt(fabs(M(ipar_z0,ipar_z0)));
	mod_dax[imodule] = sqrt(fabs(M(ipar_ax,ipar_ax)));
	mod_day[imodule] = sqrt(fabs(M(ipar_ay,ipar_ay)));
	mod_daz[imodule] = sqrt(fabs(M(ipar_az,ipar_az)));
	
      }
      
      int iparx = 3*imodule;
      int ipary = 3*imodule+1;
      int iparz = 3*imodule+2;
      
      if( refmod >= 0 && refmod < nmodules && imodule > refmod ){
	// ipar_x0 = 6*(imodule-1);
	// ipar_y0 = 6*(imodule-1)+1;
	// ipar_z0 = 6*(imodule-1)+2;
	// ipar_ax = 6*(imodule-1)+3;
	// ipar_ay = 6*(imodule-1)+4;
	// ipar_az = 6*(imodule-1)+5;

	iparx = 3*(imodule-1);
	ipary = 3*(imodule-1)+1;
	iparz = 3*(imodule-1)+2;
      }

      //now the solution represents the CHANGE in each parameter from the previous iteration or, on the first iteration the starting solution:
      
      if( imodule != refmod ){
	
	
	if( offsetsonlyflag != 0 ){
	  mod_x0[imodule] = Solution_posonly(iparx) + prev_x0[imodule];
	  mod_y0[imodule] = Solution_posonly(ipary) + prev_y0[imodule];
	  mod_z0[imodule] = Solution_posonly(iparz) + prev_z0[imodule];
	  
	  mod_ax[imodule] = prev_ax[imodule];
	  mod_ay[imodule] = prev_ay[imodule];
	  mod_az[imodule] = prev_az[imodule];
	}

	if( rotationsonlyflag != 0 ){
	  mod_x0[imodule] = prev_x0[imodule];
	  mod_y0[imodule] = prev_y0[imodule];
	  mod_z0[imodule] = prev_z0[imodule];

	  mod_ax[imodule] = Solution_rotonly(iparx) + prev_ax[imodule];
	  mod_ay[imodule] = Solution_rotonly(ipary) + prev_ay[imodule];
	  mod_az[imodule] = Solution_rotonly(iparz) + prev_az[imodule];
	  
	}

	
	
	
      }

      for( int ipar=0; ipar<6; ipar++ ){
	startpar[imodule*6+ipar] = 0.0;
      }
    }
    
    cout << "ending solution: " << endl;
    for( map<int,double>::iterator imod=mod_x0.begin(); imod!=mod_x0.end(); ++imod ){
      int module = imod->first;
      cout << "Module " << module << ": (x0,y0,z0,ax,ay,az)=("
	   << mod_x0[module] << ", " << mod_y0[module] << ", " << mod_z0[module] << ", "
	   << mod_ax[module] << ", " << mod_ay[module] << ", " << mod_az[module] << ")" << endl;
    }
    for( map<int,double>::iterator imod=mod_x0.begin(); imod!=mod_x0.end(); ++imod ){
      int module = imod->first;
      cout << "(Change from previous)/sigma: (dx0,dy0,dz0,dax,day,daz)=("
	   << (mod_x0[module]-prev_x0[module])/mod_dx0[module] << ", "
	   << (mod_y0[module]-prev_y0[module])/mod_dy0[module] << ", "
	   << (mod_z0[module]-prev_z0[module])/mod_dz0[module] << ", "
	   << (mod_ax[module]-prev_ax[module])/mod_dax[module] << ", "
	   << (mod_ay[module]-prev_ay[module])/mod_day[module] << ", "
	   << (mod_az[module]-prev_az[module])/mod_daz[module] << ")" << endl;
    }
  }
  
  if( (offsetsonlyflag == 0 && rotationsonlyflag == 0) ){
    
    TMinuit *ExtraFit = new TMinuit( 6*nmodules );

    ExtraFit->SetFCN( CHI2_FCN );

    int ierflg=0;
  
    for( int mod=0; mod<nmodules; mod++ ){
      TString sparname;

      ExtraFit->mnparm( 6*mod,   sparname.Format( "mod%d_dx0", mod ), 0.0, mod_dx0[mod],0,0,ierflg);
      ExtraFit->mnparm( 6*mod+1, sparname.Format( "mod%d_dy0", mod ), 0.0, mod_dy0[mod],0,0,ierflg);
      ExtraFit->mnparm( 6*mod+2, sparname.Format( "mod%d_dz0", mod ), 0.0, mod_dz0[mod],0,0,ierflg);
      ExtraFit->mnparm( 6*mod+3, sparname.Format( "mod%d_dax", mod ), 0.0, mod_dax[mod],0,0,ierflg);
      ExtraFit->mnparm( 6*mod+4, sparname.Format( "mod%d_day", mod ), 0.0, mod_day[mod],0,0,ierflg);
      ExtraFit->mnparm( 6*mod+5, sparname.Format( "mod%d_daz", mod ), 0.0, mod_daz[mod],0,0,ierflg);
    }

  
    double arglist[10];
    arglist[0]=1;
    ExtraFit->mnexcm("SET ERR",arglist,1,ierflg);

    arglist[0] = 5000;
    arglist[1] = 1.;
  
    ExtraFit->mnexcm("MIGRAD",arglist,2,ierflg);

    for( int mod=0; mod<nmodules; mod++ ){
      double dummy;
      double dx0,dy0,dz0,dax,day,daz;
      ExtraFit->GetParameter(6*mod,dx0,dummy);
      ExtraFit->GetParameter(6*mod+1,dy0,dummy);
      ExtraFit->GetParameter(6*mod+2,dz0,dummy);
      ExtraFit->GetParameter(6*mod+3,dax,dummy);
      ExtraFit->GetParameter(6*mod+4,day,dummy);
      ExtraFit->GetParameter(6*mod+5,daz,dummy);
      mod_x0[mod] += dx0;
      mod_y0[mod] += dy0;
      mod_z0[mod] += dz0;
      mod_ax[mod] += dax;
      mod_ay[mod] += day;
      mod_az[mod] += daz;
    
    }

  }
  
  ofstream outfile(outputfilename);

  //outfile << "mod_x0 ";
  TString x0line = "mod_x0 ", y0line="mod_y0 ", z0line = "mod_z0 ";
  TString axline = "mod_ax ", ayline="mod_ay ", azline = "mod_az ";
  for( map<int,double>::iterator imod=mod_x0.begin(); imod!=mod_x0.end(); ++imod ){
    int module = imod->first;
    TString stemp;
    stemp.Form( " %15.6g", mod_x0[module] );
    x0line += stemp;
    stemp.Form( " %15.6g", mod_y0[module] );
    y0line += stemp;
    stemp.Form( " %15.6g", mod_z0[module] );
    z0line += stemp;
    stemp.Form( " %15.6g", mod_ax[module] );
    axline += stemp;
    stemp.Form( " %15.6g", mod_ay[module] );
    ayline += stemp;
    stemp.Form( " %15.6g", mod_az[module] );
    azline += stemp;
  }
  outfile << x0line << endl;
  outfile << y0line << endl;
  outfile << z0line << endl;
  outfile << axline << endl;
  outfile << ayline << endl;
  outfile << azline << endl;

  
  
  elist->Delete();

}
