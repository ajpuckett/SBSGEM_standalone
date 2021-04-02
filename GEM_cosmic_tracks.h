//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr  2 12:36:52 2021 by ROOT version 6.22/08
// from TTree Tout/INFN GEM 4-layer cosmic tracks
// found on file: temp3805_allevents.root
//////////////////////////////////////////////////////////

#ifndef GEM_cosmic_tracks_h
#define GEM_cosmic_tracks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class GEM_cosmic_tracks {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           EventID;
   Int_t           Ntracks;
   Double_t        CALOsum;
   Int_t           NGOODSCINT;
   Double_t        TrackXp;
   Double_t        TrackYp;
   Double_t        TrackX;
   Double_t        TrackY;
   Double_t        TrackChi2NDF;
   Int_t           TrackNhits;
   Int_t           Nlayers;
   Int_t           Ncluster[4];   //[Nlayers]
   Int_t           HitModule[4];   //[TrackNhits]
   Int_t           HitLayer[4];   //[TrackNhits]
   Double_t        HitXlocal[4];   //[TrackNhits]
   Double_t        HitYlocal[4];   //[TrackNhits]
   Double_t        HitXglobal[4];   //[TrackNhits]
   Double_t        HitYglobal[4];   //[TrackNhits]
   Double_t        HitZglobal[4];   //[TrackNhits]
   Double_t        HitXresid[4];   //[TrackNhits]
   Double_t        HitYresid[4];   //[TrackNhits]
   Double_t        HitXresidE[4];   //[TrackNhits]
   Double_t        HitYresidE[4];   //[TrackNhits]
   Double_t        HitSigX[4];   //[TrackNhits]
   Double_t        HitSigY[4];   //[TrackNhits]
   Double_t        HitXmom[4];   //[TrackNhits]
   Double_t        HitYmom[4];   //[TrackNhits]
   Double_t        HitADCX[4];   //[TrackNhits]
   Double_t        HitADCY[4];   //[TrackNhits]
   Double_t        HitADCasym[4];   //[TrackNhits]
   Double_t        HitTmean[4];   //[TrackNhits]
   Double_t        HitdT[4];   //[TrackNhits]
   Double_t        HitCorrCoeff[4];   //[TrackNhits]
   Double_t        StripMaxCorrCoeff[4];   //[TrackNhits]
   Int_t           HitNstripX[4];   //[TrackNhits]
   Int_t           HitNstripY[4];   //[TrackNhits]
   Int_t           HitXstripMax[4];   //[TrackNhits]
   Int_t           HitYstripMax[4];   //[TrackNhits]
   Int_t           HitXstripLo[4];   //[TrackNhits]
   Int_t           HitYstripLo[4];   //[TrackNhits]
   Int_t           HitXstripHi[4];   //[TrackNhits]
   Int_t           HitYstripHi[4];   //[TrackNhits]

   // List of branches
   TBranch        *b_EventID;   //!
   TBranch        *b_Ntracks;   //!
   TBranch        *b_CALOsum;   //!
   TBranch        *b_NGOODSCINT;   //!
   TBranch        *b_TrackXp;   //!
   TBranch        *b_TrackYp;   //!
   TBranch        *b_TrackX;   //!
   TBranch        *b_TrackY;   //!
   TBranch        *b_TrackChi2NDF;   //!
   TBranch        *b_TrackNhits;   //!
   TBranch        *b_Nlayers;   //!
   TBranch        *b_Ncluster;   //!
   TBranch        *b_HitModule;   //!
   TBranch        *b_HitLayer;   //!
   TBranch        *b_HitXlocal;   //!
   TBranch        *b_HitYlocal;   //!
   TBranch        *b_HitXglobal;   //!
   TBranch        *b_HitYglobal;   //!
   TBranch        *b_HitZglobal;   //!
   TBranch        *b_HitXresid;   //!
   TBranch        *b_HitYresid;   //!
   TBranch        *b_HitXresidE;   //!
   TBranch        *b_HitYresidE;   //!
   TBranch        *b_HitSigX;   //!
   TBranch        *b_HitSigY;   //!
   TBranch        *b_HitXmom;   //!
   TBranch        *b_HitYmom;   //!
   TBranch        *b_HitADCX;   //!
   TBranch        *b_HitADCY;   //!
   TBranch        *b_HitADCasym;   //!
   TBranch        *b_HitTmean;   //!
   TBranch        *b_HitdT;   //!
   TBranch        *b_HitCorrCoeff;   //!
   TBranch        *b_StripMaxCorrCoeff;   //!
   TBranch        *b_HitNstripX;   //!
   TBranch        *b_HitNstripY;   //!
   TBranch        *b_HitXstripMax;   //!
   TBranch        *b_HitYstripMax;   //!
   TBranch        *b_HitXstripLo;   //!
   TBranch        *b_HitYstripLo;   //!
   TBranch        *b_HitXstripHi;   //!
   TBranch        *b_HitYstripHi;   //!

   GEM_cosmic_tracks(TTree *tree=0);
   virtual ~GEM_cosmic_tracks();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GEM_cosmic_tracks_cxx
GEM_cosmic_tracks::GEM_cosmic_tracks(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("temp3805_allevents.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("temp3805_allevents.root");
      }
      f->GetObject("Tout",tree);

   }
   Init(tree);
}

GEM_cosmic_tracks::~GEM_cosmic_tracks()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GEM_cosmic_tracks::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GEM_cosmic_tracks::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void GEM_cosmic_tracks::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventID", &EventID, &b_EventID);
   fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
   fChain->SetBranchAddress("CALOsum", &CALOsum, &b_CALOsum);
   fChain->SetBranchAddress("NGOODSCINT", &NGOODSCINT, &b_NGOODSCINT);
   fChain->SetBranchAddress("TrackXp", &TrackXp, &b_TrackXp);
   fChain->SetBranchAddress("TrackYp", &TrackYp, &b_TrackYp);
   fChain->SetBranchAddress("TrackX", &TrackX, &b_TrackX);
   fChain->SetBranchAddress("TrackY", &TrackY, &b_TrackY);
   fChain->SetBranchAddress("TrackChi2NDF", &TrackChi2NDF, &b_TrackChi2NDF);
   fChain->SetBranchAddress("TrackNhits", &TrackNhits, &b_TrackNhits);
   fChain->SetBranchAddress("Nlayers", &Nlayers, &b_Nlayers);
   fChain->SetBranchAddress("Ncluster", Ncluster, &b_Ncluster);
   fChain->SetBranchAddress("HitModule", HitModule, &b_HitModule);
   fChain->SetBranchAddress("HitLayer", HitLayer, &b_HitLayer);
   fChain->SetBranchAddress("HitXlocal", HitXlocal, &b_HitXlocal);
   fChain->SetBranchAddress("HitYlocal", HitYlocal, &b_HitYlocal);
   fChain->SetBranchAddress("HitXglobal", HitXglobal, &b_HitXglobal);
   fChain->SetBranchAddress("HitYglobal", HitYglobal, &b_HitYglobal);
   fChain->SetBranchAddress("HitZglobal", HitZglobal, &b_HitZglobal);
   fChain->SetBranchAddress("HitXresid", HitXresid, &b_HitXresid);
   fChain->SetBranchAddress("HitYresid", HitYresid, &b_HitYresid);
   fChain->SetBranchAddress("HitXresidE", HitXresidE, &b_HitXresidE);
   fChain->SetBranchAddress("HitYresidE", HitYresidE, &b_HitYresidE);
   fChain->SetBranchAddress("HitSigX", HitSigX, &b_HitSigX);
   fChain->SetBranchAddress("HitSigY", HitSigY, &b_HitSigY);
   fChain->SetBranchAddress("HitXmom", HitXmom, &b_HitXmom);
   fChain->SetBranchAddress("HitYmom", HitYmom, &b_HitYmom);
   fChain->SetBranchAddress("HitADCX", HitADCX, &b_HitADCX);
   fChain->SetBranchAddress("HitADCY", HitADCY, &b_HitADCY);
   fChain->SetBranchAddress("HitADCasym", HitADCasym, &b_HitADCasym);
   fChain->SetBranchAddress("HitTmean", HitTmean, &b_HitTmean);
   fChain->SetBranchAddress("HitdT", HitdT, &b_HitdT);
   fChain->SetBranchAddress("HitCorrCoeff", HitCorrCoeff, &b_HitCorrCoeff);
   fChain->SetBranchAddress("StripMaxCorrCoeff", StripMaxCorrCoeff, &b_StripMaxCorrCoeff);
   fChain->SetBranchAddress("HitNstripX", HitNstripX, &b_HitNstripX);
   fChain->SetBranchAddress("HitNstripY", HitNstripY, &b_HitNstripY);
   fChain->SetBranchAddress("HitXstripMax", HitXstripMax, &b_HitXstripMax);
   fChain->SetBranchAddress("HitYstripMax", HitYstripMax, &b_HitYstripMax);
   fChain->SetBranchAddress("HitXstripLo", HitXstripLo, &b_HitXstripLo);
   fChain->SetBranchAddress("HitYstripLo", HitYstripLo, &b_HitYstripLo);
   fChain->SetBranchAddress("HitXstripHi", HitXstripHi, &b_HitXstripHi);
   fChain->SetBranchAddress("HitYstripHi", HitYstripHi, &b_HitYstripHi);
   Notify();
}

Bool_t GEM_cosmic_tracks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GEM_cosmic_tracks::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GEM_cosmic_tracks::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GEM_cosmic_tracks_cxx
