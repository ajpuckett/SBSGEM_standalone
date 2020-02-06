//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 15 11:27:01 2020 by ROOT version 6.18/00
// from TTree GEMHit/Hit list
// found on file: ../HallADATA/run288.root
//////////////////////////////////////////////////////////

#ifndef GEMHit_tree_HallAtest_h
#define GEMHit_tree_HallAtest_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class GEMHit_tree_HallAtest {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           evtID;
   Int_t           nch;
   Int_t           strip[660];   //[nch]
   Int_t           detID[660];   //[nch]
   Int_t           planeID[660];   //[nch]
   Int_t           adc0[660];   //[nch]
   Int_t           adc1[660];   //[nch]
   Int_t           adc2[660];   //[nch]
   Int_t           adc3[660];   //[nch]
   Int_t           adc4[660];   //[nch]
   Int_t           adc5[660];   //[nch]
   Int_t           nfadc;
   Int_t           fCH[20];   //[nfadc]
   Int_t           ftimting[20];   //[nfadc]
   Int_t           fadc[20];   //[nfadc]
   Int_t           ntdc;
   Int_t           tCH[19];   //[ntdc]
   Double_t        ttiming[19];   //[ntdc]

   // List of branches
   TBranch        *b_evtID;   //!
   TBranch        *b_nch;   //!
   TBranch        *b_strip;   //!
   TBranch        *b_detID;   //!
   TBranch        *b_planeID;   //!
   TBranch        *b_adc0;   //!
   TBranch        *b_adc1;   //!
   TBranch        *b_adc2;   //!
   TBranch        *b_adc3;   //!
   TBranch        *b_adc4;   //!
   TBranch        *b_adc5;   //!
   TBranch        *b_nfadc;   //!
   TBranch        *b_fCH;   //!
   TBranch        *b_ftimting;   //!
   TBranch        *b_fadc;   //!
   TBranch        *b_ntdc;   //!
   TBranch        *b_tCH;   //!
   TBranch        *b_ttiming;   //!

   GEMHit_tree_HallAtest(TTree *tree=0);
   virtual ~GEMHit_tree_HallAtest();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GEMHit_tree_HallAtest_cxx
GEMHit_tree_HallAtest::GEMHit_tree_HallAtest(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../HallADATA/run288.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../HallADATA/run288.root");
      }
      f->GetObject("GEMHit",tree);

   }
   Init(tree);
}

GEMHit_tree_HallAtest::~GEMHit_tree_HallAtest()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GEMHit_tree_HallAtest::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GEMHit_tree_HallAtest::LoadTree(Long64_t entry)
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

void GEMHit_tree_HallAtest::Init(TTree *tree)
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

   fChain->SetBranchAddress("evtID", &evtID, &b_evtID);
   fChain->SetBranchAddress("nch", &nch, &b_nch);
   fChain->SetBranchAddress("strip", strip, &b_strip);
   fChain->SetBranchAddress("detID", detID, &b_detID);
   fChain->SetBranchAddress("planeID", planeID, &b_planeID);
   fChain->SetBranchAddress("adc0", adc0, &b_adc0);
   fChain->SetBranchAddress("adc1", adc1, &b_adc1);
   fChain->SetBranchAddress("adc2", adc2, &b_adc2);
   fChain->SetBranchAddress("adc3", adc3, &b_adc3);
   fChain->SetBranchAddress("adc4", adc4, &b_adc4);
   fChain->SetBranchAddress("adc5", adc5, &b_adc5);
   fChain->SetBranchAddress("nfadc", &nfadc, &b_nfadc);
   fChain->SetBranchAddress("fCH", fCH, &b_fCH);
   fChain->SetBranchAddress("ftimting", ftimting, &b_ftimting);
   fChain->SetBranchAddress("fadc", fadc, &b_fadc);
   fChain->SetBranchAddress("ntdc", &ntdc, &b_ntdc);
   fChain->SetBranchAddress("tCH", tCH, &b_tCH);
   fChain->SetBranchAddress("ttiming", ttiming, &b_ttiming);
   Notify();
}

Bool_t GEMHit_tree_HallAtest::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GEMHit_tree_HallAtest::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GEMHit_tree_HallAtest::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GEMHit_tree_HallAtest_cxx
