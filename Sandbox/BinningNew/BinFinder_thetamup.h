//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb  1 21:27:21 2018 by ROOT version 6.06/06
// from TTree mom_tree/mom_tree
// found on file: ./ubxsecana_output_bnbcosmic_mcc8.6.root
//////////////////////////////////////////////////////////

#ifndef BinFinder_thetamup_h
#define BinFinder_thetamup_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class BinFinder_thetamup {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        mom_tree_true;
   Double_t        mom_tree_mcs;
   Bool_t          mom_tree_contained;
   Bool_t          mom_tree_selected;

   // List of branches
   TBranch        *b_mom_tree_true;   //!
   TBranch        *b_mom_tree_mcs;   //!
   TBranch        *b_mom_tree_contained;   //!
   TBranch        *b_mom_tree_selected;   //!

   BinFiner_thetamup(TTree *tree=0);
   virtual ~BinFiner_thetamup();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void DoFit(TH1D*, std::string, std::string, double *, bool save_plot);
   void SmearingMatrix();
};

#endif

#ifdef BinFiner_thetamup_cxx
BinFiner_thetamup::BinFiner_thetamup(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../../Files/Output/ubxsecana_output_bnbcosmic_mcc8.9_500k.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../../Files/Output/ubxsecana_output_bnbcosmic_mcc8.9_500k.root");
      }
      f->GetObject("mom_tree",tree);

   }
   Init(tree);
}

BinFiner_thetamup::~BinFiner_thetamup()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BinFiner_thetamup::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BinFiner_thetamup::LoadTree(Long64_t entry)
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

void BinFiner_thetamup::Init(TTree *tree)
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

   fChain->SetBranchAddress("mom_tree_true", &mom_tree_true, &b_mom_tree_true);
   fChain->SetBranchAddress("mom_tree_mcs", &mom_tree_mcs, &b_mom_tree_mcs);
   fChain->SetBranchAddress("mom_tree_contained", &mom_tree_contained, &b_mom_tree_contained);
   fChain->SetBranchAddress("mom_tree_selected", &mom_tree_selected, &b_mom_tree_selected);
   Notify();
}

Bool_t BinFiner_thetamup::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BinFiner_thetamup::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BinFiner_thetamup::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BinFiner_thetamup_cxx
