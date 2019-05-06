//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Feb 10 19:54:33 2019 by ROOT version 6.06/08
// from TTree true_reco_tree/true_reco_tree
// found on file: /uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb10/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22.root
//////////////////////////////////////////////////////////

#ifndef BinFinder_h
#define BinFinder_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class BinFinder {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        mom_true;
   Double_t        mom_mcs;
   Double_t        pmom_true;
   Double_t        pmom_reco;
   Bool_t          contained;
   Bool_t          selected;
   Double_t        angle_true;
   Double_t        angle_reco;
   Double_t        event_weight;
   vector<string>  *wgtsnames_genie_multisim;
   vector<double>  *wgts_genie_multisim;
   vector<string>  *wgtsnames_extra_syst;
   vector<double>  *wgts_extra_syst;
   vector<string>  *wgtsnames_flux_multisim;
   vector<double>  *wgts_flux_multisim;

   // List of branches
   TBranch        *b_mom_true;   //!
   TBranch        *b_mom_mcs;   //!
   TBranch        *b_pmom_true;   //!
   TBranch        *b_pmom_reco;   //!
   TBranch        *b_contained;   //!
   TBranch        *b_selected;   //!
   TBranch        *b_angle_true;   //!
   TBranch        *b_angle_reco;   //!
   TBranch        *b_event_weight;   //!
   TBranch        *b_wgtsnames_genie_multisim;   //!
   TBranch        *b_wgts_genie_multisim;   //!
   TBranch        *b_wgtsnames_extra_syst;   //!
   TBranch        *b_wgts_extra_syst;   //!
   TBranch        *b_wgtsnames_flux_multisim;   //!
   TBranch        *b_wgts_flux_multisim;   //!

   BinFinder(TTree *tree=0);
   virtual ~BinFinder();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void DoFit(TH1D*, std::string, std::string, double *, bool save_plot);

};

#endif

#ifdef BinFinder_cxx
BinFinder::BinFinder(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb10/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb10/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22.root");
      }
      f->GetObject("true_reco_tree",tree);

   }
   Init(tree);
}

BinFinder::~BinFinder()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BinFinder::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BinFinder::LoadTree(Long64_t entry)
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

void BinFinder::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   wgtsnames_genie_multisim = 0;
   wgts_genie_multisim = 0;
   wgtsnames_extra_syst = 0;
   wgts_extra_syst = 0;
   wgtsnames_flux_multisim = 0;
   wgts_flux_multisim = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mom_true", &mom_true, &b_mom_true);
   fChain->SetBranchAddress("mom_mcs", &mom_mcs, &b_mom_mcs);
   fChain->SetBranchAddress("pmom_true", &pmom_true, &b_pmom_true);
   fChain->SetBranchAddress("pmom_reco", &pmom_reco, &b_pmom_reco);
   fChain->SetBranchAddress("contained", &contained, &b_contained);
   fChain->SetBranchAddress("selected", &selected, &b_selected);
   fChain->SetBranchAddress("angle_true", &angle_true, &b_angle_true);
   fChain->SetBranchAddress("angle_reco", &angle_reco, &b_angle_reco);
   fChain->SetBranchAddress("event_weight", &event_weight, &b_event_weight);
   fChain->SetBranchAddress("wgtsnames_genie_multisim", &wgtsnames_genie_multisim, &b_wgtsnames_genie_multisim);
   fChain->SetBranchAddress("wgts_genie_multisim", &wgts_genie_multisim, &b_wgts_genie_multisim);
   fChain->SetBranchAddress("wgtsnames_extra_syst", &wgtsnames_extra_syst, &b_wgtsnames_extra_syst);
   fChain->SetBranchAddress("wgts_extra_syst", &wgts_extra_syst, &b_wgts_extra_syst);
   fChain->SetBranchAddress("wgtsnames_flux_multisim", &wgtsnames_flux_multisim, &b_wgtsnames_flux_multisim);
   fChain->SetBranchAddress("wgts_flux_multisim", &wgts_flux_multisim, &b_wgts_flux_multisim);
   Notify();
}

Bool_t BinFinder::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BinFinder::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BinFinder::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BinFinder_cxx
