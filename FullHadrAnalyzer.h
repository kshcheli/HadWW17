//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct  8 21:28:41 2018 by ROOT version 6.10/09
// from TChain ntp1/
//////////////////////////////////////////////////////////

#ifndef FullHadrAnalyzer_h
#define FullHadrAnalyzer_h


//#define SIG
//#define PYT
//#define DAT


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

using namespace std;
// Header file for the classes stored in the TTree if any.
#include "vector"
//#include "vector"
//#include "vector"

class FullHadrAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *jet_pt;
   vector<float>   *jet_energy;
   vector<float>   *jet_phi;
   vector<float>   *jet_eta;
   vector<float>   *jet_mass;
   vector<float>   *jet_tau1;
   vector<float>   *jet_tau2;
   vector<float>   *jet_corrmass;
   vector<float>   *jet_vertexz;
   vector<float>   *dijet_mass;
   vector<float>   *dijet_pt;
   vector<float>   *dijet_y;
   vector<float>   *dijet_phi;
   vector<float>   *dijet_dphi;
   vector<float>   *pps_track_x;
   vector<float>   *pps_track_y;
   vector<int>     *pps_track_rpid;
   vector<string>  *hlt;
#ifndef DAT
   vector<float>   *gen_proton_px;
   vector<float>   *gen_proton_py;
   vector<float>   *gen_proton_pz;
   vector<float>   *gen_proton_xi;
   vector<float>   *gen_proton_t;
   vector<float>   *gen_jet_pt;
   vector<float>   *gen_jet_eta;
   vector<float>   *gen_jet_phi;
   vector<float>   *gen_jet_energy;
   vector<float>   *gen_dijet_mass;
   vector<float>   *gen_dijet_y;
   vector<float>   *jet_jer_res;
   vector<float>   *jet_jer_sf;
   vector<float>   *jet_jer_sfup;
   vector<float>   *jet_jer_sfdown;
   Int_t           gen_n_e;
   Int_t           gen_n_mu;
   Int_t           gen_n_tau;
   UInt_t          mc_pu_trueinteractions_;
#endif
  Float_t         pileupWeight;
 
   UInt_t          nVertices;
   Int_t           run;
   Long64_t        event;
   Int_t           lumiblock;





   // List of branches
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_energy;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_mass;   //!
   TBranch        *b_jet_tau1;   //!
   TBranch        *b_jet_tau2;   //!
   TBranch        *b_jet_corrmass;   //!
   TBranch        *b_jet_vertexz;   //!
   TBranch        *b_dijet_mass;   //!
   TBranch        *b_dijet_pt;   //!
   TBranch        *b_dijet_y;   //!
   TBranch        *b_dijet_phi;   //!
   TBranch        *b_dijet_dphi;   //!
   TBranch        *b_pps_track_x;   //!
   TBranch        *b_pps_track_y;   //!
   TBranch        *b_pps_track_rpid;   //!
   TBranch        *b_hlt;   //!
#ifndef DAT
   TBranch        *b_gen_proton_px;   //!
   TBranch        *b_gen_proton_py;   //!
   TBranch        *b_gen_proton_pz;   //!
   TBranch        *b_gen_proton_xi;   //!
   TBranch        *b_gen_proton_t;   //!
   TBranch        *b_gen_jet_pt;   //!
   TBranch        *b_gen_jet_eta;   //!
   TBranch        *b_gen_jet_phi;   //!
   TBranch        *b_gen_jet_energy;   //!
   TBranch        *b_gen_dijet_mass;   //!
   TBranch        *b_gen_dijet_y;   //!
   TBranch        *b_jet_jer_res;   //!
   TBranch        *b_jet_jer_sf;   //!
   TBranch        *b_jet_jer_sfup;   //!
   TBranch        *b_jet_jer_sfdown;   //!
   TBranch        *b_gen_n_e;   //!
   TBranch        *b_gen_n_mu;   //!
   TBranch        *b_gen_n_tau;   //!
   TBranch        *b_mc_pu_trueinteractions;   //!
#endif
   TBranch        *b_pileupWeight;   //!
   TBranch        *b_nVertices;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!




   FullHadrAnalyzer(TTree *tree=0);
   virtual ~FullHadrAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef FullHadrAnalyzer_cxx
FullHadrAnalyzer::FullHadrAnalyzer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("ntp1",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("ntp1","");
#ifdef SIG 
      chain->Add("/eos/cms/store/user/jjhollar/ExclWWHadronic2017Analysis/WWhadronic_ExclusiveWW_AllDecays_a0W2point5e-6_xi1to30percent_ntuplesv2.root/demo/ntp1");
#endif
      
#ifdef DAT
//      chain->Add("/eos/cms/store/user/kshcheli/ExclWWHadronic2017Analysis/Data17/JetHT/ExclWWjets_Run2017B-31Mar2018-v1_all.root/demo/ntp1");
      chain->Add("/eos/cms/store/user/kshcheli/ExclWWHadronic2017Analysis/Data17/JetHT/ExclWWjets_Run2017C-31Mar2018-v1_all.root/demo/ntp1");
//      chain->Add("/eos/cms/store/user/kshcheli/ExclWWHadronic2017Analysis/Data17/JetHT/ExclWWjets_Run2017D-31Mar2018-v1_all.root/demo/ntp1");
//      chain->Add("/eos/cms/store/user/kshcheli/ExclWWHadronic2017Analysis/Data17/JetHT/ExclWWjets_Run2017E-31Mar2018-v1_all.root/demo/ntp1");
//      chain->Add("/eos/cms/store/user/kshcheli/ExclWWHadronic2017Analysis/Data17/JetHT/ExclWWjets_Run2017F-31Mar2018-v1_all.root/demo/ntp1");
#endif

#ifdef PYT
      chain->Add("/eos/cms/store/user/jjhollar/ExclWWHadronic2017Analysis/WWhadronic_QCDPt170to300_Pythia8_merge_12Apr2018_ntuplesv2.root/demo/ntp1");
      chain->Add("/eos/cms/store/user/jjhollar/ExclWWHadronic2017Analysis/WWhadronic_QCDPt300to470_Pythia8_merge_12Apr2018_ntuplesv2.root/demo/ntp1");
      chain->Add("/eos/cms/store/user/jjhollar/ExclWWHadronic2017Analysis/WWhadronic_QCDPt470to600_Pythia8_merge_12Apr2018_ntuplesv2.root/demo/ntp1");
      chain->Add("/eos/cms/store/user/jjhollar/ExclWWHadronic2017Analysis/WWhadronic_QCDPt600to800_Pythia8_merge_12Apr2018_ntuplesv2.root/demo/ntp1");
      chain->Add("/eos/cms/store/user/jjhollar/ExclWWHadronic2017Analysis/WWhadronic_QCDPt800to1000_Pythia8_merge_12Apr2018_ntuplesv2.root/demo/ntp1");
      chain->Add("/eos/cms/store/user/jjhollar/ExclWWHadronic2017Analysis/WWhadronic_QCDPt1000to1400_Pythia8_merge_12Apr2018_ntuplesv2.root/demo/ntp1");
#endif
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

FullHadrAnalyzer::~FullHadrAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FullHadrAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FullHadrAnalyzer::LoadTree(Long64_t entry)
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

void FullHadrAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jet_pt = 0;
   jet_energy = 0;
   jet_phi = 0;
   jet_eta = 0;
   jet_mass = 0;
   jet_tau1 = 0;
   jet_tau2 = 0;
   jet_corrmass = 0;
   jet_vertexz = 0;
   dijet_mass = 0;
   dijet_pt = 0;
   dijet_y = 0;
   dijet_phi = 0;
   dijet_dphi = 0;
   pps_track_x = 0;
   pps_track_y = 0;
   pps_track_rpid = 0;
   hlt = 0;

#ifndef DAT
   gen_proton_px = 0;
   gen_proton_py = 0;
   gen_proton_pz = 0;
   gen_proton_xi = 0;
   gen_proton_t = 0;
   gen_jet_pt = 0;
   gen_jet_eta = 0;
   gen_jet_phi = 0;
   gen_jet_energy = 0;
   gen_dijet_mass = 0;
   gen_dijet_y = 0;
   jet_jer_res = 0;
   jet_jer_sf = 0;
   jet_jer_sfup = 0;
   jet_jer_sfdown = 0;
#endif

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_energy", &jet_energy, &b_jet_energy);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_mass", &jet_mass, &b_jet_mass);
   fChain->SetBranchAddress("jet_tau1", &jet_tau1, &b_jet_tau1);
   fChain->SetBranchAddress("jet_tau2", &jet_tau2, &b_jet_tau2);
   fChain->SetBranchAddress("jet_corrmass", &jet_corrmass, &b_jet_corrmass);
   fChain->SetBranchAddress("jet_vertexz", &jet_vertexz, &b_jet_vertexz);
   fChain->SetBranchAddress("dijet_mass", &dijet_mass, &b_dijet_mass);
   fChain->SetBranchAddress("dijet_pt", &dijet_pt, &b_dijet_pt);
   fChain->SetBranchAddress("dijet_y", &dijet_y, &b_dijet_y);
   fChain->SetBranchAddress("dijet_phi", &dijet_phi, &b_dijet_phi);
   fChain->SetBranchAddress("dijet_dphi", &dijet_dphi, &b_dijet_dphi);
   fChain->SetBranchAddress("pps_track_x", &pps_track_x, &b_pps_track_x);
   fChain->SetBranchAddress("pps_track_y", &pps_track_y, &b_pps_track_y);
   fChain->SetBranchAddress("pps_track_rpid", &pps_track_rpid, &b_pps_track_rpid);
   fChain->SetBranchAddress("hlt", &hlt, &b_hlt);
#ifndef DAT
   fChain->SetBranchAddress("gen_proton_px", &gen_proton_px, &b_gen_proton_px);
   fChain->SetBranchAddress("gen_proton_py", &gen_proton_py, &b_gen_proton_py);
   fChain->SetBranchAddress("gen_proton_pz", &gen_proton_pz, &b_gen_proton_pz);
   fChain->SetBranchAddress("gen_proton_xi", &gen_proton_xi, &b_gen_proton_xi);
   fChain->SetBranchAddress("gen_proton_t", &gen_proton_t, &b_gen_proton_t);
   fChain->SetBranchAddress("gen_jet_pt", &gen_jet_pt, &b_gen_jet_pt);
   fChain->SetBranchAddress("gen_jet_eta", &gen_jet_eta, &b_gen_jet_eta);
   fChain->SetBranchAddress("gen_jet_phi", &gen_jet_phi, &b_gen_jet_phi);
   fChain->SetBranchAddress("gen_jet_energy", &gen_jet_energy, &b_gen_jet_energy);
   fChain->SetBranchAddress("gen_dijet_mass", &gen_dijet_mass, &b_gen_dijet_mass);
   fChain->SetBranchAddress("gen_dijet_y", &gen_dijet_y, &b_gen_dijet_y);
   fChain->SetBranchAddress("jet_jer_res", &jet_jer_res, &b_jet_jer_res);
   fChain->SetBranchAddress("jet_jer_sf", &jet_jer_sf, &b_jet_jer_sf);
   fChain->SetBranchAddress("jet_jer_sfup", &jet_jer_sfup, &b_jet_jer_sfup);
   fChain->SetBranchAddress("jet_jer_sfdown", &jet_jer_sfdown, &b_jet_jer_sfdown);
   fChain->SetBranchAddress("gen_n_e", &gen_n_e, &b_gen_n_e);
   fChain->SetBranchAddress("gen_n_mu", &gen_n_mu, &b_gen_n_mu);
   fChain->SetBranchAddress("gen_n_tau", &gen_n_tau, &b_gen_n_tau);

   fChain->SetBranchAddress("mc_pu_trueinteractions_", &mc_pu_trueinteractions_, &b_mc_pu_trueinteractions);
#endif
   fChain->SetBranchAddress("pileupWeight", &pileupWeight, &b_pileupWeight);
   fChain->SetBranchAddress("nVertices", &nVertices, &b_nVertices);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   Notify();
}

Bool_t FullHadrAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FullHadrAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FullHadrAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FullHadrAnalyzer_cxx
