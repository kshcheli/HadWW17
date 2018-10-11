#define FullHadrAnalyzer_cxx

#define SIG
//#define PYT
//#define DAT

#include "FullHadrAnalyzer.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TH1.h>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSpline.h>
#include <TGraph.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TF1.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>



void FullHadrAnalyzer::Loop()
{

double w=1; TString fi="";
TString fnm;fnm="out/TestOutput";
#ifdef DAT
  cout << " ~~~~~   Running on the data ~~~~~~~~" << endl;
  fnm+="_dat";
  w=1.;
#endif

#ifdef SIG
  cout << " ~~~~~   Running on the signal MC ~~~~~~~~" << endl;
  fnm+="_sig";
  w=1.;
#endif

#ifdef PYT
  Long64_t NSKIM[]={156218, 14999535, 24828195,65529198,33670567, 19630824}; // number of skimmed events per each pt bin
  cout << " ~~~~~   Running on the QCD MC ~~~~~~~~" << endl;
  fnm+="_pyt";
  
  int NFILES=6;
   double sigma[6] = {103500.0, 6830.0 , 642.1 , 185.9 , 32.293 , 9.4183 }; // pb 
   int N[6] = {29829920, 65665222, 27384428, 65665222, 33677480, 19631814 }; // Total number of events per pt bin sample
   double scale[6] = {0., 0., 0., 0., 0., 0.};

  for(int i = 0; i < NFILES; i++){
   scale[i] = sigma[i]/N[i];
   cout <<  "Pythia xsection scale "<< scale[i] << endl;
  }
#endif
// ____ UTILITY VARS ___________

TRandom3 randomSrc; 
int nis[6];
Long64_t jump[]={0, 156218+1, 156218+14999535+1, 156218+14999535+24828195+1, 156218+14999535+24828195+65529198+1, 156218+14999535+24828195+65529198+33670567+1}; // for navigation 

//________________________HIST BOOKING ___________________

 TH1::SetDefaultSumw2();
 TH2::SetDefaultSumw2();
   map<string,TH1D*> h1;
   map<string,TH2D*> h2;
// TH1D

 h1["ptb"] = new TH1D("1", "pt_balance" , 100 , 0. , 5.);
 h1["ptb_my"] = new TH1D("2", "pt_balance my" , 100 , 0. , 5.);
 h1["ptb_bc"] = new TH1D("3", "pt_balance bc" , 100 , 0. , 5.);

 h1["tau21DDT1"] = new TH1D("4", "tau21DDT1" , 100 , 0. , 2.);
 h1["tau21DDT2"] = new TH1D("5", "tau21DDT2" , 100 , 0. , 2.);

 h1["pt1"] = new TH1D("6", "pt1" , 500 , 0. , 2000.);
 h1["pt2"] = new TH1D("7", "pt2" , 500 , 0. , 2000.);

 h1["m1"] = new TH1D("8", "m1" , 50 , 0. , 500.);
 h1["m2"] = new TH1D("9", "m2" , 50 , 0. , 500.);

 h1["eta1"] = new TH1D("10", "eta1" ,100., -3 , 3. );
 h1["eta2"] = new TH1D("11", "eta2" ,100., -3 , 3. );

 h1["e1"] = new TH1D("12", "e1" , 500 , 0. , 2000.);
 h1["e2"] = new TH1D("13", "e2" , 500 , 0. , 2000.);

 h1["phi1"] = new TH1D("14", "phi1" ,100, -4 , 4. );
 h1["phi2"] = new TH1D("15", "phi2" ,100, -4 , 4. );

 h1["prm1"] = new TH1D("16", "prm1" , 100. , 0. , 200.);
 h1["prm2"] = new TH1D("17", "prm2" , 100. , 0. , 200.);


 h1["dphi"] = new TH1D("18", "dphi" , 100. , 0. , 1.);
 h1["mdj"] = new TH1D("19", "mdj" , 1000. , 0. , 5000.);

 h1["mdjJER"] = new TH1D("20", "mdj JER" , 400. , 700. , 2100.);
 h1["mdjNOJER"] = new TH1D("21", "mdj NO JER" , 400. , 700. , 2100.);


 h1["pileupWeight"] = new TH1D("22", "puWeight" , 1000 , -2. , 2.);

// TH2D
 h2["test2"] = new TH2D("test", "test" , 100 , 0. , 30., 100, 0, 30);
// ______________________ LOOP START_____________________
 if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries /* frac*/;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


 if (jentry%100000 == 0){ cout <<  jentry << " done out of "<< nentries << endl; }
//___________________  SUPPRESS (SEMI-)LEPTONIC IF SIGNAL ______________
#ifdef SIG
if((gen_n_e+gen_n_mu+gen_n_tau)>0)continue;
w=pileupWeight;


cout << "Pileup weight: "<< pileupWeight << endl;
#endif

//___________________ QCD WEIGHTING ____________________

#ifdef PYT
int is=-1;


 TString current_file_name;
 current_file_name = (TString)((fChain->GetCurrentFile())->GetName());

//cout << current_file_name << endl;

    if      (current_file_name.Contains("170to300") ) { fi="_0"; w = scale[0];is=0; }
    else if (current_file_name.Contains("300to470") ) { fi="_1"; w = scale[1];is=1; }
    else if (current_file_name.Contains("470to600") ) { fi="_2"; w = scale[2];is=2; }
    else if (current_file_name.Contains("600to800") ) { fi="_3"; w = scale[3];is=3; }
    else if (current_file_name.Contains("800to1000") ) { fi="_4"; w = scale[4];is=4; }
    else if (current_file_name.Contains("1000to1400") ) { fi="_5"; w = scale[5];is=5; }
    else {cout << "Weight for Pythia missing!" << endl;break;}

w*=pileupWeight;
fi="";

if(jentry >= jump[is]+((int)(NSKIM[is]*0.1)) ) // use 10% of stats from each sample
 {
            if(is!=5){jentry=jump[is+1];} 
            if(is==5){break;} 
 }

#endif

//_________________ JER __________________________
double C_JER=1;
int a=0;int b=1;
#ifndef DAT
   TLorentzVector recojtmp;
   TLorentzVector gjtmp;
bool matched=false;
int indmatched=-1;

  double first =-999;int second = -999.; 
  int firsti =-999;int secondi = -999.; 
  

 for(int ij=0;ij<jet_pt->size();ij++ ){

   matched=false;
   indmatched=-1;
   C_JER=1;
   recojtmp.SetPtEtaPhiE((*jet_pt)[ij],(*jet_eta)[ij],(*jet_phi)[ij],(*jet_energy)[ij]); 

// Look for matching gen-level jet for given reco jet
   for(int ig=0; ig<gen_jet_pt->size(); ig++){ // loop over gen jets
      gjtmp.SetPtEtaPhiE((*gen_jet_pt)[ig],(*gen_jet_eta)[ig],(*gen_jet_phi)[ig],(*gen_jet_energy)[ig]); 
      if( (recojtmp.DeltaR(gjtmp) < (0.8/2.)) && (fabs(recojtmp.Pt() - gjtmp.Pt())<(3*(*jet_jer_res)[ij]*recojtmp.Pt()) ) ){matched=true; indmatched=ig;}// 0.8 is cone radius
    }

   if(matched){
              C_JER = 1 + ((*jet_jer_sf)[ij] - 1 )*( (recojtmp.Pt() - (*gen_jet_pt)[indmatched]) / recojtmp.Pt() );
             }
   else      {
               C_JER = 1 + randomSrc.Gaus(0, (*jet_jer_res)[ij]) * (sqrt( max( (*jet_jer_sf)[ij]*(*jet_jer_sf)[ij] - 1., 0.)   ) ) ;
             }

//cout << "C_JER "<< C_JER << endl;
//cout << "before " << (*jet_pt)[ij]<< endl;
if(C_JER<0) C_JER=0;
//  if(C_JER>0.){
  (*jet_pt)[ij]=(*jet_pt)[ij]*C_JER;
  (*jet_energy)[ij]=(*jet_energy)[ij]*C_JER;
 //}


     if ( (*jet_pt)[ij] > first) 
        { 
            second = first; secondi=firsti;
            first = (*jet_pt)[ij]; firsti=ij;
        } 
    else if ((*jet_pt)[ij] > second) 
        {           
            second = (*jet_pt)[ij]; secondi=ij;
        } 

//cout << "after " << (*jet_pt)[ij]<< endl;
  }
//if(secondi<firsti)cout << " ---------------------------> "<<endl;
//cout << "Leading subleading: "<< firsti << "  "<< secondi << endl;
a=firsti; b=secondi;
#endif

//_______________________    VARS ____________________________________
double tau21_DDT_a = 0; 
double tau21_DDT_b = 0;
TLorentzVector lj;
TLorentzVector slj; 
TLorentzVector dj;

//_________________________  CUTS ____________________________
bool NJCUT = false;

bool DETACUT = false;
bool ETACUTS = false;
bool PTCUTS = false;
bool DIJETMCUT = false;


bool B2GCUTS = false;

bool WTAGBOTH = false;
bool WTAGBOTHMY = false;
bool PRUNEDTAG = false;

bool LEADPRUNEDTAG = false;
bool SLEADPRUNEDTAG = false;


NJCUT=(jet_pt->size()>1);

if(NJCUT){

  lj.SetPtEtaPhiE((*jet_pt)[a],	(*jet_eta)[a],(*jet_phi)[a],(*jet_energy)[a]);
                    double rho_sh_a = log((lj.M()*lj.M())/lj.Pt()/1.) ; // 1GeV -- mu
                    tau21_DDT_a =  ((*jet_tau2)[a]/(*jet_tau1)[a]) - (-0.082)*rho_sh_a;


  slj.SetPtEtaPhiE((*jet_pt)[b],(*jet_eta)[b],(*jet_phi)[b],(*jet_energy)[b]);
                    double rho_sh_b = log ((slj.M()*slj.M())/slj.Pt()/1.) ; // 1GeV -- mu
                    tau21_DDT_b = ((*jet_tau2)[b]/(*jet_tau1)[b]) - (-0.082)*rho_sh_b;


dj = lj+slj;
DIJETMCUT = (dj.M()>1126.);
PTCUTS    = (lj.Pt()>200. && slj.Pt()>200.);
ETACUTS   = (fabs(lj.Eta())<2.5 && fabs(slj.Eta()) < 2.5);
DETACUT  = ( fabs(lj.Eta()-slj.Eta()) < 1.3 );

B2GCUTS=(NJCUT && PTCUTS && DIJETMCUT && ETACUTS && DETACUT);


WTAGBOTH = ( (tau21_DDT_a)<0.75 && (tau21_DDT_b)<0.75 && (*jet_corrmass)[a]>65. && (*jet_corrmass)[a]<105. && (*jet_corrmass)[b]>65. && (*jet_corrmass)[b]<105. );
WTAGBOTHMY = ( (tau21_DDT_a)<0.65 && (tau21_DDT_b)<0.65 && (*jet_corrmass)[a]>68. && (*jet_corrmass)[a]<88. && (*jet_corrmass)[b]>68. && (*jet_corrmass)[b]<88. );
PRUNEDTAG = ((*jet_corrmass)[a]>65. && (*jet_corrmass)[a]<105. && (*jet_corrmass)[b]>65. && (*jet_corrmass)[b]<105.);

LEADPRUNEDTAG = ( (*jet_corrmass)[a]>55. && (*jet_corrmass)[a]<215. );
SLEADPRUNEDTAG = ( (*jet_corrmass)[b]>55. && (*jet_corrmass)[b]<215. );

}
//________________________   HISTOGRAMS ETC ___________________________ 


    
     
if(B2GCUTS){

   h1["ptb_bc"]->Fill( (*jet_pt)[a] / (*jet_pt)[b], w);

   h1["dphi"]->Fill(1 - lj.DeltaPhi(slj)/3.14159, w);
   h1["mdj"]->Fill(dj.M(), w);

   h1["mdjJER"]->Fill(dj.M(), w);
   h1["mdjNOJER"]->Fill( (*dijet_mass)[0], w);

   h1["m1"]->Fill(lj.M(), w);
   h1["m2"]->Fill(slj.M(), w);

   h1["pt1"]->Fill(lj.Pt(), w);
   h1["pt2"]->Fill(slj.Pt(), w);

   h1["eta1"]->Fill(lj.Eta(), w);
   h1["eta2"]->Fill(slj.Eta(), w);

   h1["phi1"]->Fill(lj.Phi(), w);
   h1["phi2"]->Fill(slj.Phi(), w);

   h1["e1"]->Fill(lj.E(), w);
   h1["e2"]->Fill(slj.E(), w);


   h1["prm1"]->Fill((*jet_corrmass)[a], w);
   h1["prm2"]->Fill((*jet_corrmass)[b], w);


cout << "W inside loop "<< w<<endl;
 h1["pileupWeight"]->Fill(w); 
   if(WTAGBOTH && (1 - lj.DeltaPhi(slj)/3.14159)<0.1)h1["ptb"]->Fill( (*jet_pt)[a] / (*jet_pt)[b], w);
   if(WTAGBOTHMY)h1["ptb_my"]->Fill( (*jet_pt)[a] / (*jet_pt)[b], w);

   if(LEADPRUNEDTAG)h1["tau21DDT1"]->Fill(tau21_DDT_a, w);
   if(SLEADPRUNEDTAG)h1["tau21DDT2"]->Fill(tau21_DDT_b, w);


      }//b2gcuts
 }// event loop





TFile* f = new TFile(fnm+fi+"xchecks.root", "RECREATE");
for(map<string,TH1D*>::iterator it_histo = h1.begin(); it_histo != h1.end(); ++it_histo)
     (*it_histo).second->Write();
for(map<string,TH2D*>::iterator it_histo = h2.begin(); it_histo != h2.end(); ++it_histo)
     (*it_histo).second->Write();
f->Close();

}// Loop end





int run()
 {

 FullHadrAnalyzer m;
 m.Loop();
        
 return 0;
 }
   