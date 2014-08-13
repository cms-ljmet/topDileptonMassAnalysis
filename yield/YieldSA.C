#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
//#include "Objects/interface/MVAtree.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <cmath>

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooHistPdf.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"

const double PR_MU = 0.860;
const double FR_MU = 0.188;

const double PR_EL = 0.859;
const double FR_EL = 0.157;


float npSF[3] = {1.,1.,1.};
float dyEstimate[3] = {-1.,-1.,-1.};
float dySF[3] = {1.,1.,1.};
float mcSF = 1.;
int sign;

bool tex = true;

using namespace std;

typedef pair<int,float> IFPair;
typedef pair<float,float> FFPair;

enum Histos_t {AMWT, HT, MET, MINMLB, ST, NBJET, JETPT, RCN, CHHAD, NEHAD, HADDIFF, NJET, MLL};


	void Print4Vector(TLorentzVector const &v1){
	  printf("| %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f |", v1.Pt(), v1.Eta(), v1.Phi(),  v1.E() ,v1.M(), v1.Px() , v1.Py() , v1.Pz());
	}


class MVAtree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
  int run, lumi, event;
   Float_t        weight;
   Float_t        met;
   Int_t 	   nBTag, nJets, nLeptons;
   Float_t        lepton1_eta, lepton2_eta, lepton3_eta, lepton4_eta;
   Float_t        lepton1_pt, lepton2_pt, lepton3_pt, lepton4_pt;
   Float_t        lepton1_phi, lepton2_phi, lepton3_phi, lepton4_phi;
   Float_t        lepton1_energy, lepton2_energy, lepton3_energy, lepton4_energy;
   Float_t        jet1_pt, jet2_pt, jet3_pt, jet4_pt,jet5_pt, jet6_pt, jet7_pt, jet8_pt;
   Float_t        jet1_eta, jet2_eta, jet3_eta, jet4_eta,jet5_eta, jet6_eta, jet7_eta, jet8_eta;
   Float_t        jet1_phi, jet2_phi, jet3_phi, jet4_phi,jet5_phi, jet6_phi, jet7_phi, jet8_phi;
   Float_t        jet1_energy, jet2_energy, jet3_energy, jet4_energy,jet5_energy, jet6_energy, jet7_energy, jet8_energy;
   Float_t         isolation1, isolation2, isolation3, isolation4;
   Int_t           eventType;
   Int_t           jets_btag;
   Float_t         jet1_RCN, jet2_RCN, jet3_RCN, jet5_RCN, jet6_RCN, jet7_RCN, jet4_RCN, jet8_RCN;
   Int_t           jet1_nChHad, jet1_nNeHad, jet2_nChHad, jet2_nNeHad, jet3_nChHad, jet3_nNeHad, jet4_nChHad, jet4_nNeHad, jet5_nChHad, jet5_nNeHad, jet6_nChHad, jet6_nNeHad, jet7_nChHad, jet7_nNeHad, jet8_nChHad, jet8_nNeHad;

   Float_t        minDRJetLepton, minDRBJetLepton;
   Float_t        minMLB, maxMLB;
   Float_t        ht;
   Float_t        leptonSumMass;
   Float_t        CATopJetsBtag_disc1, CATopJetsBtag_disc2;
   Float_t        CATopJets1_pt, CATopJets2_pt;
   Float_t        CATopJets1_eta, CATopJets2_eta;
   Float_t        CATopJets1_phi, CATopJets2_phi;
   Float_t        CATopJets1_energy, CATopJets2_energy;
   Float_t        CAWBtag_disc1, CAWBtag_disc2;
   Float_t        CAWJets1_pt, CAWJets2_pt;
   Float_t        CAWJets1_eta, CAWJets2_eta;
   Float_t        CAWJets1_phi, CAWJets2_phi;
   Float_t        CAWJets1_energy, CAWJets2_energy;
   Int_t           nCAWJets, nCATopJets;
   Float_t         amwt;
   Float_t         amwtPeakWeight;
   Float_t        leptonSum, jetSum, leptonMETSum, leptonJetsSum, jetsMETSum, leptonJetsMETSum;
   Int_t           isMuon1, isMuon2, isMuon3, isMuon4;

   // List of branches
   TBranch        *b_weight, *b_PF_met_pt, *b_nBTag, *b_nJets, *b_nLeptons, *b_lepton1_eta, *b_lepton2_eta, *b_lepton3_eta, *b_lepton4_eta, *b_lepton1_pt, *b_lepton2_pt, *b_lepton3_pt, *b_lepton4_pt, *b_lepton1_phi, *b_lepton2_phi, *b_lepton3_phi, *b_lepton4_phi, *b_lepton1_energy, *b_lepton2_energy, *b_lepton3_energy, *b_lepton4_energy, *b_jet1_pt, *b_jet2_pt, *b_jet3_pt, *b_jet4_pt, *b_jet1_eta, *b_jet2_eta, *b_jet3_eta, *b_jet4_eta, *b_jet1_phi, *b_jet2_phi, *b_jet3_phi, *b_jet4_phi, *b_jet1_energy, *b_jet2_energy, *b_jet3_energy, *b_jet4_energy, *b_minDRJetLepton, *b_minDRBJetLepton, *b_minMLB, *b_HT, *b_leptonSumMass, *b_CATopJetsBtag_disc1, *b_CATopJetsBtag_disc2, *b_CATopJets1_pt, *b_CATopJets2_pt, *b_CATopJets1_eta, *b_CATopJets2_eta, *b_CATopJets1_phi, *b_CATopJets2_phi, *b_CATopJets1_energy, *b_CATopJets2_energy, *b_CAWBtag_disc1, *b_CAWBtag_disc2, *b_CAWJets1_pt, *b_CAWJets2_pt, *b_CAWJets1_eta, *b_CAWJets2_eta, *b_CAWJets1_phi, *b_CAWJets2_phi, *b_CAWJets1_energy, *b_CAWJets2_energy, *b_nCAWJets, *b_nCATopJets, *b_amwt, *b_amwtPeakWeight, *b_leptonSum, *b_jetSum, *b_leptonMETSum, *b_leptonJetsSum, *b_jetsMETSum, *b_leptonJetsMETSum, *b_isMuon1, *b_isMuon2, *b_isMuon3, *b_isMuon4;   //!
   TBranch        *b_jet5_pt, *b_jet5_eta, *b_jet5_phi, *b_jet5_energy, *b_jet6_pt, *b_jet6_eta, *b_jet6_phi, *b_jet6_energy, *b_jet7_pt, *b_jet7_eta, *b_jet7_phi, *b_jet7_energy, *b_jet8_pt, *b_jet8_eta, *b_jet8_phi, *b_jet8_energy, *b_isolation1, *b_isolation2, *b_isolation3, *b_isolation4, *b_eventType, *b_maxMLB, *b_jets_btag,
      *b_run, *b_lumi,*b_event   ;   //!
   TBranch        *b_jet1_nChHad, *b_jet1_nNeHad, *b_jet1_RCN, *b_jet2_nChHad, *b_jet2_nNeHad, *b_jet2_RCN, *b_jet3_nChHad, *b_jet3_nNeHad, *b_jet3_RCN, *b_jet4_nChHad, *b_jet4_nNeHad, *b_jet4_RCN, *b_jet5_nChHad, *b_jet5_nNeHad, *b_jet5_RCN, *b_jet6_nChHad, *b_jet6_nNeHad, *b_jet6_RCN, *b_jet7_nChHad, *b_jet7_nNeHad, *b_jet7_RCN, *b_jet8_nChHad, *b_jet8_nNeHad, *b_jet8_RCN;   //!

   MVAtree(TTree *tree=0);
   MVAtree(const TString fileName);
   void setChannel(int ch){channel=ch;}
   void setBJets(int min=0, int max=100) {minBJets=min;maxBJets=max;}
   void setJets(int min=0, int max=100) {minJets=min;maxJets=max;}
   void setMETcut(float m, float max) {metCut=m;metCutMax=max;}
   void setHTcut(float m, float max) {htCut=m;htCutMax=max;}
   void setAMWTcut(float m) {amwtCut=m;}
   void setjetSumCut(float m) {jetSumCut=m;}
   void setleptonMETSumCut(float m) {leptonMETSumCut=m;}
   void setleptonJetsSumCut(float m) {leptonJetsSumCut=m;}
   void setjetsMETSumCut(float m) {jetsMETSumCut=m;}
   void setleptonJetsMETSumCut(float m) {leptonJetsMETSumCut=m;}
   void setleptonSumCut(float m) {leptonSumCut=m;}
   void setEventType(int e){et=e;}
   void setCAW(bool e){useCAW=e;}
   void setminMlbcut(float m){minMlbCut=m;}
   void setLeptonSumMass  (float min=0.0, float max=1000.0) {minleptonSumMass=min;maxleptonSumMass=max;}
   void setLeptonSumMassEM(float min=0.0, float max=1000.0) {minleptonSumMassEM=min;maxleptonSumMassEM=max;}
   void clearSF(){sf[0] = sf[1] = sf[2] = 1.;}
   void setSF(float * a){sf[0] = a[0];sf[1] = a[1]; sf[2] = a[2];}
   float amwtCut;

   virtual ~MVAtree();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual float     Loop() {return LoopDouble().second;}
   pair<int,float>     LoopDouble();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
//    vector<TH1F*> histos;

  TH1F *newPlot;
  void plots(int type, string name)
    {doPlots=true;plotType=type, sampleName=name;}
  void noPlots()
    {doPlots=false;}

  private:
    int channel, minBJets, minJets, maxJets, et, maxBJets;
    float amwtUpperCut, amwtLowerCut, metCut, htCut, minMlbCut, htCutMax;
    float jetSumCut, leptonSumCut, leptonMETSumCut, leptonJetsSumCut, jetsMETSumCut, leptonJetsMETSumCut;
    float minleptonSumMass, maxleptonSumMass, minleptonSumMassEM, maxleptonSumMassEM, metCutMax;
    bool useCAW;
    int countLepton(int pdgid);
    inline float dr(float eta1, float phi1, float eta2, float phi2)
     {return sqrt((eta1-eta2)*(eta1-eta2)+(phi1-phi2)*(phi1-phi2));}
    inline bool isTag(int i){ return ((jets_btag && (1<<i))>>i );}
    string sampleName;
    bool doPlots;
   int plotType;

   float sf[3];
};


MVAtree::MVAtree(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("final.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("final.root");
      }
      f->GetObject("MVA",tree);

   }
   Init(tree);
}

MVAtree::MVAtree(const TString fileName) : fChain(0)
{
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName);
   if (!f) f = new TFile(fileName);
   TTree* tree = (TTree*)gDirectory->Get("MVA");
   Init(tree);
}

MVAtree::~MVAtree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MVAtree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MVAtree::LoadTree(Long64_t entry)
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

void MVAtree::Init(TTree *tree)
{
  et = 3;
  channel = -1;
  minBJets = minJets = 0;
  maxBJets = maxJets = 100;
  htCut=0.;
  htCutMax = -1;
  minMlbCut=-15;
  minleptonSumMass = 10000;
  maxleptonSumMass = 0;
  minleptonSumMassEM = 10000;
  maxleptonSumMassEM = 0;
  amwtUpperCut = 2000.;
  amwtLowerCut = 0.;
  metCut=0.;
  metCutMax = 9999999.;
  jetSumCut = 0.;
  leptonMETSumCut = 0.;
  leptonJetsSumCut = 0.;
  jetsMETSumCut = 0.;
  leptonJetsMETSumCut = 0.;
  leptonSumCut=0.;
  useCAW = false;
  sf[0] = sf[1] = sf[2] = 1.;
  newPlot = 0;
  sampleName = "";
  doPlots=false;
  amwtCut=00.;
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

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("met", &met, &b_PF_met_pt);
   fChain->SetBranchAddress("nBTag", &nBTag, &b_nBTag);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("nLeptons", &nLeptons, &b_nLeptons);
   fChain->SetBranchAddress("lepton1_eta", &lepton1_eta, &b_lepton1_eta);
   fChain->SetBranchAddress("lepton2_eta", &lepton2_eta, &b_lepton2_eta);
   fChain->SetBranchAddress("lepton3_eta", &lepton3_eta, &b_lepton3_eta);
   fChain->SetBranchAddress("lepton4_eta", &lepton4_eta, &b_lepton4_eta);
   fChain->SetBranchAddress("lepton1_pt", &lepton1_pt, &b_lepton1_pt);
   fChain->SetBranchAddress("lepton2_pt", &lepton2_pt, &b_lepton2_pt);
   fChain->SetBranchAddress("lepton3_pt", &lepton3_pt, &b_lepton3_pt);
   fChain->SetBranchAddress("lepton4_pt", &lepton4_pt, &b_lepton4_pt);
   fChain->SetBranchAddress("lepton1_phi", &lepton1_phi, &b_lepton1_phi);
   fChain->SetBranchAddress("lepton2_phi", &lepton2_phi, &b_lepton2_phi);
   fChain->SetBranchAddress("lepton3_phi", &lepton3_phi, &b_lepton3_phi);
   fChain->SetBranchAddress("lepton4_phi", &lepton4_phi, &b_lepton4_phi);
   fChain->SetBranchAddress("lepton1_energy", &lepton1_energy, &b_lepton1_energy);
   fChain->SetBranchAddress("lepton2_energy", &lepton2_energy, &b_lepton2_energy);
   fChain->SetBranchAddress("lepton3_energy", &lepton3_energy, &b_lepton3_energy);
   fChain->SetBranchAddress("lepton4_energy", &lepton4_energy, &b_lepton4_energy);
   fChain->SetBranchAddress("jet1_pt", &jet1_pt, &b_jet1_pt);
   fChain->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
   fChain->SetBranchAddress("jet3_pt", &jet3_pt, &b_jet3_pt);
   fChain->SetBranchAddress("jet4_pt", &jet4_pt, &b_jet4_pt);
   fChain->SetBranchAddress("jet1_eta", &jet1_eta, &b_jet1_eta);
   fChain->SetBranchAddress("jet2_eta", &jet2_eta, &b_jet2_eta);
   fChain->SetBranchAddress("jet3_eta", &jet3_eta, &b_jet3_eta);
   fChain->SetBranchAddress("jet4_eta", &jet4_eta, &b_jet4_eta);
   fChain->SetBranchAddress("jet1_phi", &jet1_phi, &b_jet1_phi);
   fChain->SetBranchAddress("jet2_phi", &jet2_phi, &b_jet2_phi);
   fChain->SetBranchAddress("jet3_phi", &jet3_phi, &b_jet3_phi);
   fChain->SetBranchAddress("jet4_phi", &jet4_phi, &b_jet4_phi);
   fChain->SetBranchAddress("jet1_energy", &jet1_energy, &b_jet1_energy);
   fChain->SetBranchAddress("jet2_energy", &jet2_energy, &b_jet2_energy);
   fChain->SetBranchAddress("jet3_energy", &jet3_energy, &b_jet3_energy);
   fChain->SetBranchAddress("jet4_energy", &jet4_energy, &b_jet4_energy);
   fChain->SetBranchAddress("jet5_pt", &jet5_pt, &b_jet5_pt);
   fChain->SetBranchAddress("jet5_eta", &jet5_eta, &b_jet5_eta);
   fChain->SetBranchAddress("jet5_phi", &jet5_phi, &b_jet5_phi);
   fChain->SetBranchAddress("jet5_energy", &jet5_energy, &b_jet5_energy);
   fChain->SetBranchAddress("jet6_pt", &jet6_pt, &b_jet6_pt);
   fChain->SetBranchAddress("jet6_eta", &jet6_eta, &b_jet6_eta);
   fChain->SetBranchAddress("jet6_phi", &jet6_phi, &b_jet6_phi);
   fChain->SetBranchAddress("jet6_energy", &jet6_energy, &b_jet6_energy);
   fChain->SetBranchAddress("jet7_pt", &jet7_pt, &b_jet7_pt);
   fChain->SetBranchAddress("jet7_eta", &jet7_eta, &b_jet7_eta);
   fChain->SetBranchAddress("jet7_phi", &jet7_phi, &b_jet7_phi);
   fChain->SetBranchAddress("jet7_energy", &jet7_energy, &b_jet7_energy);
   fChain->SetBranchAddress("jet8_pt", &jet8_pt, &b_jet8_pt);
   fChain->SetBranchAddress("jet8_eta", &jet8_eta, &b_jet8_eta);
   fChain->SetBranchAddress("jet8_phi", &jet8_phi, &b_jet8_phi);
   fChain->SetBranchAddress("jet8_energy", &jet8_energy, &b_jet8_energy);
   fChain->SetBranchAddress("jets_btag", &jets_btag, &b_jets_btag);
   fChain->SetBranchAddress("jet1_nChHad", &jet1_nChHad, &b_jet1_nChHad);
   fChain->SetBranchAddress("jet1_nNeHad", &jet1_nNeHad, &b_jet1_nNeHad);
   fChain->SetBranchAddress("jet1_RCN", &jet1_RCN, &b_jet1_RCN);
   fChain->SetBranchAddress("jet2_nChHad", &jet2_nChHad, &b_jet2_nChHad);
   fChain->SetBranchAddress("jet2_nNeHad", &jet2_nNeHad, &b_jet2_nNeHad);
   fChain->SetBranchAddress("jet2_RCN", &jet2_RCN, &b_jet2_RCN);
   fChain->SetBranchAddress("jet3_nChHad", &jet3_nChHad, &b_jet3_nChHad);
   fChain->SetBranchAddress("jet3_nNeHad", &jet3_nNeHad, &b_jet3_nNeHad);
   fChain->SetBranchAddress("jet3_RCN", &jet3_RCN, &b_jet3_RCN);
   fChain->SetBranchAddress("jet4_nChHad", &jet4_nChHad, &b_jet4_nChHad);
   fChain->SetBranchAddress("jet4_nNeHad", &jet4_nNeHad, &b_jet4_nNeHad);
   fChain->SetBranchAddress("jet4_RCN", &jet4_RCN, &b_jet4_RCN);
   fChain->SetBranchAddress("jet5_nChHad", &jet5_nChHad, &b_jet5_nChHad);
   fChain->SetBranchAddress("jet5_nNeHad", &jet5_nNeHad, &b_jet5_nNeHad);
   fChain->SetBranchAddress("jet5_RCN", &jet5_RCN, &b_jet5_RCN);
   fChain->SetBranchAddress("jet6_nChHad", &jet6_nChHad, &b_jet6_nChHad);
   fChain->SetBranchAddress("jet6_nNeHad", &jet6_nNeHad, &b_jet6_nNeHad);
   fChain->SetBranchAddress("jet6_RCN", &jet6_RCN, &b_jet6_RCN);
   fChain->SetBranchAddress("jet7_nChHad", &jet7_nChHad, &b_jet7_nChHad);
   fChain->SetBranchAddress("jet7_nNeHad", &jet7_nNeHad, &b_jet7_nNeHad);
   fChain->SetBranchAddress("jet7_RCN", &jet7_RCN, &b_jet7_RCN);
   fChain->SetBranchAddress("jet8_nChHad", &jet8_nChHad, &b_jet8_nChHad);
   fChain->SetBranchAddress("jet8_nNeHad", &jet8_nNeHad, &b_jet8_nNeHad);
   fChain->SetBranchAddress("jet8_RCN", &jet8_RCN, &b_jet8_RCN);
   fChain->SetBranchAddress("minDRJetLepton", &minDRJetLepton, &b_minDRJetLepton);
   fChain->SetBranchAddress("minDRBJetLepton", &minDRBJetLepton, &b_minDRBJetLepton);
   fChain->SetBranchAddress("minMLB", &minMLB, &b_minMLB);
   fChain->SetBranchAddress("maxMLB", &maxMLB, &b_maxMLB);
   fChain->SetBranchAddress("HT", &ht, &b_HT);
   fChain->SetBranchAddress("leptonSumMass", &leptonSumMass, &b_leptonSumMass);
   fChain->SetBranchAddress("CATopJetsBtag_disc1", &CATopJetsBtag_disc1, &b_CATopJetsBtag_disc1);
   fChain->SetBranchAddress("CATopJetsBtag_disc2", &CATopJetsBtag_disc2, &b_CATopJetsBtag_disc2);
   fChain->SetBranchAddress("CATopJets1_pt", &CATopJets1_pt, &b_CATopJets1_pt);
   fChain->SetBranchAddress("CATopJets2_pt", &CATopJets2_pt, &b_CATopJets2_pt);
   fChain->SetBranchAddress("CATopJets1_eta", &CATopJets1_eta, &b_CATopJets1_eta);
   fChain->SetBranchAddress("CATopJets2_eta", &CATopJets2_eta, &b_CATopJets2_eta);
   fChain->SetBranchAddress("CATopJets1_phi", &CATopJets1_phi, &b_CATopJets1_phi);
   fChain->SetBranchAddress("CATopJets2_phi", &CATopJets2_phi, &b_CATopJets2_phi);
   fChain->SetBranchAddress("CATopJets1_energy", &CATopJets1_energy, &b_CATopJets1_energy);
   fChain->SetBranchAddress("CATopJets2_energy", &CATopJets2_energy, &b_CATopJets2_energy);
   fChain->SetBranchAddress("CAWBtag_disc1", &CAWBtag_disc1, &b_CAWBtag_disc1);
   fChain->SetBranchAddress("CAWBtag_disc2", &CAWBtag_disc2, &b_CAWBtag_disc2);
   fChain->SetBranchAddress("CAWJets1_pt", &CAWJets1_pt, &b_CAWJets1_pt);
   fChain->SetBranchAddress("CAWJets2_pt", &CAWJets2_pt, &b_CAWJets2_pt);
   fChain->SetBranchAddress("CAWJets1_eta", &CAWJets1_eta, &b_CAWJets1_eta);
   fChain->SetBranchAddress("CAWJets2_eta", &CAWJets2_eta, &b_CAWJets2_eta);
   fChain->SetBranchAddress("CAWJets1_phi", &CAWJets1_phi, &b_CAWJets1_phi);
   fChain->SetBranchAddress("CAWJets2_phi", &CAWJets2_phi, &b_CAWJets2_phi);
   fChain->SetBranchAddress("CAWJets1_energy", &CAWJets1_energy, &b_CAWJets1_energy);
   fChain->SetBranchAddress("CAWJets2_energy", &CAWJets2_energy, &b_CAWJets2_energy);
   fChain->SetBranchAddress("nCAWJets", &nCAWJets, &b_nCAWJets);
   fChain->SetBranchAddress("nCATopJets", &nCATopJets, &b_nCATopJets);
   fChain->SetBranchAddress("amwt", &amwt, &b_amwt);
   fChain->SetBranchAddress("amwtPeakWeight", &amwtPeakWeight, &b_amwtPeakWeight);
   fChain->SetBranchAddress("leptonSum", &leptonSum, &b_leptonSum);
   fChain->SetBranchAddress("jetSum", &jetSum, &b_jetSum);
   fChain->SetBranchAddress("leptonMETSum", &leptonMETSum, &b_leptonMETSum);
   fChain->SetBranchAddress("leptonJetsSum", &leptonJetsSum, &b_leptonJetsSum);
   fChain->SetBranchAddress("jetsMETSum", &jetsMETSum, &b_jetsMETSum);
   fChain->SetBranchAddress("leptonJetsMETSum", &leptonJetsMETSum, &b_leptonJetsMETSum);
   fChain->SetBranchAddress("isMuon1", &isMuon1, &b_isMuon1);
   fChain->SetBranchAddress("isMuon2", &isMuon2, &b_isMuon2);
   fChain->SetBranchAddress("isMuon3", &isMuon3, &b_isMuon3);
   fChain->SetBranchAddress("isMuon4", &isMuon4, &b_isMuon4);
   fChain->SetBranchAddress("isolation1", &isolation1, &b_isolation1);
   fChain->SetBranchAddress("isolation2", &isolation2, &b_isolation2);
   fChain->SetBranchAddress("isolation3", &isolation3, &b_isolation3);
   fChain->SetBranchAddress("isolation4", &isolation4, &b_isolation4);
   fChain->SetBranchAddress("eventType", &eventType, &b_eventType);
   Notify();
}

Bool_t MVAtree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MVAtree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

pair<int,float> MVAtree::LoopDouble( )
// 	float leptonSumCut, float htCut, float leptonMETSumCut, float leptonJetsSumCut,
// 	float jetsMETSumCut, float leptonJetsMETSumCut, float amwtCut)
{
   if (fChain == 0) {
     cout << "Chain does not exist\n";
     return IFPair(-100000,-100000);
   }
int verbose = 0;
   Long64_t nentries = fChain->GetEntriesFast();
   if (verbose) cout << "Entries : "<<nentries<<endl;
   if (verbose) cout << channel << " "<<minBJets<<endl;


   if (doPlots) {
     //if (newPlot!=0) delete newPlot;
     string chName;

     switch( channel ) {
       case 0: chName = "MuMu"; break;
       case 1: chName = "ElMu"; break;
       case 2: chName = "ElEl"; break;
       default: chName = "All" ; break;
     }

    char name[300];
        Double_t xbins[1000];

    switch( plotType ) {
      case AMWT: 
	sprintf(name,"%s_%s_%s","amwt", sampleName.c_str(),chName.c_str());
	newPlot = new TH1F(name, name, 500, 100., 600.);
        break;
      case HT:
	for (Int_t i=0;i<=32;i++) xbins[i]=i*25;
	xbins[33]=1000;
	xbins[34]=1200;
	sprintf(name,"%s_%s_%s","HT", sampleName.c_str(),chName.c_str());
//	newPlot = new TH1F(name, name, 40,  0.,  1200.);
	newPlot = new TH1F(name, name, 34,  xbins);
        break;
      case MET: 
	sprintf(name,"%s_%s_%s","MET", sampleName.c_str(),chName.c_str());
	newPlot = new TH1F(name, name,200,  0., 1000.);
        break;
      case MINMLB: // ../DrawOneHistoWithSyst OS_23J_3.root  Top All OS_23J HT 1 minMlb  1 "min(m_{lb})"

	for (Int_t i=0;i<=16;i++) xbins[i]=i*25;
//	xbins[17]=500;
	xbins[17]=600;
// 	xbins[19]=700;
// 	xbins[20]=800;
      
	sprintf(name,"%s_%s_%s","minMlb", sampleName.c_str(),chName.c_str());
// 	newPlot = new TH1F(name, name,32,  0.,  800.);
 	newPlot = new TH1F(name, name,17, xbins);
        break;
      case ST: //../DrawOneHistoWithSyst OS.root Top All OS_23J leptonJetsMETSum 1 "ST [GeV]"
	for (Int_t i=0;i<=28;i++) xbins[i]=i*50;
	xbins[29]=1400;
	xbins[30]=1600;
	sprintf(name,"%s_%s_%s","leptonJetsMETSum", sampleName.c_str(),chName.c_str());
// 	newPlot = new TH1F(name, name,40,  0., 2000.);
	newPlot = new TH1F(name, name,30,  xbins);
        break;
      case NBJET:
	sprintf(name,"%s_%s_%s","nBJets", sampleName.c_str(),chName.c_str());
	newPlot = new TH1F(name, name,11,  -1.5, 9.5);
        break;
      case JETPT:
	sprintf(name,"%s_%s_%s","jetPt", sampleName.c_str(),chName.c_str());
	newPlot = new TH1F(name, name,200,  0., 2000);
        break;
      case RCN:
	sprintf(name,"%s_%s_%s","AK5JetRCN", sampleName.c_str(),chName.c_str());
	newPlot = new TH1F(name, name,500,   0.,  100.);
        break;
      case CHHAD:
	sprintf(name,"%s_%s_%s","AK5JetnChHad", sampleName.c_str(),chName.c_str());
	newPlot = new TH1F(name, name,100,   0.,  100.);
        break;
      case NEHAD:
	sprintf(name,"%s_%s_%s","AK5JetnNeuHad", sampleName.c_str(),chName.c_str());
	newPlot = new TH1F(name, name,100,   0.,  100.);
        break;
      case HADDIFF:
	sprintf(name,"%s_%s_%s","AK5JetnHadDiff", sampleName.c_str(),chName.c_str());
	newPlot = new TH1F(name, name,201,   -100,  100.);
        break;
      case NJET:
	sprintf(name,"%s_%s_%s","nJets", sampleName.c_str(),chName.c_str());
	newPlot = new TH1F(name, name,15,   -0.5,  14.5);
      case MLL:
	sprintf(name,"%s_%s_%s","mll", sampleName.c_str(),chName.c_str());
	newPlot = new TH1F(name, name,100,0.,  500.);
        break;
   }


  }




   Long64_t nbytes = 0, nb = 0;
   float selected = 0.;
   int events=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

// Channel check:
    int thisChannel;
    if (countLepton(13)==0) thisChannel=2;
    else if (countLepton(11)==0) thisChannel=0;
    else thisChannel=1;
    if ((channel>-1)&&(channel<3)&& (channel!=thisChannel)) continue;

    if (et!=eventType) continue;

    vector<TLorentzVector> leptonLV;
    vector<int> leptonF;
    if (nLeptons>0) {TLorentzVector lv;lv.SetPxPyPzE(lepton1_pt*cos(lepton1_phi), lepton1_pt*sin(lepton1_phi), lepton1_pt*TMath::SinH(lepton1_eta), lepton1_energy); leptonLV.push_back(lv); leptonF.push_back(isMuon1);}
    if (nLeptons>1) {TLorentzVector lv;lv.SetPxPyPzE(lepton2_pt*cos(lepton2_phi), lepton2_pt*sin(lepton2_phi), lepton2_pt*TMath::SinH(lepton2_eta), lepton2_energy); leptonLV.push_back(lv); leptonF.push_back(isMuon2);}
    if (nLeptons>2) {TLorentzVector lv;lv.SetPxPyPzE(lepton3_pt*cos(lepton3_phi), lepton3_pt*sin(lepton3_phi), lepton3_pt*TMath::SinH(lepton3_eta), lepton3_energy); leptonLV.push_back(lv); leptonF.push_back(isMuon3);}
    if (nLeptons>3) {TLorentzVector lv;lv.SetPxPyPzE(lepton4_pt*cos(lepton4_phi), lepton4_pt*sin(lepton4_phi), lepton4_pt*TMath::SinH(lepton4_eta), lepton4_energy); leptonLV.push_back(lv); leptonF.push_back(isMuon4);}


    vector<TLorentzVector> jetLV;
    vector<int> jetCh, jetNe;
    vector<float> jetRCN;
    if (nJets>0) {TLorentzVector lv;lv.SetPxPyPzE(jet1_pt*cos(jet1_phi), jet1_pt*sin(jet1_phi), jet1_pt*TMath::SinH(jet1_eta), jet1_energy); jetLV.push_back(lv); jetNe.push_back(jet1_nNeHad);jetCh.push_back(jet1_nChHad);jetRCN.push_back(jet1_RCN);}
    if (nJets>1) {TLorentzVector lv;lv.SetPxPyPzE(jet2_pt*cos(jet2_phi), jet2_pt*sin(jet2_phi), jet2_pt*TMath::SinH(jet2_eta), jet2_energy); jetLV.push_back(lv); jetNe.push_back(jet2_nNeHad);jetCh.push_back(jet2_nChHad);jetRCN.push_back(jet2_RCN);}
    if (nJets>2) {TLorentzVector lv;lv.SetPxPyPzE(jet3_pt*cos(jet3_phi), jet3_pt*sin(jet3_phi), jet3_pt*TMath::SinH(jet3_eta), jet3_energy); jetLV.push_back(lv); jetNe.push_back(jet3_nNeHad);jetCh.push_back(jet3_nChHad);jetRCN.push_back(jet3_RCN);}
    if (nJets>3) {TLorentzVector lv;lv.SetPxPyPzE(jet4_pt*cos(jet4_phi), jet4_pt*sin(jet4_phi), jet4_pt*TMath::SinH(jet4_eta), jet4_energy); jetLV.push_back(lv); jetNe.push_back(jet4_nNeHad);jetCh.push_back(jet4_nChHad);jetRCN.push_back(jet4_RCN);}
    if (nJets>4) {TLorentzVector lv;lv.SetPxPyPzE(jet5_pt*cos(jet5_phi), jet5_pt*sin(jet5_phi), jet5_pt*TMath::SinH(jet5_eta), jet5_energy); jetLV.push_back(lv); jetNe.push_back(jet5_nNeHad);jetCh.push_back(jet5_nChHad);jetRCN.push_back(jet5_RCN);}
    if (nJets>5) {TLorentzVector lv;lv.SetPxPyPzE(jet6_pt*cos(jet6_phi), jet6_pt*sin(jet6_phi), jet6_pt*TMath::SinH(jet6_eta), jet6_energy); jetLV.push_back(lv); jetNe.push_back(jet6_nNeHad);jetCh.push_back(jet6_nChHad);jetRCN.push_back(jet6_RCN);}
    if (nJets>6) {TLorentzVector lv;lv.SetPxPyPzE(jet7_pt*cos(jet7_phi), jet7_pt*sin(jet7_phi), jet7_pt*TMath::SinH(jet7_eta), jet7_energy); jetLV.push_back(lv); jetNe.push_back(jet7_nNeHad);jetCh.push_back(jet7_nChHad);jetRCN.push_back(jet7_RCN);}
    if (nJets>7) {TLorentzVector lv;lv.SetPxPyPzE(jet8_pt*cos(jet8_phi), jet8_pt*sin(jet8_phi), jet8_pt*TMath::SinH(jet8_eta), jet8_energy); jetLV.push_back(lv); jetNe.push_back(jet8_nNeHad);jetCh.push_back(jet8_nChHad);jetRCN.push_back(jet8_RCN);}

// Cuts:

    int myJets = 0;
    for (int ij = 0; ij<min(8,nJets);++ij){

      if ((ij==0) &&(jet1_pt> jetSumCut)) ++myJets;
      if ((ij==1) &&(jet2_pt> jetSumCut)) ++myJets;
      if ((ij==2) &&(jet3_pt> jetSumCut)) ++myJets;
      if ((ij==3) &&(jet4_pt> jetSumCut)) ++myJets;
      if ((ij==4) &&(jet5_pt> jetSumCut)) ++myJets;
      if ((ij==5) &&(jet6_pt> jetSumCut)) ++myJets;
      if ((ij==6) &&(jet7_pt> jetSumCut)) ++myJets;
      if ((ij==7) &&(jet8_pt> jetSumCut)) ++myJets;
    }
//        cout << jet1_pt <<" "<<jet2_pt<<" "<<jetSumCut<<endl;
//     cout << myJets<<nJets<<endl;
    int myBJets = nBTag;
    double myHT = ht;
    bool useJet[8] = {false, false, false, false, false, false, false, false};
      TLorentzVector jet, caw;
//if (nCAWJets>0) continue;
    if (useCAW) {
      bool cutPtCAW = false;
      myJets=nJets;
      int rem=0;
      for (int iw = 0; iw<min(2,nCAWJets);++iw){
        if (iw==0) caw.SetPxPyPzE(CAWJets1_pt*cos(CAWJets1_phi), CAWJets1_pt*sin(CAWJets1_phi), CAWJets1_pt*TMath::SinH(CAWJets1_eta), CAWJets1_energy);
        if (iw==1) caw.SetPxPyPzE(CAWJets2_pt*cos(CAWJets2_phi), CAWJets2_pt*sin(CAWJets2_phi), CAWJets2_pt*TMath::SinH(CAWJets2_eta), CAWJets2_energy);
       //cout << "CAW "<< iw<<" "<< caw.Pt()<<" "<< CAWJets1_pt<<" "<<myHT<<endl;
        if (caw.Pt()>200.) cutPtCAW=true;
	myHT+=caw.Pt();
	myJets+=2;
	for (int ij = 0; ij<min(8,nJets);++ij){
          if (ij==0) jet.SetPxPyPzE(jet1_pt*cos(jet1_phi), jet1_pt*sin(jet1_phi), jet1_pt*TMath::SinH(jet1_eta), jet1_energy);
          if (ij==1) jet.SetPxPyPzE(jet2_pt*cos(jet2_phi), jet2_pt*sin(jet2_phi), jet2_pt*TMath::SinH(jet2_eta), jet2_energy);
          if (ij==2) jet.SetPxPyPzE(jet3_pt*cos(jet3_phi), jet3_pt*sin(jet3_phi), jet3_pt*TMath::SinH(jet3_eta), jet3_energy);
          if (ij==3) jet.SetPxPyPzE(jet4_pt*cos(jet4_phi), jet4_pt*sin(jet4_phi), jet4_pt*TMath::SinH(jet4_eta), jet4_energy);
          if (ij==4) jet.SetPxPyPzE(jet5_pt*cos(jet5_phi), jet5_pt*sin(jet5_phi), jet5_pt*TMath::SinH(jet5_eta), jet5_energy);
          if (ij==5) jet.SetPxPyPzE(jet6_pt*cos(jet6_phi), jet6_pt*sin(jet6_phi), jet6_pt*TMath::SinH(jet6_eta), jet6_energy);
          if (ij==6) jet.SetPxPyPzE(jet7_pt*cos(jet7_phi), jet7_pt*sin(jet7_phi), jet7_pt*TMath::SinH(jet7_eta), jet7_energy);
          if (ij==7) jet.SetPxPyPzE(jet8_pt*cos(jet8_phi), jet8_pt*sin(jet8_phi), jet8_pt*TMath::SinH(jet8_eta), jet8_energy);
//       cout << ij<<" "<< iw<<" "<< jet.DeltaR(caw)<<" "<< myJets<<" "<<ht<<endl;
	  if (jet.DeltaR(caw) < 0.8) {
       //cout << "removed "<< iw<<" "<< caw.Pt()<<" "<< jet.Pt()<<" "<<myHT<<endl;
	    myHT-=jet.Pt();
	    --myJets;
	    ++rem;
	    continue;
	  }
	  useJet[ij]=true;
	}
      }
      if ((nCAWJets>0) && (!cutPtCAW)) continue;
      //cout << myJets<<" "<< nJets<<" "<< nCAWJets<<" "<< rem<<" "<< myHT<<" "<<ht<<endl;
      if (myJets<nCAWJets) exit(1);
    }

    if ((nBTag<minBJets)||(nBTag>maxBJets)) continue;
    if ((myJets<minJets)||(myJets>maxJets)) continue;
//     if ((myJets==2 && nBTag<1) || (myJets==3 && nBTag<2)) continue;
    if ( (sign==-1) && (thisChannel==0 || thisChannel==2) &&(met<metCut)) continue;
    if ( (sign==-1) && (thisChannel==0 || thisChannel==2) &&(met>metCutMax)) continue;
//     if (met<metCut) continue;
//     if (met>metCutMax) continue;

    if (myHT<htCut) continue;
    if ((htCutMax>2.)&&myHT>htCutMax) continue;
    if (minMLB<minMlbCut) continue;

//     if (jetSum	   <jetSumCut	       ) continue;
    if (leptonMETSum    <leptonMETSumCut    ) continue;
  if (verbose>1) cout << "c \n";

    if (leptonSum   <leptonSumCut   ) continue;
    if (leptonJetsSum   <leptonJetsSumCut   ) continue;
    if (jetsMETSum	   <jetsMETSumCut      ) continue;
    if (leptonJetsMETSum<leptonJetsMETSumCut) continue;
    if (verbose>1) cout << "d \n";
    if ( (sign==-1) && (thisChannel==0 || thisChannel==2) &&
    (leptonSumMass>minleptonSumMass)&&(leptonSumMass<maxleptonSumMass)) continue;
    if ( (sign==-1) && (thisChannel==1) &&
    (leptonSumMass>minleptonSumMassEM)&&(leptonSumMass<maxleptonSumMassEM)) continue;

//     if (amwtPeakWeight<25) continue;
    if (amwt<amwtCut) continue;

//      bool getJet = false;
//      for (unsigned int il = 0;il<jetLV.size();il++ ){
// 	if ((abs(jetLV[il].Eta())>0.9) && (abs(jetLV[il].Eta())<1.6) && (jetCh[il]-jetNe[il])>40) getJet=true;
// }if (getJet) continue;

      selected+=(weight*sf[thisChannel]);
      ++events;

   if (doPlots) {
     bool getJet = false;
     for (unsigned int il = 0;il<jetLV.size();il++ )
	      if ((abs(jetLV[il].Eta())>0.9) && (abs(jetLV[il].Eta())<1.6)) getJet = true;
     switch( plotType ) {
      case 0: if (amwtPeakWeight>0.) newPlot->Fill(amwt,(weight*sf[thisChannel])); break;
      case 1: newPlot->Fill(myHT,(weight*sf[thisChannel])); break;
      case 2: //if (!getJet) 
      newPlot->Fill(met,(weight*sf[thisChannel])); break;
      case 3: newPlot->Fill(minMLB,(weight*sf[thisChannel])); break;
      case 4: newPlot->Fill(leptonJetsMETSum,(weight*sf[thisChannel])); break;
      case NBJET: {
      newPlot->Fill(nBTag,(weight*sf[thisChannel]));
      newPlot->Fill(nBTag,(weight*sf[thisChannel])); break;}
      case JETPT:  	 for (unsigned int il = 0;il<jetLV.size();il++ ){
	      if ((abs(jetLV[il].Eta())>0.9) && (abs(jetLV[il].Eta())<1.6))
	      newPlot->Fill(jetLV[il].Pt(),(weight*sf[thisChannel]));} break;
      case RCN:for (unsigned int il = 0;il>jetLV.size();il++ ){
	      if (true)
	      newPlot->Fill(jetRCN[il],(weight*sf[thisChannel]));} break;
        break;
      case CHHAD:for (unsigned int il = 0;il>jetLV.size();il++ ){
	      if (true)
	      newPlot->Fill(jetCh[il],(weight*sf[thisChannel]));}
        break;
      case NEHAD:for (unsigned int il = 0;il>jetLV.size();il++ ){
	      if (true)
	      newPlot->Fill(jetNe[il],(weight*sf[thisChannel]));}
        break;
      case HADDIFF:for (unsigned int il = 0;il<jetLV.size();il++ ){
	      if (true)
	      newPlot->Fill(jetCh[il]-jetNe[il],(weight*sf[thisChannel]));}
        break;
      case NJET:newPlot->Fill(myJets,(weight*sf[thisChannel])); break; 
      case MLL:newPlot->Fill(leptonSumMass,(weight*sf[thisChannel])); break; 
     }

   }
   }

   return pair<int,float>(events,selected);
}

vector<MVAtree*> backgroundsTree;
typedef map< int, vector<MVAtree*> > TreeMap;
typedef map< int, vector<MVAtree*> >::iterator TreeMapIt;
TreeMap signalTreeMap;
MVAtree* signalTree;
vector<MVAtree*> dataTree;
MVAtree* nonPrompt;
vector<string> backgroundsName;
vector<double> backgroundsUnc;
map< int, vector<string> > signalName;
vector<string> dataName;
vector<float> backgroundsNbr;
vector<float> signalNbr;
map< int, vector<float> > signalBR, signalNorm;
map< string, float > initialNumber;
map< int, float > xSection, newXSection;
int numTags;

char *getTagString() {
    char* tag = new char[3];
    sprintf(tag,"_%iB",numTags);
    return tag;
}

void setSumCuts(
	float jetSumCut 	 ,
	float leptonMETSumCut	 ,
	float leptonJetsSumCut   ,
	float jetsMETSumCut	 ,
	float leptonJetsMETSumCut,
	float leptonSumCut, bool verbose = true) {
  if (verbose) {
    cout << "--> Require jetSumCut	   "<<jetSumCut 	 <<endl;
    cout << "--> Require leptonMETSumCut     "<<leptonMETSumCut	 <<endl;
    cout << "--> Require leptonJetsSumCut    "<<leptonJetsSumCut   <<endl;
    cout << "--> Require jetsMETSumCut	   "<<jetsMETSumCut	 <<endl;
    cout << "--> Require leptonJetsMETSumCut "<<leptonJetsMETSumCut<<endl;
    cout << "--> Require leptonSumCut 	   "<<leptonSumCut 	 <<endl;
  }

  for (TreeMapIt itS=signalTreeMap.begin();itS != signalTreeMap.end();++itS) {
    for (unsigned int i=0;i<itS->second.size();++i) {
      signalTree= itS->second[i];
      signalTree->setjetSumCut 	 (jetSumCut	     );
      signalTree->setleptonMETSumCut    (leptonMETSumCut    );
      signalTree->setleptonJetsSumCut   (leptonJetsSumCut   );
      signalTree->setjetsMETSumCut      (jetsMETSumCut	  );
      signalTree->setleptonJetsMETSumCut(leptonJetsMETSumCut);
      signalTree->setleptonSumCut       (leptonSumCut  );
    }
  }
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setjetSumCut 	         (jetSumCut	     );
    backgroundsTree[i]->setleptonMETSumCut	 (leptonMETSumCut    );
    backgroundsTree[i]->setleptonJetsSumCut	 (leptonJetsSumCut   );
    backgroundsTree[i]->setjetsMETSumCut	 (jetsMETSumCut      );
    backgroundsTree[i]->setleptonJetsMETSumCut   (leptonJetsMETSumCut);
    backgroundsTree[i]->setleptonSumCut (leptonSumCut);  }
  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->setjetSumCut 	 	 (jetSumCut	     );
    dataTree[i]->setleptonMETSumCut	 (leptonMETSumCut    );
    dataTree[i]->setleptonJetsSumCut	 (leptonJetsSumCut   );
    dataTree[i]->setjetsMETSumCut	 (jetsMETSumCut      );
    dataTree[i]->setleptonJetsMETSumCut  (leptonJetsMETSumCut);
    dataTree[i]->setleptonSumCut (leptonSumCut);  }
  if (nonPrompt!=0) {
    nonPrompt->setjetSumCut 	 	 (jetSumCut	     );
    nonPrompt->setleptonMETSumCut        (leptonMETSumCut    );
    nonPrompt->setleptonJetsSumCut       (leptonJetsSumCut   );
    nonPrompt->setjetsMETSumCut	 (jetsMETSumCut      );
    nonPrompt->setleptonJetsMETSumCut  (leptonJetsMETSumCut);
    nonPrompt->setleptonSumCut (leptonSumCut);
  }
}

double getTotFakesSS(double f1, double f2, double p1, double p2,
double nt00, double nt01,double nt10,double nt11){

  double det = 1.0 / ((p1 - f1)*(p2 - f2));
  double hNff=(nt00*p1*p2 - nt10*(1-p1)*p2 - nt01 *p1*(1-p2) + nt11*(1-p1)*(1-p2))*det*f1*f2;
  double hNpf=(-nt00*f1*p2 + nt10*(1-f1)*p2 + nt01 *f1*(1-p2) - nt11*(1-f1)*(1-p2))*det*p1*f2;
  double hNfp=(-nt00*p1*f2 + nt10*(1-p1)*f2 + nt01 *p1*(1-f2) - nt11*(1-p1)*(1-f2))*det*f1*p2;
  return hNff+hNpf+hNfp;
}
double getTotFakesSSErr(double f1, double f2, double p1, double p2,
double nt00, double nt01,double nt10,double nt11){

  double det = 1.0 / ((p1 - f1)*(p2 - f2));
  double a = nt00*pow((p1*p2*det*f1*f2 -f1*p2*det*p1*f2 -p1*f2*det*f1*p2),2) +
	     nt10*pow((-(1-p1)*p2*det*f1*f2 + (1-f1)*p2*det*p1*f2 + (1-p1)*f2*det*f1*p2) , 2) +
	     nt01*pow((p1*(1-p2)*det*f1*f2 + f1*(1-p2)*det*p1*f2 + p1*(1-f2)*det*f1*p2), 2) +
	     nt11*pow(((1-p1)*(1-p2)*det*f1*f2 - (1-f1)*(1-p2)*det*p1*f2 - (1-p1)*(1-f2)*det*f1*p2),2);
  return sqrt(a);
}

double getTotFakesTRI(double f, double p,
  double nt0, double nt1,double nt2,double nt3){

  double det = 1.0/((p-f)*(p-f)*(p-f));

  double hNfff = (nt0*p*p*p -p*p*(1-p)*nt1 +p*(1-p)*(1-p)*nt2 -(1-p)*(1-p)*(1-p)*nt3)*det;
  double hNffp = (-3*p*p*f*nt0 +(2*p*f*(1-p)+p*p*(1-f))*nt1
      -(f*(1-p)*(1-p) + 2*p*(1-p)*(1-f))*nt2 + 3*(1-p)*(1-p)*(1-f)*nt3)*det;
  double hNfpp = (3*p*f*f*nt0 -(f*f*(1-p) + 2*p*f*(1-f))*nt1
      + (2*f*(1-p)*(1-f) + p*(1-f)*(1-f))*nt2 -3*(1-p)*(1-f)*(1-f)*nt3)*det;
//   cout <<  "check fr "<<hNfff <<" "<< hNffp<<" "<<hNfpp<<" "<< hNfff + hNffp+hNfpp<<" "<<det<<endl;
// cout << nt0<<" "<< nt1<<" "<< nt2<<" "<< nt3<<endl;
  //return hNfff + hNffp+hNfpp;
  return hNfff*f*f*f + hNffp*f*f*p +hNfpp*f*p*p;
  
}
double getTotFakesTRIErr(double f, double p,
  double nt0, double nt1,double nt2,double nt3){

  double det = 1.0/((p-f)*(p-f)*(p-f));

  double a =	nt0* pow( ( p*p*p*f*f*f -3*p*p*f*f*f*p +3*p*f*f*f*p*p) ,2) +
  		nt1* pow( ( -p*p*(1-p)*f*f*f + (2*p*f*(1-p)+p*p*(1-f))*f*f*p -(f*f*(1-p) + 2*p*f*(1-f))*f*p*p) ,2)+
		nt2* pow( ( p*(1-p)*(1-p)*f*f*f -(f*(1-p)*(1-p)*f*f*p + 2*p*(1-p)*(1-f))+ (2*f*(1-p)*(1-f) + p*(1-f)*(1-f))*f*p*p) ,2) +
		nt3* pow( ( -(1-p)*(1-p)*(1-p)*f*f*f + 3*(1-p)*(1-p)*(1-f)*f*f*p-3*(1-p)*(1-f)*(1-f)*f*p*p) ,2);
//   cout <<  "check err "<<( p*p*p -3*p*p*f +3*p*f*f)*det <<" "<< 
//   pow( ( p*p*p -3*p*p*f* +3*p*f*f) ,2) <<" "<<  nt0* pow( ( p*p*p -3*p*p*f +3*p*f*f) ,2) <<" "<< 
//   sqrt(a)*det<<" "<<det<<endl;
//   cout <<  "check err "<< p*p*p*det <<" "<<-3*p*p*f*det <<" "<< 3*p*f*f*det<<endl;
//   cout <<  "check err "<< p*p*p*det -3*p*p*f*det + 3*p*f*f*det<<endl;
//   cout <<  "check err "<< p*p*p*det -3*p*p*f*det + 3*p*f*f*det<<endl;
//   cout <<  "check err "<< ( p*p*p -3*p*p*f +3*p*f*f)*det<<endl;
//   cout <<  "check err "<< ( p*p*p -3*p*p*f* +3*p*f*f)*det<<endl;
// cout << nt0<<" "<< nt1<<" "<< nt2<<" "<< nt3<<endl;
// cout <<  "\ncheck err "<< p <<" "<<f<<" "<< (nt0* ( p*p*p -3*p*p*f +3*p*f*f)  +
//   		nt1* ( -p*p*(1-p) + (2*p*f*(1-p)+p*p*(1-f)) -(f*f*(1-p) + 2*p*f*(1-f))) +
// 		nt2* ( p*(1-p)*(1-p) -(f*(1-p)*(1-p) + 2*p*(1-p)*(1-f))+ (2*f*(1-p)*(1-f) + p*(1-f)*(1-f)))  +
// 		nt3* ( -(1-p)*(1-p)*(1-p) + 3*(1-p)*(1-p)*(1-f)-3*(1-p)*(1-f)*(1-f)))*det <<" "<<
//  sqrt(a)*det<<endl;
  return sqrt(a)*det;
}

void setCAW(bool b, bool verbose = true) {
  if (verbose) cout << "--> Require CAW "<<b<<endl;
  for (TreeMapIt itS=signalTreeMap.begin();itS != signalTreeMap.end();++itS) {
    for (unsigned int i=0;i<itS->second.size();++i) {
      itS->second[i]->setCAW(b);;
    }
  }
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setCAW(b);
  }
  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->setCAW(b);
  }
  if (nonPrompt!=0) nonPrompt->setCAW(b);
}


void setJets(int min=0, int max=100, bool verbose = false) {
  if (verbose) cout << "--> Require between "<<min<< " and "<<max<< " jets \n";
  for (TreeMapIt itS=signalTreeMap.begin();itS != signalTreeMap.end();++itS) {
    for (unsigned int i=0;i<itS->second.size();++i) {
      itS->second[i]->setJets(min, max);
    }
  }
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setJets(min, max);
  }
  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->setJets(min, max);
  }
  if (nonPrompt!=0) nonPrompt->setJets(min, max);
}
void setBJets(int min=0, int max=100, bool verbose = true) {
  if (verbose) cout << "--> Require at between "<<min<< " and "<<max<< " b-tagged jets \n";
  for (TreeMapIt itS=signalTreeMap.begin();itS != signalTreeMap.end();++itS) {
    for (unsigned int i=0;i<itS->second.size();++i) {
      itS->second[i]->setBJets(min, max);
    }
  }
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setBJets(min, max);
  }
  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->setBJets(min, max);
  }
  if (nonPrompt!=0) nonPrompt->setBJets(min, max);
}

void setLeptonSumMass(float min=10000.0, float max=0.0, bool verbose = false) {
  if (verbose) cout << "--> Remove EE/MM events with lepton invairant mass between "<<min<< " and "<<max<< " GeV \n";
  for (TreeMapIt itS=signalTreeMap.begin();itS != signalTreeMap.end();++itS) {
    for (unsigned int i=0;i<itS->second.size();++i) {
      itS->second[i]->setLeptonSumMass(min, max);
    }
  }
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setLeptonSumMass(min, max);
  }
  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->setLeptonSumMass(min, max);
  }
  if (nonPrompt!=0) nonPrompt->setLeptonSumMass(min, max);
}
void setLeptonSumMassEM(float min=10000.0, float max=0.0, bool verbose = false) {
  if (verbose) cout << "--> Remove EM events with lepton invairant mass between "<<min<< " and "<<max<< " GeV \n";
  for (TreeMapIt itS=signalTreeMap.begin();itS != signalTreeMap.end();++itS) {
    for (unsigned int i=0;i<itS->second.size();++i) {
      itS->second[i]->setLeptonSumMassEM(min, max);
    }
  }
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setLeptonSumMassEM(min, max);
  }
  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->setLeptonSumMassEM(min, max);
  }
  if (nonPrompt!=0) nonPrompt->setLeptonSumMassEM(min, max);
}

void setMETcut(float min, float max, bool verbose) {
  if (verbose) cout << "--> Set MET cut between "<<min<< " and "<<max<< " GeV \n";
  for (TreeMapIt itS=signalTreeMap.begin();itS != signalTreeMap.end();++itS) {
    for (unsigned int i=0;i<itS->second.size();++i) {
      itS->second[i]->setMETcut(min, max);
    }
  }
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setMETcut(min, max);
  }
  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->setMETcut(min, max);
  }
  if (nonPrompt!=0) nonPrompt->setMETcut(min, max);
}

void setHTcut(float min, float max, bool verbose) {
  if (verbose) cout << "--> Set HT cut between "<<min<< " and "<<max<< " GeV \n";
  for (TreeMapIt itS=signalTreeMap.begin();itS != signalTreeMap.end();++itS) {
    for (unsigned int i=0;i<itS->second.size();++i) {
      itS->second[i]->setHTcut(min, max);
    }
  }
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setHTcut(min, max);
  }
  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->setHTcut(min, max);
  }
  if (nonPrompt!=0) nonPrompt->setHTcut(min, max);
}
void setAMWTcut(float b, bool verbose = true) {
  if (verbose) cout << "--> Set AMWT cut to "<<b<<endl;
  for (TreeMapIt itS=signalTreeMap.begin();itS != signalTreeMap.end();++itS) {
    for (unsigned int i=0;i<itS->second.size();++i) {
      itS->second[i]->setAMWTcut(b);;
    }
  }
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setAMWTcut(b);
  }
  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->setAMWTcut(b);
  }
  if (nonPrompt!=0) nonPrompt->setAMWTcut(b);
}

void setminMlbcut(float b, bool verbose = true) {
  if (verbose) cout << "--> Set minMlb cut to "<<b<<endl;
  for (TreeMapIt itS=signalTreeMap.begin();itS != signalTreeMap.end();++itS) {
    for (unsigned int i=0;i<itS->second.size();++i) {
      itS->second[i]->setminMlbcut(b);;
    }
  }
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setminMlbcut(b);
  }
  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->setminMlbcut(b);
  }
  if (nonPrompt!=0) nonPrompt->setminMlbcut(b);
}

FFPair nonPB(int channel) {
  nonPrompt->setChannel(channel);
  nonPrompt->setEventType(0);
  float nt00 = nonPrompt->Loop();
  nonPrompt->setEventType(1);
  float nt01 = nonPrompt->Loop();
  nonPrompt->setEventType(2);
  float nt10 = nonPrompt->Loop();
  nonPrompt->setEventType(3);
  float nt11 = nonPrompt->Loop();
//   cout << "nt "<< nt00<<" "<< nt01<<" "<< nt10<<" "<< nt11<<endl;
  if (sign==1) {
//   cout << getTotFakesSS(FR_MU, FR_MU, PR_MU, PR_MU, nt00, nt01, nt10, nt11) <<" "<<
// 	 getTotFakesSSErr(FR_MU, FR_MU, PR_MU, PR_MU, nt00, nt01, nt10, nt11) <<" "<<
//     getTotFakesSS(FR_EL, FR_MU, PR_EL, PR_MU, nt00, nt01, nt10, nt11)<<" "<<
//     getTotFakesSSErr(FR_EL, FR_MU, PR_EL, PR_MU, nt00, nt01, nt10, nt11)<<" "<<
//     getTotFakesSS(FR_EL, FR_EL, PR_EL, PR_EL,nt00, nt01, nt10, nt11)<<" "<<
//     getTotFakesSSErr(FR_EL, FR_EL, PR_EL, PR_EL,nt00, nt01, nt10, nt11)<<" "<<endl;
    float fra, frb, pra, prb;
    if (channel==0) {fra= frb = FR_MU; pra = prb = PR_MU;}
    else if (channel==1) {fra = FR_EL; frb = FR_MU; pra = PR_EL; prb = PR_MU;}
    else {fra= frb = FR_EL; pra = prb = PR_EL;}
    return FFPair(
    	getTotFakesSS(fra, frb, pra, prb, nt00, nt01, nt10, nt11),
    	getTotFakesSSErr(fra, frb, pra, prb, nt00, nt01, nt10, nt11) );
  } else if (sign==0) {
//   cout << getTotFakesTRI(FR_MU, PR_MU, nt00, nt01, nt10, nt11)  <<" "<<
// 	  getTotFakesTRIErr(FR_MU, PR_MU, nt00, nt01, nt10, nt11)  <<" "<<
//     	  getTotFakesTRI(FR_MU, PR_EL, nt00, nt01, nt10, nt11)  <<" "<<
//     	  getTotFakesTRIErr(FR_MU, PR_EL, nt00, nt01, nt10, nt11)  <<" "<<
//     	  getTotFakesTRI(FR_EL, PR_EL, nt00, nt01, nt10, nt11)  <<" "<<
//     	  getTotFakesTRIErr(FR_EL, PR_EL, nt00, nt01, nt10, nt11)  <<" "<<endl;

    float fr, pr;
    if (channel==0) {fr = FR_MU; pr = PR_MU;}
    else if (channel==1) {fr = FR_MU; pr = PR_EL;}
    else {fr = FR_EL; pr = PR_EL;}
    return FFPair(
	getTotFakesTRI(fr, pr, nt00, nt01, nt10, nt11),
	getTotFakesTRIErr(fr, pr, nt00, nt01, nt10, nt11) );
  }
  return FFPair(-1,-1);
}

void nonPB(float* totalFake, bool applySF = true) {
  for (int j=0;j<3;++j) {
    totalFake[j] = nonPB(j).first*(applySF?npSF[j]:1.0);
  }
}

void nonPB(pair<float,float>* totalFake, bool applySF = true) {
//   for (int j=0;j<3;++j) {
//     totalFake[j] = nonPB(j)*(applySF?npSF[j]:1.0);
//   }
}


float data(int channel = -1)
{
  dataTree[0]->setChannel(channel);
  return dataTree[0]->Loop();
}

float background(int channel = -1) {
  float totalBackground = 0.;
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setChannel(channel);
    totalBackground+=(backgroundsTree[i]->Loop()*mcSF*
    	(sign==-1&&i==(backgroundsTree.size()-1)&&(channel!=-1)?dySF[channel]:1.0));
  }
  return totalBackground;
}


void test(int mass, bool verbose = false) {
  float totalSignal = 0.;
  signalNbr.clear();
  backgroundsNbr.clear();
  float totalBackground = 0.;
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setChannel(-1);
    backgroundsNbr.push_back( backgroundsTree[i]->Loop()*mcSF);
    if (verbose) cout << backgroundsName[i] <<" "<< backgroundsNbr[i]<<endl;
    totalBackground+=backgroundsNbr[i];
  }

//   for (unsigned int i=0;i<signalTree.size();++i) {
//     signalNbr.push_back( signalTree[i]->Loop());
//     if (verbose) cout << signalName[i] <<" "<< signalNbr[i]<<"\t"<< signalNbr[i]/totalBackground<<"\t"<<signalNbr[i]/sqrt(signalNbr[i]+totalBackground)<<endl;
//   }
  if (nonPrompt!=0) {
    float totalFake[3];
    nonPB(totalFake);
    if (verbose) cout << "Non-Prompt\t"<<fixed<<setprecision(2)
	 << totalFake[0]+totalFake[1]+totalFake[2] << "\n";
    totalBackground+=totalFake[0]+totalFake[1]+totalFake[2];
  }

  for (unsigned int i=0;i<signalTreeMap[mass].size();++i) {
    signalTreeMap[mass][i]->setChannel(-1);
    signalNbr.push_back( signalTreeMap[mass][i]->Loop()/500./signalNorm[mass][i]);
//     cout << signalName[mass][i] <<" "<< signalNbr[i]<< " "<<fixed<<setprecision(3)
//        << signalNbr[i]/sqrt(totalBackground)<<endl;
    totalSignal+=signalNbr[i]*signalBR[mass][i];
  }



  dataTree[0]->setChannel(-1);
  float dataNbr = dataTree[0]->Loop();
  if (verbose) cout << "Data "<< dataNbr<<endl;

  cout << "Signal: "<< totalSignal <<"\t - Background: "
       << totalBackground
        <<"\t - Data: "       << dataNbr
       <<"\t - S/sqrt(B): "
<<fixed<<setprecision(4)
       << totalSignal/sqrt(totalBackground)<<endl;//<<"\t"

//   cout << totalSignal <<"\t\t"
//        << totalBackground<<"\t\t"
//        << dataNbr<<"\t\t"
//        << totalSignal/sqrt(totalBackground)<<endl;//<<"\t"
//        << totalSignal/sqrt(totalSignal+totalBackground)<<endl;

}


void getPromptSF(bool verbose = false) {
  for ( int j=0;j<3;++j) {
    backgroundsNbr.clear();

    float totalBackground = 0.;
    for (unsigned int i=0;i<backgroundsTree.size();++i) {
      backgroundsTree[i]->setChannel(j);
      backgroundsNbr.push_back( backgroundsTree[i]->Loop()*mcSF);
      if (verbose) cout << backgroundsName[i] <<" "<< backgroundsNbr[i]<<endl;
      totalBackground+=backgroundsNbr[i];
    }

    float totalFake = nonPB(j).first;
    if (verbose) cout << "Non-Prompt\t"<<fixed<<setprecision(2)
	 << totalFake << "\n";

    dataTree[0]->setChannel(j);
    float dataNbr = dataTree[0]->Loop();
    if (verbose) cout << "Data "<< dataNbr<<endl;

    cout << "Channel: " <<j<<endl;
    cout << " Data                : "<<dataNbr<<endl;
    cout << " Prompt Background   : "<<totalBackground<<endl;
    cout << " Non-Prompt estimate : "<<totalFake<<endl;
    cout << " Non-Prompt SF       : "<<(dataNbr-totalBackground)/totalFake<<endl;
    npSF[j]=(dataNbr-totalBackground)/totalFake;
  }
}
void getDYSF(bool verbose = false) {
  for ( int j=0;j<3;++j) {
    if (j==1) continue;
    backgroundsNbr.clear();

    float totalBackground = 0.;
    for (unsigned int i=0;i<backgroundsTree.size()-1;++i) {
      backgroundsTree[i]->setChannel(j);
      backgroundsNbr.push_back( backgroundsTree[i]->Loop()*mcSF);
      if (verbose) cout << backgroundsName[i] <<" "<< backgroundsNbr[i]<<endl;
      totalBackground+=backgroundsNbr[i];
    }

    backgroundsTree.back()->setChannel(j);
    float dy = backgroundsTree.back()->Loop()*mcSF;
    if (verbose) cout << "Drell Yan\t"<<fixed<<setprecision(2)
	 << dy << "\n";

    dataTree[0]->setChannel(j);
    float dataNbr = dataTree[0]->Loop();
    if (verbose) cout << "Data "<< dataNbr<<endl;

    cout << "Channel: " <<j<<endl;
    cout << " Data             : "<<dataNbr<<endl;
    cout << " Other Background : "<<totalBackground<<endl;
    cout << " DY estimate      : "<<dy<<endl;
    cout << " DY SF            : "<<(dataNbr-totalBackground)/dy<<endl;
    npSF[j]=(dataNbr-totalBackground)/dy;
  }
}

pair<float,float> getDYSFRoutinSlow(int ch)
{
  cout <<"\n\n DY SF estimate using the R_out/in method channel "<<ch<<"\n\n";
  
  // DY signal cuts:
  backgroundsTree.back()->setChannel(ch);

  setLeptonSumMass();
  double nevent_MC     = backgroundsTree.back()->Loop()*mcSF;
  setLeptonSumMass(76,106, false);
  double nevent_Out_MC = backgroundsTree.back()->Loop()*mcSF;
  double nevent_In_MC  = nevent_MC-nevent_Out_MC;

  double Routin = nevent_Out_MC/nevent_In_MC;
  cout <<"MC: \n";
  cout << "nevent total MC " << nevent_MC << endl;
  cout << "nevent in zpeak  " << nevent_In_MC << endl;
  cout << "nevent out zpeak " << nevent_Out_MC << endl;
  cout << "RAW Routin in MC " << Routin << endl;

  // Data signal cuts:

  dataTree[0]->setChannel(ch);
  setLeptonSumMass();
  double nevent_Data     = dataTree[0]->Loop();
  setLeptonSumMass(76,106, false);
  double nevent_Out_Data = dataTree[0]->Loop();
  double nevent_In_Data  = nevent_Data-nevent_Out_Data;

  cout <<"Data: \n";
  cout << "nevent total Data " << nevent_Data << endl;
  cout << "nevent in zpeak  " << nevent_In_Data << endl;
  cout << "nevent out zpeak " << nevent_Out_Data << endl;


  //Count N_in e-mu
  dataTree[0]->setChannel(1);
  setLeptonSumMassEM();
  double nevent_Data_EM     = dataTree[0]->Loop();
  setLeptonSumMassEM(76,106, false);
  double nevent_Out_Data_EM = dataTree[0]->Loop();
  double nevent_In_Data_EM  = nevent_Data_EM - nevent_Out_Data_EM;
  setLeptonSumMassEM();

  cout << "nevent total Data EM " << nevent_Data_EM << endl;
  cout << "nevent Out   Data EM    " << nevent_Out_Data_EM<< endl;
  cout << "nevent In    Data EM    " << nevent_In_Data_EM<< endl;

  //***********************
  //Calculation kee, kmumu
  //***********************
  
//   setLeptonSumMass();
//   dataTree[0]->setChannel(0);
//   double highMM = dataTree[0]->Loop();
//   dataTree[0]->setChannel(2);
//   double highEE = dataTree[0]->Loop();
// 
//   setLeptonSumMass(76,106, false);
//   dataTree[0]->setChannel(0);
//   double outhighMM = dataTree[0]->Loop();
//   dataTree[0]->setChannel(2);
//   double outhighEE = dataTree[0]->Loop();

  setLeptonSumMass();
  setMETcut(0., 99999., false);
  setHTcut(0, -1, false);
  setminMlbcut(-50, false);
  setSumCuts(0., 0.,0.,0.,0.,0., false);
  dataTree[0]->setChannel(0);
  double allMM = dataTree[0]->Loop();
  dataTree[0]->setChannel(2);
  double allEE = dataTree[0]->Loop();

  setLeptonSumMass(76,106, false);
  dataTree[0]->setChannel(0);
  double outMM = dataTree[0]->Loop();
  dataTree[0]->setChannel(2);
  double outEE = dataTree[0]->Loop();


  double lowInMM = allMM - outMM;// - highMM - outMM + outhighMM;
  double lowInEE = allEE - outEE;// - highEE - outEE + outhighEE;

  cout << "Calculation k MM " << lowInMM<<" "<<allMM <<" "<< outMM <<"\n ";//<<highMM <<" "<<  outhighMM<<endl;
  cout << "Calculation k EE " << lowInEE<<" "<<allEE <<" "<< outEE <<"\n ";//<<highEE <<" "<<  outhighEE<<endl;
  double kFact = 0;
  if (ch==0) kFact = sqrt(lowInMM/lowInEE);
  if (ch==2) kFact = sqrt(lowInEE/lowInMM);

  cout << "kFact " << kFact << endl;

  // Rin/out from MC:
  setBJets(0,0);
  setminMlbcut(-50, false);
//    setSumCuts(0., 0.,0.,0.,300.,0., false);
  dataTree[0]->setChannel(ch);
  backgroundsTree.back()->setChannel(ch);
  backgroundsTree.front()->setChannel(ch);
  setLeptonSumMass();
  setMETcut(0., 99999., false);
  double All_Data = dataTree[0]->Loop();
  double All_MC = backgroundsTree.back()->Loop();
//   double All_tt = backgroundsTree.front()->Loop();
  setMETcut(10, 99999., false);
  double HighData = dataTree[0]->Loop();
  double HighMC = backgroundsTree.back()->Loop();
//   double Hightt= backgroundsTree.front()->Loop();
  setLeptonSumMass(76,106, false);
  double OutHi_Data = dataTree[0]->Loop();
  double OutHi_MC = backgroundsTree.back()->Loop();
//   double OutHi_tt = backgroundsTree.front()->Loop();
  setMETcut(0., 99999., false);
  double Out_Data = dataTree[0]->Loop();
  double Out_MC = backgroundsTree.back()->Loop();
//   double Out_tt = backgroundsTree.front()->Loop();

  double OutLow_Data = Out_Data - OutHi_Data;
  double InLow_Data  = All_Data - HighData - Out_Data + OutHi_Data;

  double OutLow_MC = Out_MC - OutHi_MC;
  double InLow_MC  = All_MC - HighMC - Out_MC + OutHi_MC;
//   double OutLow_tt = Out_tt - OutHi_tt;
//   double InLow_tt  = All_tt - Hightt - Out_tt + OutHi_tt;

  double CFR = (OutLow_Data/InLow_Data)/(OutLow_MC/InLow_MC);
  cout << "Calculation R Data " << All_Data <<" "<< HighData <<" "<< Out_Data <<" "<< OutHi_Data<<endl;
  cout << "Calculation R MC   " << All_MC <<" "<< HighMC <<" "<< Out_MC <<" "<< OutHi_MC<<endl;
//   cout << "Calculation CFR tt " << All_tt <<" "<< Hightt <<" "<< Out_tt <<" "<< OutHi_tt<<endl;
//   cout << "Calculation CFR  " << OutLow_Data <<" "<< InLow_Data <<" "<< OutLow_tt <<" "<< InLow_tt<<endl;
  cout << "R out/in Data control region " << (OutLow_Data/InLow_Data) <<" +/-"<< sqrt(OutLow_Data+OutLow_Data*OutLow_Data/InLow_Data)/InLow_Data<<endl;
  cout << "R out/in MC   control region " << (OutLow_MC/InLow_MC) <<" +/-"<< sqrt(OutLow_MC+OutLow_MC*OutLow_MC/InLow_MC)/InLow_MC<<endl;

  cout << "CFR " << CFR << endl;
  cout << "Corected R out/in MC" << Routin*CFR << endl;


  double theEstimate = Routin*CFR*(nevent_In_Data - 0.5*nevent_In_Data_EM*kFact);
  double theEstimateData = (OutLow_Data/InLow_Data)*(nevent_In_Data - 0.5*nevent_In_Data_EM*kFact);

  cout << "estimated (Data) Z contamination using R out/in MC  , Channel " << ch << " is " << theEstimate << endl;//" +/- \n";// << stat_error << " (stat) " << endl;
  cout << "estimated (Data) Z contamination using R out/in Data, Channel " << ch << " is " << theEstimateData << endl;//<< " +/- \n";// << stat_error << " (stat) " << endl;
 // ofile << "ZSF is  " << channel << " is " <<  theEstimate/ << " +/- " << stat_error << " (stat) " << endl;
  cout << "expected (MC) Z contamination for channel    " << ch << " is " << nevent_Out_MC << endl;//<<  " (stat) \n";//  << nevent_OutErr_MC << endl;
  cout << "The SF: "<<theEstimateData/nevent_Out_MC<<endl;

  return FFPair(theEstimateData, theEstimateData/nevent_Out_MC);

}

pair<float,float> getDYSFRoutin(int ch)
{
  cout <<"\n\n DY SF estimate using the R_out/in method channel "<<ch<<"\n\n";
  
  // DY signal cuts:
  backgroundsTree.back()->setChannel(ch);

  setLeptonSumMass();
  double nevent_MC     = backgroundsTree.back()->Loop()*mcSF;
  setLeptonSumMass(76,106, false);
  double nevent_Out_MC = backgroundsTree.back()->Loop()*mcSF;
  double nevent_In_MC  = nevent_MC-nevent_Out_MC;

  double Routin = nevent_Out_MC/nevent_In_MC;
  cout <<"MC: \n";
  cout << "nevent total MC " << nevent_MC << endl;
  cout << "nevent in zpeak  " << nevent_In_MC << endl;
  cout << "nevent out zpeak " << nevent_Out_MC << endl;
  cout << "RAW Routin in MC " << Routin << endl;


  // Data signal cuts:

  dataTree[0]->setChannel(ch);
  setLeptonSumMass();
  double nevent_Data     = dataTree[0]->Loop();
  setLeptonSumMass(76,106, false);
  double nevent_Out_Data = dataTree[0]->Loop();
  double nevent_In_Data  = nevent_Data-nevent_Out_Data;

  cout <<"Data: \n";
  cout << "nevent total Data " << nevent_Data << endl;
  cout << "nevent in zpeak  " << nevent_In_Data << endl;
  cout << "nevent out zpeak " << nevent_Out_Data << endl;


  //Count N_in e-mu
  dataTree[0]->setChannel(1);
  setLeptonSumMassEM();
  double nevent_Data_EM     = dataTree[0]->Loop();
  setLeptonSumMassEM(76,106, false);
  double nevent_Out_Data_EM = dataTree[0]->Loop();
  double nevent_In_Data_EM  = nevent_Data_EM - nevent_Out_Data_EM;
  setLeptonSumMassEM();

  cout << "nevent total Data EM " << nevent_Data_EM << endl;
  cout << "nevent In Data EM    " << nevent_In_Data_EM<< endl;

  //***********************
  //Calculation kee, kmumu
  //***********************
  
//   setLeptonSumMass();
//   dataTree[0]->setChannel(0);
//   double highMM = dataTree[0]->Loop();
//   dataTree[0]->setChannel(2);
//   double highEE = dataTree[0]->Loop();
// 
//   setLeptonSumMass(76,106, false);
//   dataTree[0]->setChannel(0);
//   double outhighMM = dataTree[0]->Loop();
//   dataTree[0]->setChannel(2);
//   double outhighEE = dataTree[0]->Loop();

  setLeptonSumMass();
  setMETcut(0., 99999., false);
   setHTcut(0, -1, false);
  setminMlbcut(-50, false);
  setSumCuts(0., 0.,0.,0.,0.,0., false);
  dataTree[0]->setChannel(0);
  double allMM = dataTree[0]->Loop();
  dataTree[0]->setChannel(2);
  double allEE = dataTree[0]->Loop();

  setLeptonSumMass(76,106, false);
  dataTree[0]->setChannel(0);
  double outMM = dataTree[0]->Loop();
  dataTree[0]->setChannel(2);
  double outEE = dataTree[0]->Loop();


  double lowInMM = allMM - outMM;// - highMM - outMM + outhighMM;
  double lowInEE = allEE - outEE;// - highEE - outEE + outhighEE;

  cout << "Calculation k MM " << lowInMM<<" "<<allMM <<" "<< outMM <<"\n ";//<<highMM <<" "<<  outhighMM<<endl;
  cout << "Calculation k EE " << lowInEE<<" "<<allEE <<" "<< outEE <<"\n ";//<<highEE <<" "<<  outhighEE<<endl;
  double kFact = 0;
  if (ch==0) kFact = sqrt(lowInMM/lowInEE);
  if (ch==2) kFact = sqrt(lowInEE/lowInMM);

  cout << "kFact " << kFact << endl;

  // Rin/out from data in control regio, 0 btags, 0<MET<10 GeV:
  setBJets(0,0);
  setminMlbcut(-50, false);
  dataTree[0]->setChannel(ch);
  setLeptonSumMass();
  setMETcut(0., 10., false);
  double AllLowData = dataTree[0]->Loop();
  setLeptonSumMass(76,106, false);
  
  setMETcut(0., 10., false);
  double OutLow_Data = dataTree[0]->Loop();
  double InLow_Data  = AllLowData - OutLow_Data;//All_Data - HighData -Out_Data + OutHi_Data;

  cout << "Calculation R Data " << AllLowData <<" "<< OutLow_Data<<" "<< OutLow_Data <<" "<< InLow_Data<<endl;
  cout << "R out/in Data control region " << (OutLow_Data/InLow_Data) <<" +/-"<< sqrt(OutLow_Data+OutLow_Data*OutLow_Data/InLow_Data)/InLow_Data<<endl;

  double theEstimateData = (OutLow_Data/InLow_Data)*(nevent_In_Data - 0.5*nevent_In_Data_EM*kFact);
  double stat_error = sqrt(1/OutLow_Data+1/InLow_Data +1/nevent_In_Data) * (OutLow_Data/InLow_Data)*(nevent_In_Data);

  cout << "estimated (Data) Z contamination using R out/in Data " << ch << " is " << theEstimateData << " +/- " << stat_error << " (stat) " << endl;
  cout << "The SF: "<<theEstimateData/nevent_Out_MC<<endl;
  return FFPair(theEstimateData, theEstimateData/nevent_Out_MC);
}

pair<float,float> getDYestimate(int ch, double R)
{
  cout <<"\n\n DY SF estimate using R_out/in of "<<R<<" in channel "<<ch<<"\n\n";
  
  // DY signal cuts:
  backgroundsTree.back()->setChannel(ch);

  setLeptonSumMass();
  double nevent_MC     = backgroundsTree.back()->Loop()*mcSF;
  setLeptonSumMass(76,106, false);
  double nevent_Out_MC = backgroundsTree.back()->Loop()*mcSF;
  double nevent_In_MC  = nevent_MC-nevent_Out_MC;

  cout <<"MC: \n";
  cout << "nevent total MC " << nevent_MC << endl;
  cout << "nevent in zpeak  " << nevent_In_MC << endl;
  cout << "nevent out zpeak " << nevent_Out_MC << endl;


  // Data signal cuts:

  dataTree[0]->setChannel(ch);
  setLeptonSumMass();
  double nevent_Data     = dataTree[0]->Loop();
  setLeptonSumMass(76,106, false);
  double nevent_Out_Data = dataTree[0]->Loop();
  double nevent_In_Data  = nevent_Data-nevent_Out_Data;

  cout <<"Data: \n";
  cout << "nevent total Data " << nevent_Data << endl;
  cout << "nevent in zpeak  " << nevent_In_Data << endl;
  cout << "nevent out zpeak " << nevent_Out_Data << endl;


  //Count N_in e-mu
  dataTree[0]->setChannel(1);
  setLeptonSumMassEM();
  double nevent_Data_EM     = dataTree[0]->Loop();
  setLeptonSumMassEM(76,106, false);
  double nevent_Out_Data_EM = dataTree[0]->Loop();
  double nevent_In_Data_EM  = nevent_Data_EM - nevent_Out_Data_EM;
  setLeptonSumMassEM();

  cout << "nevent total Data EM " << nevent_Data_EM << endl;
  cout << "nevent In Data EM    " << nevent_In_Data_EM<< endl;

  //***********************
  //Calculation kee, kmumu
  //***********************
  
//   setLeptonSumMass();
//   dataTree[0]->setChannel(0);
//   double highMM = dataTree[0]->Loop();
//   dataTree[0]->setChannel(2);
//   double highEE = dataTree[0]->Loop();
// 
//   setLeptonSumMass(76,106, false);
//   dataTree[0]->setChannel(0);
//   double outhighMM = dataTree[0]->Loop();
//   dataTree[0]->setChannel(2);
//   double outhighEE = dataTree[0]->Loop();

  setLeptonSumMass();
  setMETcut(0., 99999., false);
   setHTcut(0, -1, false);
  setminMlbcut(-50, false);
  setSumCuts(0., 0.,0.,0.,0.,0., false);
  dataTree[0]->setChannel(0);
  double allMM = dataTree[0]->Loop();
  dataTree[0]->setChannel(2);
  double allEE = dataTree[0]->Loop();

  setLeptonSumMass(76,106, false);
  dataTree[0]->setChannel(0);
  double outMM = dataTree[0]->Loop();
  dataTree[0]->setChannel(2);
  double outEE = dataTree[0]->Loop();


  double lowInMM = allMM - outMM;// - highMM - outMM + outhighMM;
  double lowInEE = allEE - outEE;// - highEE - outEE + outhighEE;

  cout << "Calculation k MM " << lowInMM<<" "<<allMM <<" "<< outMM <<"\n ";//<<highMM <<" "<<  outhighMM<<endl;
  cout << "Calculation k EE " << lowInEE<<" "<<allEE <<" "<< outEE <<"\n ";//<<highEE <<" "<<  outhighEE<<endl;
  double kFact = 0;
  if (ch==0) kFact = sqrt(lowInMM/lowInEE);
  if (ch==2) kFact = sqrt(lowInEE/lowInMM);


  double theEstimateData = R*(nevent_In_Data - 0.5*nevent_In_Data_EM*kFact);

  cout << "estimated (Data) Z contamination using R out/in Data " << ch << " is " << theEstimateData << endl;
  cout << "The SF: "<<theEstimateData/nevent_Out_MC<<endl;
  

  return FFPair(theEstimateData, theEstimateData/nevent_Out_MC);


  
}

#define prntVal(val,err)  cout << " & "<< right<<setw(9) <<fixed<<setprecision(2)<< val; if (tex && doUncert) cout << " $\\pm$ "<< left<<setw(8) <<fixed<<setprecision(2)<< err;
#define prntValOnly(val)  cout << " & "<< right<<setw(9) <<fixed<<setprecision(2)<< val; if (tex && doUncert) cout << setw(13) <<" ";

void yieldTable(bool doFake = false, bool doUncert = false) {
  if (tex) {
    cout <<"\\begin{table}[htp]\n";
    cout <<"\\begin{center}\n";
    cout <<"\\caption{Number of events after all cuts}\n";
    cout <<"\\begin{tabular}{|l|c|c|c|c|}\\hline\\hline\n";
    cout << "Sample \t\t& $\\mu\\mu$\t& e$\\mu$\t& ee \t& Sum \\\\\n\\hline\n";
  } else{
    cout <<"Sample\t\t\t MuMu\t ElMu\t ElEl\t All\n";
  }

    //   cout << "123456789x123456789x123456789x123456789x123456789x123456789x123456789x123456789x\n";

//   float nominal=0;
//   for (unsigned int i=0;i<signalTreeMap[mass].size();++i) {
//     float total = 0.;
//     float total_unc = 0.;
//     cout.width(16);
//     cout << left<< signalName[mass][i];
//     for ( int j=0;j<3;++j) {
//       signalTreeMap[mass][i]->setChannel(j);
//       pair<int,float> res = signalTreeMap[mass][i]->LoopDouble( );
//       float n = res.second/500./signalNorm[mass][i];
//       //float n = signalTreeMap[mass][i]->LoopDouble().first;
//       float err2 = 0.;
//       if (res.first>0.) err2= n*sqrt( 0.0036 + 1./res.first);
//       nominal+=n*signalBR[mass][i];
//       total+=n;
//       if (res.first>0.) total_unc+= n*n/res.first;
//       prntVal(n,err2);
//       //cout << " , "<< right<<setw(6) <<fixed<<setprecision(2)<< n;
// //       if (tex && doUncert) cout << " \\pm "<< left<<setw(6) <<fixed<<setprecision(2)<< n*0.06;
//     }
//     prntVal(total, sqrt(total_unc+(total*total*0.0036)));
//     if (tex) cout    <<" \\\\";
//     cout <<endl;
//   }
//   return;
//    cout << "Nominal "<< mass <<" : "<<fixed<<setprecision(4)<<  nominal<<endl;
  float totalBackground[4] = {0.,0.,0.,0.};
  float totalBackground_Error[4] = {0.,0.,0.,0.};
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    //if (i<7 && i!=3) continue;
    float total = 0.;
    float total_unc = 0.;
    float nbr, nbr_error;
    cout << setw(16)<<left<<backgroundsName[i];
    for ( int j=0;j<3;++j) {
      backgroundsTree[i]->setChannel(j);
      pair<int,float> res;
      if ((sign==-1)&&(i==(backgroundsTree.size()-1))&&(dyEstimate[0] != dyEstimate[2])) {
	if (j==0) res = pair<int,float> (10000,dyEstimate[0]);
	else if (j==2) res = pair<int,float> (10000,dyEstimate[2]);
	else res = backgroundsTree[i]->LoopDouble( );
      } else {
	res = backgroundsTree[i]->LoopDouble( );
      }
      nbr=res.second*mcSF;//*(sign==-1&&i==(backgroundsTree.size()-1)&&(j!=-1)?npSF[j]:1.0);
      nbr_error = 0;
      if (res.first>0.) {
        nbr_error = nbr*sqrt( backgroundsUnc[i]*backgroundsUnc[i] + 1./res.first);// nbr*backgroundsUnc[i];
        total_unc += (nbr*nbr/res.first);
      }
//       cout << res.first<<" "<<backgroundsUnc[i]<<" "<<sqrt( backgroundsUnc[i]*backgroundsUnc[i] + 1./res.first)<<" "
//       <<backgroundsUnc[i]*backgroundsUnc[i] <<" "
//       <<1./res.first<<" ";
      total+=nbr;
      totalBackground[j]+=nbr;
      totalBackground_Error[j]+=nbr_error*nbr_error;
      prntVal(nbr, nbr_error);
    }
    totalBackground[3]+=total;
    total_unc += (total*total*backgroundsUnc[i]*backgroundsUnc[i]);
    totalBackground_Error[3]+=total_unc;
    prntVal(total, sqrt(total_unc));
    if (tex) cout    <<" \\\\";
    cout <<endl;
  }

  if (tex) cout <<"\\hline\n";
  if (nonPrompt!=0) {
    cout << setw(16)<<"Total Prompt MC";
    for ( int j=0;j<4;++j) {
      prntVal(totalBackground[j], sqrt(totalBackground_Error[j]));
    }
    if (tex) cout    <<" \\\\";
    cout <<endl;

    float totalFake=0.0, tfe = 0.0, tf=0.0;
    cout << "Non-Prompt\t";
    for ( int j=0;j<3;++j) {
      FFPair fake = nonPB(j);
      if (fake.first>0.0) {
	float errsq = fake.second*fake.second + fake.first*fake.first*0.25;
	prntVal(fake.first, sqrt(errsq));
	//cout << " " << fake.second<<" "<<sqrt(errsq)<<" "<<sqrt(fake.first+fake.first*fake.first*0.25)<<" ";
	totalBackground[j]+=fake.first;
	totalBackground[3]+=fake.first;
	tf+=fake.first;

	totalBackground_Error[j]+=errsq;
	totalBackground_Error[3]+=fake.second*fake.second;
	tfe+=fake.second*fake.second;
      } else { prntVal(0., 0.);}
    }
    totalBackground_Error[3]+=tf*tf*0.25;
    tfe+=tf*tf*0.25;
    prntVal(tf, sqrt(tfe));
    if (tex)cout <<" \\\\\n\\hline";
    
  }

  cout << "\nTotal MC\t";
  for ( int j=0;j<4;++j) {
    prntVal(totalBackground[j], sqrt(totalBackground_Error[j]));
  }
  if (tex) cout    <<" \\\\";
  cout <<endl;

  if (tex)cout <<"\\hline\n";

  float data[4];
  for (unsigned int i=0;i<dataTree.size();++i) {
    data[3] = 0.;
    cout <<setw(16)<<"Data";
    for ( int j=0;j<3;++j) {
      dataTree[i]->setChannel(j);
      float nbr=dataTree[i]->Loop();
      data[j]=nbr;
      data[3]+=nbr;
      prntValOnly(nbr);
    }
    prntValOnly(data[3]);
    if (tex) cout    <<" \\\\";
    cout <<endl;
  }
  if (tex)cout <<"\\hline\n";
  cout << setw(16)<<"Data/MC";
  for ( int j=0;j<4;++j) {
    prntValOnly(data[j]/totalBackground[j]);
  }
  if (tex) cout    <<" \\\\";
  cout <<endl;
  if (tex){
    cout <<"\\end{tabular}\n";
    cout <<"\\end{center}\n";
    cout <<"\\end{table}\n";
  }

}

#define prntValBinomial(val,xS,init)  cout << " & "<< fixed<<setprecision(3)<< val.second/(500*xS*19.6*10) << " $\\pm$ "<< left<<setw(4) <<fixed<<setprecision(3)<< sqrt( val.first*(init-val.first)/(init*init*init) )*100;

void efficiencyTable() {

   for (int mass=500;mass<=1500;mass+=100) {
      float nominal=0.;


  cout << left<<setw(4)<< mass<<" \\GeVcc";
  for (unsigned int i=0;i<signalTreeMap[mass].size();++i) {
      cout.width(5);
      signalTreeMap[mass][i]->setChannel(3);
      pair<int,float> res = signalTreeMap[mass][i]->LoopDouble( );
      //float n = res.second/500./signalNorm[mass][i];
      float n = signalTreeMap[mass][i]->LoopDouble().second;
//       float err2 = n*sqrt( 0.0036 + 1./res.first);
      nominal+=n*signalBR[mass][i];
//       total_unc+=err2*err2;
//  cout << n << " " << initialNumber[string(signalName[mass][i])]
//  <<" " << n/initialNumber[string(signalName[mass][i])]
//  <<" " << res.second/500 
//  <<" " << (xSection[mass]*19.6*1000)
//  <<" " << res.second /500/(xSection[mass]*19.6*10)
//   <<endl;
      prntValBinomial(res,xSection[mass],initialNumber[string(signalName[mass][i])]);
      //cout << " , "<< right<<setw(6) <<fixed<<setprecision(2)<< n;
//       if (tex && doUncert) cout << " \\pm "<< left<<setw(6) <<fixed<<setprecision(2)<< n*0.06;
    }
//     prntVal(total, sqrt(total_unc));
    cout  << "  "<<nominal/(500*xSection[mass]*19.6*10) << "  "<<nominal/500*(newXSection[mass]/xSection[mass])*mcSF <<" \\\\";
    cout <<endl;
  }
 }

int nbrCutFlow = 5;

void setCutFlowCut(int j)
{

  switch( j )
  {
    case -1: setJets(0, false);	setBJets(0, 100, false); setMETcut(0., 99999., false); setHTcut(0, -1, false); break;
    case 0: setJets(2, false);	break;
    case 1: setJets(3, false);	break;
    case 2: setJets(3, false);	setBJets(1, 100, false); break;
    case 3: setJets(3, false);	setBJets(1, 100, false); setMETcut(30, 99999., false); break;
    case 4: setJets(3, false);	setBJets(1, 100, false); setMETcut(30, 99999., false); setHTcut(300, -1, false); break;

//     case -1: setJets(0, false); setBJets(0, 100, false); setMETcut(0., 99999., false); setHTcut(0, -1, false); break;
   //case 9: setLeptonSumMass(76.0, 106.0);break;
//    case 0: setJets(2, false);  break;
//    case 1: setJets(3, false);  break;
//    case 2: setJets(3, false);  setBJets(1, 100, false); break;
//    case 3: setJets(3, false);  setBJets(1, 100, false); setMETcut(30, false); break;
//    case 4: setJets(3, false);  setBJets(1, 100, false); setMETcut(30, false); setHTcut(300, false); break;
//    case 5: setJets(3, false);  setBJets(1, 100, false); setMETcut(30, false); setHTcut(300, false); setminMlbcut(170, false); break;
//    case 6: setJets(3, false);  setBJets(2, 100, false); setMETcut(30, false); setHTcut(300, false); setminMlbcut(170, false); break;
   //case 7: setJets(3, false);  setBJets(2, 100, false); setMETcut(30, false); setHTcut(0.0, false); break;

  }




}

void cutFlowTable(int channel, bool efficiency= false, bool doFake = false, bool verbose = false) {

  int precision = 1;
  if (efficiency) precision=3;
  if (tex) {
    cout <<"\\begin{table}[htp]\n";
    cout <<"\\begin{center}\n";
    cout <<"\\caption{Number of events}\n";
    cout <<"\\begin{tabular}{|l|";
    for ( int j=0;j<nbrCutFlow;++j) cout <<"c|";
    cout <<"}\\hline\\hline\n";
    cout << "Sample \t";
    for ( int j=0;j<nbrCutFlow;++j) cout <<"& cut "<< j;
    cout <<" \\\\\n\\hline\n";
  } else{
    cout <<"Sample\t\t";
    for ( int j=0;j<nbrCutFlow;++j) cout <<"\t cut "<< j;
    cout <<endl;
  }

  float initial= 1.;
  for (int mass=500;mass<=1500;mass+=100) {
    for (unsigned int i=0;i<signalTreeMap[mass].size();++i) {
      cout << signalName[mass][i]<< "\t";
      signalTreeMap[mass][i]->setChannel(channel);
      setCutFlowCut(-1);

      for ( int j=0;j<nbrCutFlow;++j) {
	setCutFlowCut(j);
	float n = signalTreeMap[mass][i]->Loop()/500./signalNorm[mass][i];
	if (efficiency && (j==0)) initial=n;
	//total+=n;
	cout <<fixed<<setprecision(precision)<< "\t& "<< n/initial;
      }
      //cout <<fixed<<setprecision(precision)<< "\t& "<< total;
      if (tex) cout    <<" \\\\";
      cout <<endl;
    }
  }

  float totalBackground[20] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->setChannel(channel);
    setCutFlowCut(-1);
    //float total = 0.;
    cout << backgroundsName[i]<< "\t";
    for ( int j=0;j<nbrCutFlow;++j) {
      setCutFlowCut(j);
      float nbr=backgroundsTree[i]->Loop()*mcSF*(sign==-1&&i==(backgroundsTree.size()-1)&&(channel!=-1)?npSF[channel]:1.0);
      if (efficiency && (j==0)) initial=nbr;
//       total+=nbr;
      totalBackground[j]+=nbr;
      cout <<fixed<<setprecision(precision) << "\t& "<< nbr/initial;
    }
//     totalBackground[3]+=total;
//     cout <<fixed<<setprecision(precision)<< "\t&  "<< total;
    if (tex) cout    <<" \\\\";
    cout <<endl;
  }

  if (nonPrompt!=0) {
    if (tex) cout <<"\\hline\n";
    cout << "Total Prompt MC\t";
    for ( int j=0;j<nbrCutFlow;++j) {
      if (efficiency) cout <<fixed<<setprecision(precision) << "\t& "<< totalBackground[j]/totalBackground[0];
      else cout <<fixed<<setprecision(precision) << "\t& "<< totalBackground[j];
    }
    if (tex) cout    <<" \\\\";
    cout <<endl;

    setCutFlowCut(-1);
    float totalFake[3];
    cout << "Non-Prompt\t";
    // totalBackground[3]+=totalFake[0]+totalFake[1]+totalFake[2];
    for ( int j=0;j<nbrCutFlow;++j) {
      setCutFlowCut(j);
      nonPB(totalFake);
      float fake=0;
      if (channel) fake=totalFake[0]+totalFake[1]+totalFake[2];
      else fake = totalFake[channel];
      if (efficiency && (j==0)) initial=fake;
      cout <<fixed<<setprecision(precision) << "\t& "<< fake/initial;
      totalBackground[j]+=fake;
    }
    cout <<endl;
//     cout <<"\t& "<<fixed<<setprecision(precision)
// 	 << totalFake[0]+totalFake[1]+totalFake[2]<< " \\\\\n";
  }
  if (tex)cout <<"\\hline\n";

  cout << "Total Bkg\t";
  for ( int j=0;j<nbrCutFlow;++j) {
    if (efficiency) cout <<fixed<<setprecision(precision) << "\t& "<< totalBackground[j]/totalBackground[0];
    else cout <<fixed<<setprecision(precision) << "\t& "<< totalBackground[j];
  }
  if (tex) cout    <<" \\\\";
  cout <<endl;

  if (tex)cout <<"\\hline\n";

  setCutFlowCut(-1);
  cout <<"Data\t\t ";
  dataTree[0]->setChannel(channel);
  for ( int j=0;j<nbrCutFlow;++j) {
    setCutFlowCut(j);
    float nbr=dataTree[0]->Loop();
    if (efficiency && (j==0)) initial=nbr;
    cout <<fixed<<setprecision(precision) << "\t& "<< nbr/initial;
  }
  if (tex) cout    <<" \\\\";
  cout <<endl;
  if (tex)cout <<"\\hline\n";
//   cout << "Data/MC : \t ";
//   for ( int j=0;j<4;++j) {
//     cout <<fixed<<setprecision(precision) << "\t& "<< data[j]/totalBackground[j];
//   }
//   if (tex) cout    <<" \\\\";
//   cout <<endl;
  if (tex){
    cout <<"\\end{tabular}\n";
    cout <<"\\end{center}\n";
    cout <<"\\end{table}\n";
  }

}



TH1F* nonPromptPlotSS(double f1, double f2, double p1, double p2){


  nonPrompt->setEventType(0);
  nonPrompt->Loop();
  TH1F* Nt00Histo = (TH1F*)nonPrompt->newPlot->Clone("nt00");
  string histoName = nonPrompt->newPlot->GetTitle();
  nonPrompt->setEventType(1);
  nonPrompt->Loop();
  TH1F* Nt01Histo = (TH1F*)nonPrompt->newPlot->Clone("nt01");
  nonPrompt->setEventType(2);
  nonPrompt->Loop();
  TH1F* Nt10Histo = (TH1F*)nonPrompt->newPlot->Clone("nt10");
  nonPrompt->setEventType(3);
  nonPrompt->Loop();
  TH1F* Nt11Histo = (TH1F*)nonPrompt->newPlot->Clone("nt11");
  delete nonPrompt->newPlot;
  nonPrompt->newPlot = 0;


    double det = 1.0 / ((p1 - f1)*(p2 - f2));


    TH1F* hNff = (TH1F*)Nt00Histo->Clone("hNff");
    hNff->Scale(p1*p2);
    hNff->Add(Nt10Histo, -1*(1-p1)*p2);
    hNff->Add(Nt01Histo, -1*p1*(1-p2));
    hNff->Add(Nt11Histo, (1-p1)*(1-p2));
    hNff->Scale(det*f1*f2);

    TH1F* hNpf = (TH1F*)Nt00Histo->Clone("hNpf");
    hNpf->Scale(-1*f1*p2);
    hNpf->Add(Nt10Histo, (1-f1)*p2);
    hNpf->Add(Nt01Histo, f1*(1-p2));
    hNpf->Add(Nt11Histo, -1*(1-f1)*(1-p2));
    hNpf->Scale(det*p1*f2);

    TH1F* hNfp = (TH1F*)Nt00Histo->Clone("hNfp");
    hNfp->Scale(-1*p1*f2);
    hNfp->Add(Nt10Histo, (1-p1)*f2);
    hNfp->Add(Nt01Histo, p1*(1-f2));
    hNfp->Add(Nt11Histo, -1*(1-p1)*(1-f2));
    hNfp->Scale(det*f1*p2);

    TH1F* goodHisto = (TH1F*)hNff->Clone(histoName.c_str());
    goodHisto->SetNameTitle(histoName.c_str(), histoName.c_str());
    goodHisto->Add(hNpf);
    goodHisto->Add(hNfp);

    delete hNff;
    delete hNpf;
    delete hNfp;
    delete Nt00Histo;
    delete Nt10Histo;
    delete Nt01Histo;
    delete Nt11Histo;
    return goodHisto;
}

TH1F* nonPromptPlotTRI(double f, double p){


  char name[300];
  nonPrompt->setEventType(0);
  nonPrompt->Loop();
  string histoName = nonPrompt->newPlot->GetName();
  sprintf(name,"nt0_%s",histoName.c_str());
  TH1F* Nt0Histo = (TH1F*)nonPrompt->newPlot->Clone(name);
  nonPrompt->setEventType(1);
  nonPrompt->Loop();
  sprintf(name,"nt1_%s",histoName.c_str());
  TH1F* Nt1Histo = (TH1F*)nonPrompt->newPlot->Clone(name);
  nonPrompt->setEventType(2);
  nonPrompt->Loop();
  sprintf(name,"nt2_%s",histoName.c_str());
  TH1F* Nt2Histo = (TH1F*)nonPrompt->newPlot->Clone(name);
  nonPrompt->setEventType(3);
  nonPrompt->Loop();
  sprintf(name,"nt3_%s",histoName.c_str());
  TH1F* Nt3Histo = (TH1F*)nonPrompt->newPlot->Clone(name);
  cout <<"a\n";
  delete nonPrompt->newPlot;
  cout << Nt0Histo->Integral()<< " ";
  cout << Nt1Histo->Integral()<< " ";
  cout << Nt2Histo->Integral()<< " ";
  cout << Nt3Histo->Integral()<< " \n";


    double det = 1.0/((p-f)*(p-f)*(p-f));

    TH1F* hNfff = (TH1F*)Nt0Histo->Clone("hNfff");
    hNfff->Scale(p*p*p);
    hNfff->Add(Nt1Histo, -p*p*(1-p));
    hNfff->Add(Nt2Histo, p*(1-p)*(1-p));
    hNfff->Add(Nt3Histo, -(1-p)*(1-p)*(1-p));
    hNfff->Scale(det);


    TH1F* hNffp = (TH1F*)Nt0Histo->Clone("hNffp");
    hNffp->Scale(-3*p*p*f);
    hNffp->Add(Nt1Histo, 2*p*f*(1-p)+p*p*(1-f));
    hNffp->Add(Nt2Histo, -(f*(1-p)*(1-p) + 2*p*(1-p)*(1-f)));
    hNffp->Add(Nt3Histo, 3*(1-p)*(1-p)*(1-f));
    hNffp->Scale(det);

    TH1F* hNfpp = (TH1F*)Nt0Histo->Clone("hNfpp");
    hNfpp->Scale(3*p*f*f);
    hNfpp->Add(Nt1Histo, -(f*f*(1-p) + 2*p*f*(1-f)));
    hNfpp->Add(Nt2Histo, (2*f*(1-p)*(1-f) + p*(1-f)*(1-f)));
    hNfpp->Add(Nt3Histo, -3*(1-p)*(1-f)*(1-f));
    hNfpp->Scale(det);

    TH1F* goodHisto = (TH1F*)hNfff->Clone(histoName.c_str());
    goodHisto->SetNameTitle(histoName.c_str(), histoName.c_str());
    goodHisto->Add(hNffp);
    goodHisto->Add(hNfpp);

    delete hNfff;
    delete hNffp;
    delete hNfpp;
    Nt0Histo->Write();
    Nt1Histo->Write();
    Nt2Histo->Write();
    Nt3Histo->Write();
    return goodHisto;
}

void fillBin(TH1F *plot, int bin, float val, int N)
{
  plot->SetBinContent(bin,val);
  if (N!=0) plot->SetBinError(bin,val/sqrt(N));  
}

void countingPlot(string sampleName, string rootFileName, int offset=0, bool singleHisto= true)
{
  int signBin = (1-sign);
  int bin;

  TFile* outFile = new TFile(rootFileName.c_str(), "RECREATE");
  TH1F *sigPlot;

  char name[300];
  for (int mass=500;mass<=1500;mass+=100) {
    for (unsigned int i=0;i<signalTreeMap[mass].size();++i) {
      if (!singleHisto) {
	sprintf(name,"yield_test_%s_All",signalName[mass][i].c_str());
	sigPlot = new TH1F(name,name, 12, 0., 12.);
      }
//       sigPlot->SetDefaultSumw2();
//       sigPlot->Sumw2();
      for ( int j=0;j<3;++j) {
	signalTreeMap[mass][i]->setChannel(j);
	pair<int,float> res = signalTreeMap[mass][i]->LoopDouble( );
	bin = signBin*3+j+offset+1;
	if (singleHisto) {
	  sprintf(name,"yield_test_%s_All_%ix",signalName[mass][i].c_str(),bin-1);
	  sigPlot = new TH1F(name,name, 1, 0., 1.);
	  bin=1;
	}
	fillBin(sigPlot, bin, res.second/500./signalNorm[mass][i], res.first);
// 	sigPlot->SetBinContent,res.second/500./signalNorm[mass][i]);
// 	if (res.first!=0) sigPlot->SetBinError(signBin*3+j+offset+1,res.second/500./signalNorm[mass][i]/sqrt(res.first));
// 	cout << name<<" "<<res.first<<" "<<res.second/500./signalNorm[mass][i]<<" "<<res.second/500./signalNorm[mass][i]/sqrt(res.first)<<endl;
	//sigPlot->Fill( signBin*3+j+offset, signalTreeMap[mass][i]->Loop()/500./signalNorm[mass][i]);
      }
    }
  }

  TH1F *topPlot, *bckgPlot,*dataPlot, *npPlot;
  if (!singleHisto) {
    topPlot = new TH1F("yield_test_TTBar_All", "yield_test_TTBar_All", 12,0.,12.);
    bckgPlot = new TH1F("yield_test_BCKG_All", "yield_test_BCKG_All", 12,0.,12.);
    dataPlot = new TH1F("yield_test_Data_All", "yield_test_Data_All", 12,0.,12.);
    if (nonPrompt!=0) npPlot = new TH1F("yield_test_FakeRate_All", "yield_test_FakeRate_All", 12,0.,12.);
  }

  for ( int j=0;j<3;++j) {
    bin = signBin*3+j+offset+1;
    if (singleHisto) {
      sprintf(name,"yield_test_TTBar_All_%ix", bin-1);
      topPlot = new TH1F(name,name, 1, 0., 1.);
      sprintf(name,"yield_test_BCKG_All_%ix", bin-1);
      bckgPlot = new TH1F(name,name, 1, 0., 1.);
      bin=1;
    }
    float total=0.;
    float totalUnc=0.;
    for (unsigned int i=0;i<backgroundsTree.size();++i) {
      backgroundsTree[i]->setChannel(j);
      if ((sign==-1)&&(i==0)) {
	pair<int,float> res = backgroundsTree[i]->LoopDouble( );
	fillBin(topPlot, bin, res.second, res.first);
      }else {
        pair<int,float> res = backgroundsTree[i]->LoopDouble();
	if ((sign==-1)&& (offset==0)&&(i==(backgroundsTree.size()-1))) {
	  if (j==0) res = pair<int,float> (41,dyEstimate[0]);
	  if (j==2) res = pair<int,float> (18,dyEstimate[2]);
	}
	total+= res.second;
	if (res.first!=0) totalUnc+=res.second*res.second/res.first;
	if (backgroundsName[i]==string("TTZ")) {
	  totalUnc+=res.second*res.second*0.16;
	  cout << backgroundsName[i]<<" adding " << totalUnc<<" "<<res.second*0.4;
	  cout <<endl;
	}
      }
//        cout << backgroundsName[i]<<" "<<j<<" : "<<res.first<<" : "<<res.second<<endl;
    }
    cout << signBin*3+j+offset<<" "<<0<<" : "<<total<<endl;
    bckgPlot->SetBinContent(bin,total);
    bckgPlot->SetBinError(bin,sqrt(totalUnc));  
  }
  bckgPlot->Write();
  if (sign==0) topPlot->Write();

  if (nonPrompt!=0) {
    float totalFake[3];
    //nonPB(totalFake);
    for ( int j=0;j<3;++j) {
      bin = signBin*3+j+offset+1;
      FFPair fake = nonPB(j);
      totalFake[j]=fake.first;
      if (singleHisto) {
	sprintf(name,"yield_test_FakeRate_All_%ix", bin-1);
	npPlot = new TH1F(name,name, 1, 0., 1.);
	bin=1;
      }
      if (totalFake[j]<0) totalFake[j]=0;
      //npPlot->Fill(bin-1, totalFake[j]);
      npPlot->SetBinContent(bin,fake.first);
      npPlot->SetBinError(bin,fake.second);  
    }
  }

  for ( int j=0;j<3;++j) {
    bin = signBin*3+j+offset+1;
    if (singleHisto) {
      sprintf(name,"yield_test_Data_All_%ix", bin-1);
      dataPlot = new TH1F(name,name, 1, 0., 1.);
      bin=1;
    }
    dataTree[0]->setChannel(j);
    dataPlot->Fill(bin-1, dataTree[0]->Loop() );
  }
  outFile->Write();
  outFile->Close();

}

void redoPlot(int plotType, string sampleName)
{

//   float massSF[3] = {0.002,0.002,0.002};

  //for (int mass=500;mass<=1500;mass+=100) {
  for (TreeMapIt mapIt = signalTreeMap.begin(); mapIt != signalTreeMap.end(); ++mapIt) {
    cout << mapIt->first<<endl;
    int mass = mapIt->first;
    for (unsigned int i=0;i<mapIt->second.size();++i) {
      signalTreeMap[mass][i]->plots(plotType, sampleName+"_"+signalName[mass][i] + getTagString());
      cout << sampleName+"_"+signalName[mass][i] + getTagString()<<endl;
//       signalTreeMap[mass][i]->setSF(massSF);
      for ( int j=-1;j<3;++j) {
	signalTreeMap[mass][i]->setChannel(j);
	signalTreeMap[mass][i]->Loop();
	signalTreeMap[mass][i]->newPlot->Write();
      }
      signalTreeMap[mass][i]->noPlots();
    }
  }

  for (unsigned int i=0;i<backgroundsTree.size();++i) {
    backgroundsTree[i]->plots(plotType, sampleName+"_"+backgroundsName[i] + getTagString());
    for ( int j=-1;j<3;++j) {
      backgroundsTree[i]->setChannel(j);
      backgroundsTree[i]->Loop();
      backgroundsTree[i]->newPlot->Write();
    }
    backgroundsTree[i]->noPlots();
  }

  if (sign==1) {
    nonPrompt->plots(plotType, sampleName+"_FakeRate");
    nonPrompt->setChannel(0);
    TH1F* mm = nonPromptPlotSS(FR_MU, FR_MU, PR_MU, PR_MU);
    nonPrompt->setChannel(1);
    TH1F* em = nonPromptPlotSS(FR_EL, FR_MU, PR_EL, PR_MU);
    nonPrompt->setChannel(2);
    TH1F* ee = nonPromptPlotSS(FR_EL, FR_EL, PR_EL, PR_EL);
    string name = ee->GetName();

    name.replace(name.end()-4,name.end(),"All");
    cout << name<<endl;

    TH1F* all = (TH1F*)mm->Clone(name.c_str());
    all->Add(em);
    all->Add(ee);
    mm->Write();
    em->Write();
    ee->Write();
    all->Write();
    nonPrompt->noPlots();
  } else if (sign==0) {
    nonPrompt->plots(plotType, sampleName+"_FakeRate");
    nonPrompt->setChannel(0);
    TH1F* mm = nonPromptPlotTRI(FR_MU, PR_MU);
    nonPrompt->setChannel(1);
    TH1F* em = nonPromptPlotTRI(FR_MU, PR_EL);
    nonPrompt->setChannel(2);
    TH1F* ee = nonPromptPlotTRI(FR_EL, PR_EL);
    string name = ee->GetName();

    name.replace(name.end()-4,name.end(),"All");

    TH1F* all = (TH1F*)mm->Clone(name.c_str());
    all->Add(em);
    all->Add(ee);
    mm->Write();
    em->Write();
    ee->Write();
    all->Write();
    nonPrompt->noPlots();
  }

  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->plots(plotType, sampleName+"_Data" + getTagString());
    for ( int j=-1;j<3;++j) {
      dataTree[i]->setChannel(j);
      dataTree[i]->Loop();
      dataTree[i]->newPlot->Write();
    }
    dataTree[i]->noPlots();
  }


}

int MVAtree::countLepton(int pdgid)
{
  int i=0;
  if (abs(isMuon1)==pdgid ) ++i;
  if (abs(isMuon2)==pdgid ) ++i;
  if (abs(isMuon3)==pdgid ) ++i;
  if (abs(isMuon4)==pdgid ) ++i;
  return i;
}

pair<float,float> getDYSFReMu()
{
//   gStyle->SetPadRightMargin(0.13);
//   gStyle->SetPadLeftMargin(0.13);
//   TCanvas *c1 = new TCanvas("c1", "plots",400,400,800,600);
//   c1->SetFillColor(10);
//   c1->SetFillStyle(4000);
//   c1->SetBorderSize(2); 
//   c1->SetLogy(0);

  string sampleName = "DYSF";
  for (unsigned int i=0;i<dataTree.size();++i) {
    dataTree[i]->plots(MLL, sampleName+"_Data");
    dataTree[i]->setChannel(1);
    dataTree[i]->Loop();
    dataTree[i]->noPlots();
  }

  TH1F *zMass_Data      = dataTree[0]->newPlot;
  zMass_Data->Draw();

 
  //*********************************************
  // get the histograms and prepare the sum of MC
  //*********************************************

  TH1F * totMC = new TH1F("totMC", "totMC", zMass_Data->GetNbinsX() , zMass_Data->GetXaxis()->GetXmin() , zMass_Data->GetXaxis()->GetXmax());
  for (unsigned int i=0;i<backgroundsTree.size()-1;++i) {
    backgroundsTree[i]->plots(MLL, sampleName+"_"+backgroundsName[i]);
    backgroundsTree[i]->setChannel(1);
    backgroundsTree[i]->Loop();
    backgroundsTree[i]->newPlot->Sumw2();
    totMC->Add( totMC, backgroundsTree[i]->newPlot, 1, 1);
    backgroundsTree[i]->noPlots();
  }
  cout << "a\n";

  //*********************
  // summ up MC templates
  //*********************
  backgroundsTree.back()->plots(MLL, sampleName+"_"+backgroundsName.back());
  backgroundsTree.back()->setChannel(1);
  backgroundsTree.back()->Loop();
  backgroundsTree.back()->newPlot->Sumw2();
  TH1F * totMC_DY = backgroundsTree.back()->newPlot;
  backgroundsTree.back()->noPlots();

    cout << "a\n";

  zMass_Data->SetMarkerStyle(20);
  zMass_Data->Draw("ep");
  totMC->SetLineColor(2);
  totMC->SetLineWidth(1.2);
  totMC->Draw("hsame");
  totMC_DY->SetLineColor(4);
  totMC_DY->SetLineWidth(1.2);
  totMC_DY->Draw("hsame");
  zMass_Data->Draw("epsame");
  
  
  
  //*********************************************
  // create the RooFit environement
  RooWorkspace *w = new RooWorkspace("w",kTRUE) ;  

  //*********************************************
  // create a RooFit variable : 
  // the variable that will be used for the fit
  //*********************************************
  RooRealVar * LL_Zmass        = new RooRealVar("LL_Zmass", "M_{e#mu}", zMass_Data->GetXaxis()->GetXmin() , zMass_Data->GetXaxis()->GetXmax());
  
  //create RooDataHisto for the data (to be fitted)
  RooDataHist* histoLL_Zmass   = new RooDataHist("histoFit_LL", "histoFit_LL",  *LL_Zmass,  zMass_Data ); 
  //create RooDataHisto for the MC DY : will be used to create template 1
  RooDataHist* histoFit_Template_DYEM        = new RooDataHist("histoFit_Template_DYEM",       "histoFit_Template_DYEM", *LL_Zmass, totMC_DY);
  //create RooDataHisto for the MC DY : will be used to create template 2
  RooDataHist* histoFit_Template_OtherLL     = new RooDataHist("histoFit_Template_OtherLL",    "histoFit_Template_OtherLL", *LL_Zmass, totMC);        

  //convert  RooDataHisto into the template 1
  RooHistPdf*  histoFit_Template_DYEM_pdf    = new RooHistPdf("histoFit_Template_DYEM_pdf",    "histoFit_Template_DYEM_pdf",    *LL_Zmass,  *histoFit_Template_DYEM);
  //convert  RooDataHisto into the template 2
  RooHistPdf*  histoFit_Template_OtherLL_pdf    = new RooHistPdf("histoFit_Template_OtherLL_pdf",    "histoFit_Template_OtherLL_pdf",    *LL_Zmass,  *histoFit_Template_OtherLL);

  //define a coefficient : it is the output of the fit
  //it represents the contribution (fraction) of one of the 2 template with resepct to the data after the fit :
  // N(event ttbar) = coeff*N(event data)
  RooRealVar coeffModel("coeffModel", "coeffModel", 0.5, 0., 1.);
  //associate the templates, the data histo and the coeff.
  RooAddPdf *pdf_sum = new RooAddPdf("pdf_sum"," test pdf sum",RooArgSet(*histoFit_Template_DYEM_pdf, *histoFit_Template_OtherLL_pdf), coeffModel);
  //do the fit.
  RooFitResult* myFitResults_all = pdf_sum->fitTo(*histoLL_Zmass, RooFit::Save()) ;

//   //create a "frame" : a kind of TCanvas used to display the result of the fit
//   RooPlot* frame = LL_Zmass->frame() ;
//   //plot the data on the frame
//   histoLL_Zmass->plotOn(frame) ;
//   //plots the templates (after the fit) on the frame
//   pdf_sum->plotOn(frame, Components(*histoFit_Template_DYEM_pdf), VisualizeError(*myFitResults_all), FillColor(kGreen), LineWidth(2)  );
//   pdf_sum->plotOn(frame, Components(*histoFit_Template_OtherLL_pdf), VisualizeError(*myFitResults_all), FillColor(kRed),   LineWidth(2)  );
//   pdf_sum->plotOn(frame, LineStyle(kDashed),                VisualizeError(*myFitResults_all), FillColor(kBlue), LineWidth(2) );
//   histoLL_Zmass->plotOn(frame) ;
  //histoFit_all_loose->plotOn(frame_loose) ;
  
  //c1->SetLogy();
  
  //frame->Minimum(0.01);
  //draw the frame
//   frame->Draw() ;
  //print the coeff after the fit
  coeffModel.Print() ;
  
  //calculate the various contributions
  cout << " the DY is " << coeffModel.getVal()*zMass_Data->Integral() << " +/- " << coeffModel.getError()*zMass_Data->Integral() << endl;
  cout << " the other " << (1-coeffModel.getVal())*zMass_Data->Integral() << " +/- " << coeffModel.getError()*zMass_Data->Integral() << endl;
  
  cout << " initial DY component " << totMC_DY->Integral() << " +/-" <<  totMC_DY->Integral() << endl;
  
  cout << "scale factor " << coeffModel.getVal()*zMass_Data->Integral()/totMC_DY->Integral() << "pm " <<
  coeffModel.getError()*zMass_Data->Integral()/totMC_DY->Integral()
   << endl;
  
     
//  TH1F * histoFitDY    = new TH1F("histoFitDY",    "histoFitDY",    0, 1, 100);
//  TH1F * histoFitOther = new TH1F("histoFitOther", "histoFitOther", 0, 1, 100);
//  TH1F * histoFitTot = new TH1F("histoFitOther", "histoFitOther", 0, 1, 100);
//   
//   histoFitDY->SetFillColor(3);
//   histoFitOther->SetFillColor(2);
//   histoFitTot->SetFillColor(4);
//   TLegend* qw = new TLegend(0.75,0.70,0.98,0.98);
//   qw->AddEntry(zMass_Data,         "Data" ,                "p");
//   qw->AddEntry(histoFitTot,    "result of the Fit" ,    "f");
//   qw->AddEntry(histoFitDY,     "DY component" ,        "f");
//   qw->AddEntry(histoFitOther,  "non DY component" ,    "f");
//   
//   
//   qw->Draw();

  return FFPair(coeffModel.getVal()*zMass_Data->Integral(), coeffModel.getVal()*zMass_Data->Integral()/totMC_DY->Integral());
}

void redoAllPlots(string sampleName, string rootFileName)
{
  TH1::AddDirectory(kFALSE);
  TFile* outFile = new TFile(rootFileName.c_str(), "RECREATE");
//   redoPlot(AMWT,sampleName);
//   redoPlot(HT,sampleName);
//   redoPlot(MINMLB,sampleName);
//   redoPlot(ST,sampleName);
//   redoPlot(NBJET,sampleName); // nBjets
//   redoPlot(RCN,sampleName); // nBjets
//   redoPlot(NEHAD,sampleName); // nBjets
//   redoPlot(CHHAD,sampleName); // nBjets
//   redoPlot(JETPT,sampleName);
//   redoPlot(NJET,sampleName); // nBjets
//   redoPlot(MLL,sampleName);
 redoPlot(AMWT,sampleName);
  
  outFile->Write();
  outFile->Close();
}
int main( int argc, const char* argv[] ){

  sign = -1;
//   if (strcmp(argv[1],"OS")==0) {
//     sign = -1;
//   } else if (strcmp(argv[1],"SS")==0) {
//     sign = +1;
//   } else if (strcmp(argv[1],"TRI")==0) {
//     sign = 0;
//   } else {
//     cout<<"Lepton number/sign must be (OS, LS, TRI)"<<endl;
//     return 0;
//   }
  TString signName("OS");
  TString dirName(argv[1]);
  cout <<argc<<endl;
  int minBJet=0, maxBJet=20;
  if (argc >2) minBJet=atoi(argv[2]);
  if (argc >3) maxBJet=atoi(argv[3]);


  bool doFakes = false;
  nonPrompt=0;
  

// For OS
  if (sign==-1) {
     int mass;
//      mass = 171; signalTreeMap[mass].push_back(new MVAtree(dirName+"/Top172_JEC_Total_down.root"));signalName[mass].push_back("TTbar_171.5");signalNorm[mass].push_back(1.);signalBR[mass].push_back(1);
//      mass = 173; signalTreeMap[mass].push_back(new MVAtree(dirName+"/Top172_JEC_Total_up.root"));signalName[mass].push_back("TTbar_173.5");signalNorm[mass].push_back(1.);signalBR[mass].push_back(1);
//      mass = 169; signalTreeMap[mass].push_back(new MVAtree(dirName+"/Top172_JEC_AbsoluteStat_down.root"));signalName[mass].push_back("TTbar_169.5");signalNorm[mass].push_back(1.);signalBR[mass].push_back(1);
//      mass = 175; signalTreeMap[mass].push_back(new MVAtree(dirName+"/Top172_JEC_AbsoluteStat_up.root"));signalName[mass].push_back("TTbar_175.5");signalNorm[mass].push_back(1.);signalBR[mass].push_back(1);

//      mass = 166; signalTreeMap[mass].push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top166.root"));signalName[mass].push_back("TTbar_166.5");signalNorm[mass].push_back(1.);signalBR[mass].push_back(1);
//      mass = 169; signalTreeMap[mass].push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top169.root"));signalName[mass].push_back("TTbar_169.5");signalNorm[mass].push_back(1.);signalBR[mass].push_back(1);
//      mass = 171; signalTreeMap[mass].push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top171.root"));signalName[mass].push_back("TTbar_171.5");signalNorm[mass].push_back(1.);signalBR[mass].push_back(1);
//      mass = 173; signalTreeMap[mass].push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top173.root"));signalName[mass].push_back("TTbar_173.5");signalNorm[mass].push_back(1.);signalBR[mass].push_back(1);
//      mass = 175; signalTreeMap[mass].push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top175.root"));signalName[mass].push_back("TTbar_175.5");signalNorm[mass].push_back(1.);signalBR[mass].push_back(1);
//      mass = 178; signalTreeMap[mass].push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top178.root"));signalName[mass].push_back("TTbar_178.5");signalNorm[mass].push_back(1.);signalBR[mass].push_back(1);

//      backgroundsTree.push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top166.root"));backgroundsName.push_back("TTbar166.5");backgroundsUnc.push_back(0.08);
//      backgroundsTree.push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top169.root"));backgroundsName.push_back("TTbar169.5");backgroundsUnc.push_back(0.08);
//      backgroundsTree.push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top171.root"));backgroundsName.push_back("TTbar171.5");backgroundsUnc.push_back(0.08);
//   TChain *ttch = new TChain("ST","ST");
//   ttch->AddFile(dirName+"/ptOutFile_MuMu_"+signName+"_2012_TTbarLept.root", TChain::kBigNumber , "MVA") ;
//   ttch->AddFile(dirName+"/ptOutFile_MuMu_"+signName+"_2012_TTbarSL.root", TChain::kBigNumber , "MVA") ;
//   backgroundsTree.push_back(new MVAtree(ttch));backgroundsName.push_back("TTbar172.5");backgroundsUnc.push_back(0.20);
     backgroundsTree.push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top172.root"));backgroundsName.push_back("TTbar_172.5");backgroundsUnc.push_back(0.08);
//      backgroundsTree.push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top173.root"));backgroundsName.push_back("TTbar173.5");backgroundsUnc.push_back(0.08);
//      backgroundsTree.push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top175.root"));backgroundsName.push_back("TTbar175.5");backgroundsUnc.push_back(0.08);
//      backgroundsTree.push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Top178.root"));backgroundsName.push_back("TTbar178.5");backgroundsUnc.push_back(0.08);

    backgroundsTree.push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_WJets.root"));backgroundsName.push_back("WJets");backgroundsUnc.push_back(0.05);
    TChain *STch = new TChain("ST","ST");
    STch->AddFile(dirName+"/ptOutFile_MuMu_"+signName+"_2012_T_tW.root", TChain::kBigNumber , "MVA") ;
    STch->AddFile(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Tbar_tW.root", TChain::kBigNumber , "MVA") ;
    backgroundsTree.push_back(new MVAtree(STch));backgroundsName.push_back("SingleTop");backgroundsUnc.push_back(0.20);
    backgroundsTree.push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_WW.root"));backgroundsName.push_back("WW");backgroundsUnc.push_back(0.25);
    backgroundsTree.push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_WZ.root"));backgroundsName.push_back("WZ");backgroundsUnc.push_back(0.25);
    backgroundsTree.push_back(new MVAtree(dirName+"/ptOutFile_MuMu_"+signName+"_2012_ZZ.root"));backgroundsName.push_back("ZZ");backgroundsUnc.push_back(0.25);

    TChain *DYch = new TChain("DY","DY");
    DYch->AddFile(dirName+"/ptOutFile_MuMu_"+signName+"_2012_ZJets.root", TChain::kBigNumber , "MVA") ;
    DYch->AddFile(dirName+"/ptOutFile_MuMu_"+signName+"_2012_DY1050.root", TChain::kBigNumber , "MVA") ;
    backgroundsTree.push_back(new MVAtree(DYch));backgroundsName.push_back("DrellYan");backgroundsUnc.push_back(0.30);

  }


// Data
  TChain *ch = new TChain("aa","aa");
  ch->AddFile(dirName+"/ptOutFile_ElEl_"+signName+"_2012_Data.root", TChain::kBigNumber , "MVA") ;
  ch->AddFile(dirName+"/ptOutFile_ElMu_"+signName+"_2012_Data.root", TChain::kBigNumber , "MVA") ;
  ch->AddFile(dirName+"/ptOutFile_MuMu_"+signName+"_2012_Data.root", TChain::kBigNumber , "MVA") ;

  dataTree.push_back(new MVAtree(ch));dataName.push_back("Data");

// Now do the analysis:

  float metCut = 40;

     setBJets(minBJet, maxBJet);
     setMETcut(metCut, 99999., false);
     setLeptonSumMass(76,106, true);

//   redoAllPlots("test", "plots_AllJet.root");
//     yieldTable(500,doFakes);
// exit(0);
// //   yieldTable(500,doFakes);
//      setBJets(minBJet, maxBJet);
//   redoAllPlots("test", "plots.root");
//   yieldTable(500,doFakes);
//   exit(1);
//      setMETcut(metCut, 99999., true);
//   yieldTable(500,doFakes);
//      setLeptonSumMass(76,106, true);
//   yieldTable(500,doFakes);

     setBJets(minBJet, maxBJet);
     setMETcut(metCut, 99999., false);
     setLeptonSumMass(76,106, true);

//  redoAllPlots("test", "plots.root");
//   yieldTable(500,doFakes);
//   npSF[1] = 1.18507;//getDYSFReMu().second;
//   backgroundsTree.back()->setSF(npSF);
//   yieldTable(500,doFakes);
// exit(1);
     npSF[0] = getDYSFRoutinSlow(0).second;
//      npSF[0] = 1.1531;//getDYSFRoutin(0).second;
//      npSF[0] = getDYSFRoutin(0).second;
     setBJets(minBJet, maxBJet);
     setMETcut(metCut, 99999., false);
     setLeptonSumMass(76,106, true);
     npSF[2] = getDYSFRoutin(2).second;
//      npSF[2] = 1.38909 ;// getDYSFRoutin(2).second;
// npSF[1] = 1.572;
     backgroundsTree.back()->setSF(npSF);


     setBJets(minBJet, maxBJet);
     setMETcut(metCut, 99999., false);
     setLeptonSumMass(76,106);
   yieldTable(500,doFakes);
//    redoAllPlots("test", "plots.root");
     exit(0);

//   cutFlowTable(-1, false, false, true);





}


