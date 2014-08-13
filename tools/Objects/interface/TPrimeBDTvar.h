#ifndef TPrimeBDTvar_cxx
#define TPrimeBDTvar_cxx

#include <TTree.h>
#include "../interface/TprimeEvent.h"
#include "../interface/treetop.h"

class BDTVar{

public:
  BDTVar(const double scale);	
  int isMuon1;
  int isMuon2;
  int isMuon3;
  int isMuon4;
  float PF_met_pt ;
  int nBTag;
  int nJets;
  int nLeptons;
  float lepton1_eta;
  float lepton2_eta;
  float lepton3_eta;
  float lepton4_eta;
  float lepton1_pt;
  float lepton2_pt;
  float lepton3_pt;
  float lepton4_pt;
  float lepton1_phi;
  float lepton2_phi;
  float lepton3_phi;
  float lepton4_phi;
  float lepton1_energy;
  float lepton2_energy;
  float lepton3_energy;
  float lepton4_energy;
  float jet1_pt,    jet5_pt;
  float jet2_pt,    jet6_pt;
  float jet3_pt,    jet7_pt;
  float jet4_pt,    jet8_pt;
  float jet1_eta,   jet5_eta;
  float jet2_eta,   jet6_eta;
  float jet3_eta,   jet7_eta;
  float jet4_eta,   jet8_eta;
  float jet1_phi,   jet5_phi;
  float jet2_phi,   jet6_phi;
  float jet3_phi,   jet7_phi;
  float jet4_phi,   jet8_phi;
  float jet1_energy,jet5_energy;
  float jet2_energy,jet6_energy;
  float jet3_energy,jet7_energy;
  float jet4_energy,jet8_energy;
  float minDRJetLepton;
  float dRJetLepton;
  float HT;
  float minMLB, maxMLB;
  float minDRBJetLepton;
  float leptonSumMass;
  float CAWBtag_disc1;
  float CAWBtag_disc2;
  float CATopJetsBtag_disc1;
  float CATopJetsBtag_disc2;
  float CATopJets1_pt;
  float CATopJets2_pt;
  float CATopJets1_eta;
  float CATopJets2_eta;
  float CATopJets1_phi;
  float CATopJets2_phi;
  float CATopJets1_energy;
  float CATopJets2_energy;
  float CAWJets1_pt;
  float CAWJets2_pt;
  float CAWJets1_eta;
  float CAWJets2_eta;
  float CAWJets1_phi;
  float CAWJets2_phi;
  float CAWJets1_energy;
  float CAWJets2_energy;
  int nCAWJets;
  int nCATopJets;
  float amwt;
  float amwtPeakWeight;
  float weight;
  float leptonSum, jetSum, leptonMETSum, leptonJetsSum, jetsMETSum, leptonJetsMETSum;
  int eventType;
  float isolation1, isolation2, isolation3, isolation4;
  int jets_btag;
  int run, lumi, event;
  int jet_nChHad[8];
  int jet_nNeHad[8];
  //,jet2_nChHad, jet2_nNeHad, jet3_nChHad, jet3_nNeHad, jet4_nChHad, jet4_nNeHad, jet5_nChHad, jet5_nNeHad, jet6_nChHad, jet6_nNeHad, jet7_nChHad, jet7_nNeHad, jet8_nChHad, jet8_nNeHad;
  float jet_RCN[8];//, jet2_RCN, jet3_RCN, jet4_RCN, jet5_RCN, jet6_RCN, jet7_RCN, jet8_RCN;
  std::vector<double> pdf;

  float CATopJets1_jetMass, CATopJets2_jetMass, CATopJets1_massDrop, CATopJets2_massDrop, CATopJets1_minMass, CATopJets2_minMass;
  int CATopJets1_nbrConst, CATopJets2_nbrConst;
  float CAWJets1_jetMass, CAWJets2_jetMass, CAWJets1_massDrop, CAWJets2_massDrop, CAWJets1_minMass, CAWJets2_minMass;
  int CAWJets1_nbrConst, CAWJets2_nbrConst;

   Double_t Q, x1, x2, pdf1, pdf2;
   int id1, id2;
  
  TTree* MVA;
  void FillBranches(TprimeEvent* teve, treetop* theTree, TString type="", bool doFakes=false);
  void WriteBranches();

  double scale_;
 };


#endif
