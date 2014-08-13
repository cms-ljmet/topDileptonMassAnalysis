//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 13 03:17:52 2012 by ROOT version 5.27/06b
// from TTree treetop/treetop
// found on file: TTbar_Aug13_all.root
//////////////////////////////////////////////////////////

#ifndef treetop_h
#define treetop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include "TLepton.h"
#include <stdio.h>


class treetop {
 public :
   struct SortByPt {
       bool operator()(const TLepton *e1, const TLepton *e2) {
	   return (e1->lv.Pt() > e2->lv.Pt());
       }
   };
   struct SortJetsByPt {
       bool operator()(const TJet *e1, const TJet *e2) {
           return (e1->lv.Pt() > e2->lv.Pt());
       }
   };


  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  vector <TLepton*> allMuons, allElecs, goodMuons, goodElecs, looseMuons, looseElecs, lntMuons, lntElecs;
  vector <TJet*> allAK5Jets, allJets; //allJets is obsolete! It's really all AK5 jets
  // The goodLeptons and leptonsForFakes are ordered in pT and overlaps are cleaned (dR<0.1)
  vector <TLepton*> goodLeptons, leptonsForFakes, looseLeptons;
  
  vector <TJet*> goodJets(const vector <TLepton*> & leptons); //goodJets is really good AK5 jets without considering CA jets

  void setLeptonPtCut(double cut ){leptonPtCut_=cut;}
  void setJetPtCut(double cut ) {jetPtCut_=cut;}



//
  Int_t           run, lumi, event;
  
  Int_t           nPV;
  Int_t           nAK5Jet;
  Int_t           nPFjet;
  Int_t           nElectrons;
  Int_t           nMuons;

  Int_t           event_flavor; // Flavour of the leading jet
  Int_t           n_Bjets, n_Cjets, n_Lightjets;
  
// Pile-up reweighting: 

  Double_t        weight_PU, weight_PU_ABC;

// MET
  Double_t        PF_met_phi;
  Double_t        PF_met_pt;
  Double_t        PF_met_px;
  Double_t        PF_met_py;
  Double_t        corr_met, corr_met_phi;
  Double_t        higgsWeight;
  bool isHiggs, withPDF, withCHN;

   // Declaration of leaf types
   Int_t           bunchXing;
   Int_t           dataEE;
   Int_t           dataEM;
   Int_t           dataMM;
   Int_t           event_CommonCalc;
   Int_t           nAllJets_CommonCalc;
   Int_t           nInteractions;
   Double_t	   nTrueInteractions;
   Int_t           nLooseMuons_CommonCalc;
   Int_t           nSelBtagJets_CommonCalc;
   Int_t           nSelElectrons_CommonCalc;
   Int_t           nSelJets_CommonCalc;
   Int_t           nTightMuons_CommonCalc;
   Int_t           run_CommonCalc;
   Int_t           trigEE;
   Int_t           trigEM;
   Int_t           trigMM;
   vector<int>     *AK5JetTBag;
   vector<int>     *elChargeConsistent;
   vector<int>     *elCharge;
   vector<int>     *elIsEBEE;
   vector<int>     *elMHits;
   vector<int>     *elNotConversion;
   vector<int>     *elQuality;
   vector<int>     *elVtxFitConv;
   vector<int>     *genID;
   vector<int>     *genIndex;
   vector<int>     *genMotherID;
   vector<int>     *genMotherIndex;
   vector<int>     *genStatus;
   vector<int>     *muCharge;
   vector<int>     *muGlobal;
   vector<int>     *muNMatchedStations;
   vector<int>     *muNTrackerLayers;
   vector<int>     *muNValMuHits;
   vector<int>     *muNValPixelHits;
   vector<int>     *muQuality;
   vector<double>  *AK5JetEnergy;
   vector<double>  *AK5JetEta;
   vector<double>  *AK5JetPhi;
   vector<double>  *AK5JetPt;
   vector<double>  *elD0;
   vector<double>  *elDZ;
   vector<double>  *elDeta;
   vector<double>  *elDphi;
   vector<double>  *elDxy;
   vector<double>  *elEnergy;
   vector<double>  *elEta;
   vector<double>  *elHoE;
   vector<double>  *elOoemoop;
   vector<double>  *elPhi;
   vector<double>  *elPt;
   vector<double>  *elRelIso;
   vector<double>  *elSihih;
   vector<double>  *genEnergy;
   vector<double>  *genEta;
   vector<double>  *genPhi;
   vector<double>  *genPt;
   vector<double>  *muChi2;
   vector<double>  *muDxy;
   vector<double>  *muDz;
   vector<double>  *muEnergy;
   vector<double>  *muEta;
   vector<double>  *muPhi;
   vector<double>  *muPt;
   vector<double>  *muRelIso;

   vector<double>  *elAEff;
   vector<double>  *elChIso;
   vector<double>  *elNhIso;
   vector<double>  *elPhIso;
   vector<double>  *elRhoIso;
   vector<double>  *muChIso;
   vector<double>  *muGIso;
   vector<double>  *muNhIso;
   vector<double>  *muPuIso;
   vector<double>  *AK5JetRCN;
   vector<int>     *AK5JetnChHad;
   vector<int>     *AK5JetnNeuHad;


   Double_t        PdfWeightAverage_MRST2006nnlo, PdfWeightAverage_NNPDF10, PdfWeightAverage_cteq66, PdfWeightMinus_MRST2006nnlo, PdfWeightMinus_NNPDF10, PdfWeightMinus_cteq66, PdfWeightPlus_MRST2006nnlo, PdfWeightPlus_NNPDF10, PdfWeightPlus_cteq66;
   vector<double>  *PdfWeightsVec_MRST2006nnlo, *PdfWeightsVec_NNPDF10, *PdfWeightsVec_cteq66;

   Double_t Q, x1, x2, pdf1, pdf2;
   int id1, id2;

   // List of branches
   TBranch        *b_bunchXing_PileUpCalc;   //!
   TBranch        *b_dataEE_DileptonCalc;   //!
   TBranch        *b_dataEM_DileptonCalc;   //!
   TBranch        *b_dataMM_DileptonCalc;   //!
   TBranch        *b_event_CommonCalc;   //!
   TBranch        *b_lumi_CommonCalc;   //!
   TBranch        *b_nAllJets_CommonCalc;   //!
   TBranch        *b_nInteractions_PileUpCalc;   //!
   TBranch        *b_nTrueInteractions_PileUpCalc;   //!
   TBranch        *b_nLooseMuons_CommonCalc;   //!
   TBranch        *b_nSelBtagJets_CommonCalc;   //!
   TBranch        *b_nSelElectrons_CommonCalc;   //!
   TBranch        *b_nSelJets_CommonCalc;   //!
   TBranch        *b_nTightMuons_CommonCalc;   //!
   TBranch        *b_run_CommonCalc;   //!
   TBranch        *b_trigEE_DileptonCalc;   //!
   TBranch        *b_trigEM_DileptonCalc;   //!
   TBranch        *b_trigMM_DileptonCalc;   //!
   TBranch        *b_corr_met_DileptonCalc;   //!
   TBranch        *b_corr_met_phi_DileptonCalc;   //!
   TBranch        *b_met_DileptonCalc;   //!
   TBranch        *b_met_phi_DileptonCalc;   //!
   TBranch        *b_weight_PU_PileUpCalc;   //!
   TBranch        *b_AK5JetTBag_DileptonCalc;   //!
   TBranch        *b_elChargeConsistent_DileptonCalc;   //!
   TBranch        *b_elCharge_DileptonCalc;   //!
   TBranch        *b_elIsEBEE_DileptonCalc;   //!
   TBranch        *b_elMHits_DileptonCalc;   //!
   TBranch        *b_elNotConversion_DileptonCalc;   //!
   TBranch        *b_elQuality_DileptonCalc;   //!
   TBranch        *b_elVtxFitConv_DileptonCalc;   //!
   TBranch        *b_genID_DileptonCalc;   //!
   TBranch        *b_genIndex_DileptonCalc;   //!
   TBranch        *b_genMotherID_DileptonCalc;   //!
   TBranch        *b_genMotherIndex_DileptonCalc;   //!
   TBranch        *b_genStatus_DileptonCalc;   //!
   TBranch        *b_muCharge_DileptonCalc;   //!
   TBranch        *b_muGlobal_DileptonCalc;   //!
   TBranch        *b_muNMatchedStations_DileptonCalc;   //!
   TBranch        *b_muNTrackerLayers_DileptonCalc;   //!
   TBranch        *b_muNValMuHits_DileptonCalc;   //!
   TBranch        *b_muNValPixelHits_DileptonCalc;   //!
   TBranch        *b_muQuality_DileptonCalc;   //!
   TBranch        *b_AK5JetEnergy_DileptonCalc;   //!
   TBranch        *b_AK5JetEta_DileptonCalc;   //!
   TBranch        *b_AK5JetPhi_DileptonCalc;   //!
   TBranch        *b_AK5JetPt_DileptonCalc;   //!
   TBranch        *b_elD0_DileptonCalc;   //!
   TBranch        *b_elDZ_DileptonCalc;   //!
   TBranch        *b_elDeta_DileptonCalc;   //!
   TBranch        *b_elDphi_DileptonCalc;   //!
   TBranch        *b_elDxy_DileptonCalc;   //!
   TBranch        *b_elEnergy_DileptonCalc;   //!
   TBranch        *b_elEta_DileptonCalc;   //!
   TBranch        *b_elHoE_DileptonCalc;   //!
   TBranch        *b_elOoemoop_DileptonCalc;   //!
   TBranch        *b_elPhi_DileptonCalc;   //!
   TBranch        *b_elPt_DileptonCalc;   //!
   TBranch        *b_elRelIso_DileptonCalc;   //!
   TBranch        *b_elSihih_DileptonCalc;   //!
   TBranch        *b_genEnergy_DileptonCalc;   //!
   TBranch        *b_genEta_DileptonCalc;   //!
   TBranch        *b_genPhi_DileptonCalc;   //!
   TBranch        *b_genPt_DileptonCalc;   //!
   TBranch        *b_muChi2_DileptonCalc;   //!
   TBranch        *b_muDxy_DileptonCalc;   //!
   TBranch        *b_muDz_DileptonCalc;   //!
   TBranch        *b_muEnergy_DileptonCalc;   //!
   TBranch        *b_muEta_DileptonCalc;   //!
   TBranch        *b_muPhi_DileptonCalc;   //!
   TBranch        *b_muPt_DileptonCalc;   //!
   TBranch        *b_muRelIso_DileptonCalc;   //!

   TBranch        *b_nPV_DileptonCalc;   //!
   TBranch        *b_higgsWeight_DileptonCalc;   //!
   TBranch        *b_PdfWeightAverage_MRST2006nnlo_PdfCalc, *b_PdfWeightAverage_NNPDF10_PdfCalc, *b_PdfWeightAverage_cteq66_PdfCalc, *b_PdfWeightMinus_MRST2006nnlo_PdfCalc, *b_PdfWeightMinus_NNPDF10_PdfCalc, *b_PdfWeightMinus_cteq66_PdfCalc, *b_PdfWeightPlus_MRST2006nnlo_PdfCalc, *b_PdfWeightPlus_NNPDF10_PdfCalc, *b_PdfWeightPlus_cteq66_PdfCalc;   //!
   TBranch        *b_PdfWeightsVec_MRST2006nnlo_PdfCalc, *b_PdfWeightsVec_NNPDF10_PdfCalc, *b_PdfWeightsVec_cteq66_PdfCalc;   //!
   TBranch        *b_elAEff_DileptonCalc;   //!
   TBranch        *b_elChIso_DileptonCalc;   //!
   TBranch        *b_elNhIso_DileptonCalc;   //!
   TBranch        *b_elRhoIso_DileptonCalc;   //!
   TBranch        *b_elPhIso_DileptonCalc;   //!
   TBranch        *b_muChIso_DileptonCalc;   //!
   TBranch        *b_muGIso_DileptonCalc;   //!
   TBranch        *b_muNhIso_DileptonCalc;   //!
   TBranch        *b_muPuIso_DileptonCalc;   //!
   TBranch        *b_AK5JetnChHad_DileptonCalc;   //!
   TBranch        *b_AK5JetnNeuHad_DileptonCalc;   //!
   TBranch        *b_AK5JetRCN_DileptonCalc;   //!
   TBranch *b_Q_PdfCalc, *b_id1_PdfCalc, *b_id2_PdfCalc, *b_x1_PdfCalc, *b_x2_PdfCalc, *b_pdf1_PdfCalc, *b_pdf2_PdfCalc;

  treetop(TTree *tree);
  treetop(const TString &filename);
  virtual ~treetop();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  double getEntries();
  void printEvent(int jentry);
  void printEvent();
  void checkLeptons();
  double entries;
  double leptonPtCut_, jetPtCut_;
};

#endif

#ifdef treetop_cxx

#endif // #ifdef treetop_cxx
