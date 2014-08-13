#define treetop_cxx
#include "../interface/treetop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>

void treetop::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L treetop.C
//      Root > treetop t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

//    Long64_t nentries = fChain->GetEntriesFast();
  TH1F  *h = new TH1F("hist", "his", 50, -0.5, 49.5);

//    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<1;jentry++) {
//    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      h->Fill(nTrueInteractions);
      if (ientry < 0) break;

   }
}

double treetop::getEntries()
{
  if (entries<0) entries = fChain->GetEntriesFast();
  return entries;
}; 


void treetop::printEvent(int jentry)
{
  GetEntry(jentry);
  printEvent();
}

void treetop::printEvent()
{
  cout << "run / ls / event = " << run << " / "<< lumi<< " / "<< event<< endl;
  cout << "nElectrons = " << nElectrons << endl;
  for (unsigned int i = 0;i<allElecs.size();++i ) allElecs[i]->printLepton();
  cout << "nMuons = " << nMuons << endl;  
  for (unsigned int i = 0;i<allMuons.size();++i ) allMuons[i]->printLepton();
  cout << "nAK5Jet = " << allAK5Jets.size() << endl;  
  for (unsigned int i = 0;i<allAK5Jets.size();++i ) allAK5Jets[i]->printJet();
cout <<"=======================\n";
}

void treetop::checkLeptons()
{
  goodMuons.clear();
  goodElecs.clear();
  looseMuons.clear();
  looseElecs.clear();
  goodLeptons.clear();
  lntMuons.clear();
  lntElecs.clear();
  leptonsForFakes.clear();

  for (unsigned int i = 0;i<allElecs.size();++i ) {
// cout << "Elec: "<<allElecs[i]->isDileptonLepton(leptonPtCut_)<<" "<<allElecs[i]->isLeptJetLepton(leptonPtCut_)<<" "<< allElecs[i]->goodElectron(leptonPtCut_)<<" "<<
// allElecs[i]->looseElectron(leptonPtCut_)<<endl;
// allElecs[i]->printLepton();
    if (allElecs[i]->isDileptonLepton(leptonPtCut_))   goodElecs.push_back(allElecs[i]);
    if (allElecs[i]->goodElectron(leptonPtCut_)) looseElecs.push_back(allElecs[i]);
    if (allElecs[i]->looseElectron(leptonPtCut_) && !allElecs[i]->goodElectron(leptonPtCut_)) lntElecs.push_back(allElecs[i]);
  }
  for (unsigned int i = 0;i<allMuons.size();++i ) {
// cout << "Muon: "<<allMuons[i]->isDileptonLepton(leptonPtCut_)<<" "<<allMuons[i]->isLeptJetLepton(leptonPtCut_)<<" "<< allMuons[i]->goodMuon(leptonPtCut_)<<" "<<
// allMuons[i]->looseMuon(leptonPtCut_)<<endl;
// allMuons[i]->printLepton();
    if (allMuons[i]->isDileptonLepton(leptonPtCut_))   goodMuons.push_back(allMuons[i]);
    if (allMuons[i]->goodMuon(leptonPtCut_)) looseMuons.push_back(allMuons[i]);
    if (allMuons[i]->looseMuon(leptonPtCut_) && !allMuons[i]->goodMuon(leptonPtCut_)) lntMuons.push_back(allMuons[i]);
  }
  vector <TLepton*> intLeptons;
  intLeptons.insert(intLeptons.end(), goodElecs.begin(), goodElecs.end());
  intLeptons.insert(intLeptons.end(), goodMuons.begin(), goodMuons.end());
  sort(intLeptons.begin(), intLeptons.end(), SortByPt());
  for (unsigned int i = 0;i<intLeptons.size();++i ) {
    bool pass=true;
    for (unsigned int j = 0;j<goodLeptons.size();++j ) 
      if (intLeptons[i]->lv.DeltaR(goodLeptons[j]->lv) < 0.1) pass=false; 
    if (pass) goodLeptons.push_back(intLeptons[i]);
    }

  //Repeat for leptons used in fakes study
  intLeptons.clear();
  intLeptons.insert(intLeptons.end(), goodLeptons.begin(), goodLeptons.end());
  intLeptons.insert(intLeptons.end(), lntElecs.begin(), lntElecs.end());
  intLeptons.insert(intLeptons.end(), lntMuons.begin(), lntMuons.end());

  sort(intLeptons.begin(), intLeptons.end(), SortByPt());
  
  for (unsigned int i = 0;i<intLeptons.size();++i ) {
    bool pass=true;
    for (unsigned int j = 0;j<leptonsForFakes.size();++j ) 
      if (intLeptons[i]->lv.DeltaR(leptonsForFakes[j]->lv) < 0.1) pass=false; 
    if (pass) leptonsForFakes.push_back(intLeptons[i]);
  }
  sort(leptonsForFakes.begin(), leptonsForFakes.end(), SortByPt());

}


//I keep this for compatibility, but it is confusing
vector <TJet*> treetop::goodJets(const vector <TLepton*> & leptons)
{
  vector <TJet*>  theJets;
  for (unsigned int i = 0;i<allJets.size();++i ) {
    bool answer = true;
    
    for (unsigned int j = 0;j<leptons.size();++j ) 
      {   
      if (! allJets[i]->checkJetLeptondR(*leptons[j])) answer = false;
      }
    if (answer) theJets.push_back(allJets[i]);
  }
  return theJets;

}

treetop::treetop(const TString &filename)
{
 cout << "Opening: "<<filename<<endl;
  TFile *f = new TFile(filename);
  f->cd();
  TTree *tree = (TTree*)gDirectory->Get("ljmet");
  Init(tree);
}

treetop::treetop(TTree *tree)
{
  Init(tree);
}
treetop::~treetop()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t treetop::GetEntry(Long64_t entry)
{
//Clear previous containers:
  for (unsigned int i = 0;i<allMuons.size();++i ) delete allMuons[i];
  for (unsigned int i = 0;i<allElecs.size();++i ) delete allElecs[i];
  for (unsigned int i = 0;i<allAK5Jets.size();++i ) delete allAK5Jets[i];
  allMuons.clear();
  allElecs.clear();
  allAK5Jets.clear();

  //I keep this for compatibility, but it is confusing and we should get rid of it
  for (unsigned int i = 0;i<allJets.size();++i ) delete allJets[i];
  allJets.clear();

  // Read contents of entry.
  if (!fChain) return 0;
  int stat =  fChain->GetEntry(entry);
  nMuons = muCharge->size();
  nElectrons = elPt->size();
  nAK5Jet = AK5JetTBag->size();
  nPFjet = AK5JetTBag->size(); //It's not really all PF jets; it's all AK5 jets without considering CAW

  for (int i = 0;i<nMuons;++i ) 
    allMuons.push_back(new TLepton((*muCharge)[i], (*muPt)[i],
      (*muEta)[i], (*muPhi)[i], (*muEnergy)[i],
      (*muRelIso)[i], (*muDxy)[i], (*muDz)[i],
      (*muGlobal)[i],(*muNTrackerLayers)[i],
      (*muChi2)[i], (*muNValMuHits)[i], (*muNValPixelHits)[i], (*muNMatchedStations)[i], (*muQuality)[i]));

  for (int i = 0;i<nElectrons;++i)
    allElecs.push_back(new TLepton((*elCharge)[i], (*elPt)[i],
      (*elEta)[i], (*elPhi)[i], (*elEnergy)[i],
      (*elRelIso)[i], (*elDxy)[i], (*elDZ)[i],
      (*elIsEBEE)[i], (*elQuality)[i], (*elChargeConsistent)[i],
      (*elNotConversion)[i], (*elDeta)[i], (*elDphi)[i], 
      (*elSihih)[i], (*elHoE)[i], (*elD0)[i], (*elOoemoop)[i],
      (*elMHits)[i], (*elNotConversion)[i]));

// 	cout << "--------------------------------------------\n";
  for (int i = 0;i<nAK5Jet;++i ) 
    allAK5Jets.push_back(new TJet((*AK5JetPt)[i],
      (*AK5JetEta)[i], (*AK5JetPhi)[i], (*AK5JetEnergy)[i], (bool) (*AK5JetTBag)[i],
      (*AK5JetnChHad)[i], (*AK5JetnNeuHad)[i], (*AK5JetRCN)[i]));
  //         for (int i = 0;i<nAK5Jet;++i ) allAK5Jets[i]->printJet();
// 	cout << "--------------------------------------------\n";


  //I keep this for compatibility, but it is confusing and we should get rid of it
  for (int i = 0;i<nAK5Jet;++i ) 
    allJets.push_back(new TJet((*AK5JetPt)[i],
      (*AK5JetEta)[i], (*AK5JetPhi)[i], (*AK5JetEnergy)[i], (bool) (*AK5JetTBag)[i],
      (*AK5JetnChHad)[i], (*AK5JetnNeuHad)[i], (*AK5JetRCN)[i]));


  PF_met_px = PF_met_pt*cos(PF_met_phi);
  PF_met_py = PF_met_pt*sin(PF_met_phi);

  checkLeptons();
  return stat;
}

Long64_t treetop::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void treetop::Init(TTree *tree)
{

  leptonPtCut_ = 0.;
  jetPtCut_ = 0.;

  if (!tree) return;

  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  entries = -1;

   AK5JetTBag = 0;
   elChargeConsistent = 0;
   elCharge = 0;
   elIsEBEE = 0;
   elMHits = 0;
   elNotConversion = 0;
   elQuality = 0;
   elVtxFitConv = 0;
   genID = 0;
   genIndex = 0;
   genMotherID = 0;
   genMotherIndex = 0;
   genStatus = 0;
   muCharge = 0;
   muGlobal = 0;
   muNMatchedStations = 0;
   muNTrackerLayers = 0;
   muNValMuHits = 0;
   muNValPixelHits = 0;
   muQuality = 0;
   AK5JetEnergy = 0;
   AK5JetEta = 0;
   AK5JetPhi = 0;
   AK5JetPt = 0;
   elD0 = 0;
   elDZ = 0;
   elDeta = 0;
   elDphi = 0;
   elDxy = 0;
   elEnergy = 0;
   elEta = 0;
   elHoE = 0;
   elOoemoop = 0;
   elPhi = 0;
   elPt = 0;
   elRelIso = 0;
   elSihih = 0;
   genEnergy = 0;
   genEta = 0;
   genPhi = 0;
   genPt = 0;
   muChi2 = 0;
   muDxy = 0;
   muDz = 0;
   muEnergy = 0;
   muEta = 0;
   muPhi = 0;
   muPt = 0;
   muRelIso = 0;

   elAEff = 0;
   elChIso = 0;
   elNhIso = 0;
   elPhIso = 0;
   elRhoIso = 0;
   muChIso = 0;
   muGIso = 0;
   muNhIso = 0;
   muPuIso = 0;
   withPDF=false;
   withCHN=false;
   AK5JetRCN = 0;
   AK5JetnChHad = 0;
   AK5JetnNeuHad = 0;


   PdfWeightsVec_MRST2006nnlo = 0;
   PdfWeightsVec_NNPDF10 = 0;
   PdfWeightsVec_cteq66 = 0;

//    // Set branch addresses and branch pointers
//    if (!tree) return;
//    fChain = tree;
//    fCurrent = -1;
//    fChain->SetMakeClass(1);
// 

   fChain->SetBranchAddress("run_CommonCalc", &run, &b_run_CommonCalc);
   fChain->SetBranchAddress("lumi_CommonCalc", &lumi, &b_lumi_CommonCalc);
   fChain->SetBranchAddress("event_CommonCalc", &event, &b_event_CommonCalc);

   fChain->SetBranchAddress("bunchXing_PileUpCalc", &bunchXing, &b_bunchXing_PileUpCalc);
   fChain->SetBranchAddress("nInteractions_PileUpCalc", &nInteractions, &b_nInteractions_PileUpCalc);
   fChain->SetBranchAddress("nTrueInteractions_PileUpCalc", &nTrueInteractions, &b_nTrueInteractions_PileUpCalc);
   fChain->SetBranchAddress("weight_PU_PileUpCalc", &weight_PU, &b_weight_PU_PileUpCalc);

   fChain->SetBranchAddress("dataEE_DileptonCalc", &dataEE, &b_dataEE_DileptonCalc);
   fChain->SetBranchAddress("dataEM_DileptonCalc", &dataEM, &b_dataEM_DileptonCalc);
   fChain->SetBranchAddress("dataMM_DileptonCalc", &dataMM, &b_dataMM_DileptonCalc);

   fChain->SetBranchAddress("trigEE_DileptonCalc", &trigEE, &b_trigEE_DileptonCalc);
   fChain->SetBranchAddress("trigEM_DileptonCalc", &trigEM, &b_trigEM_DileptonCalc);
   fChain->SetBranchAddress("trigMM_DileptonCalc", &trigMM, &b_trigMM_DileptonCalc);

   fChain->SetBranchAddress("corr_met_DileptonCalc", &corr_met, &b_corr_met_DileptonCalc);
   fChain->SetBranchAddress("corr_met_phi_DileptonCalc", &corr_met_phi, &b_corr_met_phi_DileptonCalc);

   fChain->SetBranchAddress("met_DileptonCalc", &PF_met_pt, &b_met_DileptonCalc);
   fChain->SetBranchAddress("met_phi_DileptonCalc", &PF_met_phi, &b_met_phi_DileptonCalc);

   //Electrons
   fChain->SetBranchAddress("elChargeConsistent_DileptonCalc", &elChargeConsistent, &b_elChargeConsistent_DileptonCalc);
   fChain->SetBranchAddress("elCharge_DileptonCalc", &elCharge, &b_elCharge_DileptonCalc);
   fChain->SetBranchAddress("elIsEBEE_DileptonCalc", &elIsEBEE, &b_elIsEBEE_DileptonCalc);
   fChain->SetBranchAddress("elMHits_DileptonCalc", &elMHits, &b_elMHits_DileptonCalc);
   fChain->SetBranchAddress("elNotConversion_DileptonCalc", &elNotConversion, &b_elNotConversion_DileptonCalc);
   fChain->SetBranchAddress("elQuality_DileptonCalc", &elQuality, &b_elQuality_DileptonCalc);
   fChain->SetBranchAddress("elVtxFitConv_DileptonCalc", &elVtxFitConv, &b_elVtxFitConv_DileptonCalc);
   fChain->SetBranchAddress("elD0_DileptonCalc", &elD0, &b_elD0_DileptonCalc);
   fChain->SetBranchAddress("elDZ_DileptonCalc", &elDZ, &b_elDZ_DileptonCalc);
   fChain->SetBranchAddress("elDeta_DileptonCalc", &elDeta, &b_elDeta_DileptonCalc);
   fChain->SetBranchAddress("elDphi_DileptonCalc", &elDphi, &b_elDphi_DileptonCalc);
   fChain->SetBranchAddress("elDxy_DileptonCalc", &elDxy, &b_elDxy_DileptonCalc);
   fChain->SetBranchAddress("elEnergy_DileptonCalc", &elEnergy, &b_elEnergy_DileptonCalc);
   fChain->SetBranchAddress("elEta_DileptonCalc", &elEta, &b_elEta_DileptonCalc);
   fChain->SetBranchAddress("elHoE_DileptonCalc", &elHoE, &b_elHoE_DileptonCalc);
   fChain->SetBranchAddress("elOoemoop_DileptonCalc", &elOoemoop, &b_elOoemoop_DileptonCalc);
   fChain->SetBranchAddress("elPhi_DileptonCalc", &elPhi, &b_elPhi_DileptonCalc);
   fChain->SetBranchAddress("elPt_DileptonCalc", &elPt, &b_elPt_DileptonCalc);
   fChain->SetBranchAddress("elSihih_DileptonCalc", &elSihih, &b_elSihih_DileptonCalc);

   fChain->SetBranchAddress("elAEff_DileptonCalc", &elAEff, &b_elAEff_DileptonCalc);
   fChain->SetBranchAddress("elChIso_DileptonCalc", &elChIso, &b_elChIso_DileptonCalc);
   fChain->SetBranchAddress("elNhIso_DileptonCalc", &elNhIso, &b_elNhIso_DileptonCalc);
   fChain->SetBranchAddress("elPhIso_DileptonCalc", &elPhIso, &b_elPhIso_DileptonCalc);
   fChain->SetBranchAddress("elRelIso_DileptonCalc", &elRelIso, &b_elRelIso_DileptonCalc);
   fChain->SetBranchAddress("elRhoIso_DileptonCalc", &elRhoIso, &b_elRhoIso_DileptonCalc);

   //Muons
   fChain->SetBranchAddress("muCharge_DileptonCalc", &muCharge, &b_muCharge_DileptonCalc);
   fChain->SetBranchAddress("muGlobal_DileptonCalc", &muGlobal, &b_muGlobal_DileptonCalc);
   fChain->SetBranchAddress("muNMatchedStations_DileptonCalc", &muNMatchedStations, &b_muNMatchedStations_DileptonCalc);
   fChain->SetBranchAddress("muNTrackerLayers_DileptonCalc", &muNTrackerLayers, &b_muNTrackerLayers_DileptonCalc);
   fChain->SetBranchAddress("muNValMuHits_DileptonCalc", &muNValMuHits, &b_muNValMuHits_DileptonCalc);
   fChain->SetBranchAddress("muNValPixelHits_DileptonCalc", &muNValPixelHits, &b_muNValPixelHits_DileptonCalc);
   fChain->SetBranchAddress("muChi2_DileptonCalc", &muChi2, &b_muChi2_DileptonCalc);
   fChain->SetBranchAddress("muDxy_DileptonCalc", &muDxy, &b_muDxy_DileptonCalc);
   fChain->SetBranchAddress("muDz_DileptonCalc", &muDz, &b_muDz_DileptonCalc);
   fChain->SetBranchAddress("muEnergy_DileptonCalc", &muEnergy, &b_muEnergy_DileptonCalc);
   fChain->SetBranchAddress("muEta_DileptonCalc", &muEta, &b_muEta_DileptonCalc);
   fChain->SetBranchAddress("muPhi_DileptonCalc", &muPhi, &b_muPhi_DileptonCalc);
   fChain->SetBranchAddress("muPt_DileptonCalc", &muPt, &b_muPt_DileptonCalc);
   fChain->SetBranchAddress("muQuality_DileptonCalc", &muQuality, &b_muQuality_DileptonCalc);

   fChain->SetBranchAddress("muChIso_DileptonCalc", &muChIso, &b_muChIso_DileptonCalc);
   fChain->SetBranchAddress("muGIso_DileptonCalc", &muGIso, &b_muGIso_DileptonCalc);
   fChain->SetBranchAddress("muNhIso_DileptonCalc", &muNhIso, &b_muNhIso_DileptonCalc);
   fChain->SetBranchAddress("muPuIso_DileptonCalc", &muPuIso, &b_muPuIso_DileptonCalc);
   fChain->SetBranchAddress("muRelIso_DileptonCalc", &muRelIso, &b_muRelIso_DileptonCalc);

   fChain->SetBranchAddress("AK5JetTBag_DileptonCalc", &AK5JetTBag, &b_AK5JetTBag_DileptonCalc);
   fChain->SetBranchAddress("AK5JetEnergy_DileptonCalc", &AK5JetEnergy, &b_AK5JetEnergy_DileptonCalc);
   fChain->SetBranchAddress("AK5JetEta_DileptonCalc", &AK5JetEta, &b_AK5JetEta_DileptonCalc);
   fChain->SetBranchAddress("AK5JetPhi_DileptonCalc", &AK5JetPhi, &b_AK5JetPhi_DileptonCalc);
   fChain->SetBranchAddress("AK5JetPt_DileptonCalc", &AK5JetPt, &b_AK5JetPt_DileptonCalc);

   fChain->SetBranchAddress("AK5JetnChHad_DileptonCalc", &AK5JetnChHad, &b_AK5JetnChHad_DileptonCalc);
   fChain->SetBranchAddress("AK5JetnNeuHad_DileptonCalc", &AK5JetnNeuHad, &b_AK5JetnNeuHad_DileptonCalc);
   fChain->SetBranchAddress("AK5JetRCN_DileptonCalc", &AK5JetRCN, &b_AK5JetRCN_DileptonCalc);

   fChain->SetBranchAddress("nPV_DileptonCalc", &nPV, &b_nPV_DileptonCalc);

   //Gen Info
   fChain->SetBranchAddress("genID_TopEventReweightCalc", &genID, &b_genID_DileptonCalc);
   fChain->SetBranchAddress("genIndex_TopEventReweightCalc", &genIndex, &b_genIndex_DileptonCalc);
   fChain->SetBranchAddress("genMotherID_TopEventReweightCalc", &genMotherID, &b_genMotherID_DileptonCalc);
   fChain->SetBranchAddress("genMotherIndex_TopEventReweightCalc", &genMotherIndex, &b_genMotherIndex_DileptonCalc);
   fChain->SetBranchAddress("genStatus_TopEventReweightCalc", &genStatus, &b_genStatus_DileptonCalc);

   fChain->SetBranchAddress("genEnergy_TopEventReweightCalc", &genEnergy, &b_genEnergy_DileptonCalc);
   fChain->SetBranchAddress("genEta_TopEventReweightCalc", &genEta, &b_genEta_DileptonCalc);
   fChain->SetBranchAddress("genPhi_TopEventReweightCalc", &genPhi, &b_genPhi_DileptonCalc);
   fChain->SetBranchAddress("genPt_TopEventReweightCalc", &genPt, &b_genPt_DileptonCalc);
   int j = fChain->SetBranchAddress("eventWeight_TopEventReweightCalc", &higgsWeight, &b_higgsWeight_DileptonCalc);
   if (j==3) {
    cout <<"This is a sample with Higgs reweight variable\n";
    isHiggs=true;
   } else {
     cout <<"This is a default sample without Higgs reweight variable\n";
     isHiggs=false;
   }

   j = fChain->SetBranchAddress("Q_PdfCalc", &Q, &b_Q_PdfCalc);
   if (j==3) {
     cout <<"This is a sample with PDF reweight variable\n";
     withPDF=true;
     fChain->SetBranchAddress("id1_PdfCalc", &id1, &b_id1_PdfCalc);
     fChain->SetBranchAddress("id2_PdfCalc", &id2, &b_id2_PdfCalc);
     fChain->SetBranchAddress("x1_PdfCalc", &x1, &b_x1_PdfCalc);
     fChain->SetBranchAddress("x2_PdfCalc", &x2, &b_x2_PdfCalc);
     fChain->SetBranchAddress("pdf1_PdfCalc", &pdf1, &b_pdf1_PdfCalc);
     fChain->SetBranchAddress("pdf2_PdfCalc", &pdf2, &b_pdf2_PdfCalc);
//      fChain->SetBranchAddress("PdfWeightPlus_MRST2006nnlo_PdfCalc", &PdfWeightPlus_MRST2006nnlo, &b_PdfWeightPlus_MRST2006nnlo_PdfCalc);
//      fChain->SetBranchAddress("PdfWeightPlus_NNPDF10_PdfCalc", &PdfWeightPlus_NNPDF10, &b_PdfWeightPlus_NNPDF10_PdfCalc);
//      fChain->SetBranchAddress("PdfWeightPlus_cteq66_PdfCalc", &PdfWeightPlus_cteq66, &b_PdfWeightPlus_cteq66_PdfCalc);
//      fChain->SetBranchAddress("PdfWeightsVec_MRST2006nnlo_PdfCalc", &PdfWeightsVec_MRST2006nnlo, &b_PdfWeightsVec_MRST2006nnlo_PdfCalc);
//      fChain->SetBranchAddress("PdfWeightsVec_NNPDF10_PdfCalc", &PdfWeightsVec_NNPDF10, &b_PdfWeightsVec_NNPDF10_PdfCalc);
//      fChain->SetBranchAddress("PdfWeightsVec_cteq66_PdfCalc", &PdfWeightsVec_cteq66, &b_PdfWeightsVec_cteq66_PdfCalc);
   } else {
     cout <<"This is a default sample without PDF reweight variable\n";
     withPDF=false;
   }


  Notify();
}

Bool_t treetop::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void treetop::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t treetop::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
