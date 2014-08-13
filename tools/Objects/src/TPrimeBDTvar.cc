#include "../interface/TPrimeBDTvar.h"


BDTVar::BDTVar(const double scale) : scale_(scale){
  MVA = new TTree("MVA","BDTInput");
  MVA->Branch("run", &run  , "run/I");
  MVA->Branch("lumi", &lumi , "lumi/I");
  MVA->Branch("event", &event, "event/I");
  MVA->Branch("weight", &weight, "weight/F");
  MVA->Branch("met",&PF_met_pt,"PF_met_pt/F");
  MVA->Branch("nBTag", &nBTag, "nBTag/I");
  MVA->Branch("nJets", &nJets, "nJets/I");
  MVA->Branch("nLeptons", &nLeptons, "nLeptons/I");
  MVA->Branch("lepton1_eta", &lepton1_eta, "lepton1_eta/F");
  MVA->Branch("lepton2_eta", &lepton2_eta, "lepton2_eta/F");
  MVA->Branch("lepton3_eta", &lepton3_eta, "lepton3_eta/F");
  MVA->Branch("lepton4_eta", &lepton4_eta, "lepton4_eta/F");
  MVA->Branch("lepton1_pt",  &lepton1_pt,  "lepton1_pt/F");
  MVA->Branch("lepton2_pt",  &lepton2_pt,  "lepton2_pt/F");
  MVA->Branch("lepton3_pt",  &lepton3_pt,  "lepton3_pt/F");
  MVA->Branch("lepton4_pt",  &lepton4_pt,  "lepton4_pt/F");
  MVA->Branch("lepton1_phi", &lepton1_phi, "lepton1_phi/F");
  MVA->Branch("lepton2_phi", &lepton2_phi, "lepton2_phi/F");
  MVA->Branch("lepton3_phi", &lepton3_phi, "lepton3_phi/F");
  MVA->Branch("lepton4_phi", &lepton4_phi, "lepton4_phi/F");
  MVA->Branch("lepton1_energy", &lepton1_energy, "lepton1_energy/F");
  MVA->Branch("lepton2_energy", &lepton2_energy, "lepton2_energy/F");
  MVA->Branch("lepton3_energy", &lepton3_energy, "lepton3_energy/F");
  MVA->Branch("lepton4_energy", &lepton4_energy, "lepton4_energy/F");
  MVA->Branch("jet1_pt", &jet1_pt, "jet1_pt/F");
  MVA->Branch("jet2_pt", &jet2_pt, "jet2_pt/F");
  MVA->Branch("jet3_pt", &jet3_pt, "jet3_pt/F");
  MVA->Branch("jet4_pt", &jet4_pt, "jet4_pt/F");
  MVA->Branch("jet1_eta", &jet1_eta, "jet1_eta/F");
  MVA->Branch("jet2_eta", &jet2_eta, "jet2_eta/F");
  MVA->Branch("jet3_eta", &jet3_eta, "jet3_eta/F");
  MVA->Branch("jet4_eta", &jet4_eta, "jet4_eta/F");
  MVA->Branch("jet1_phi", &jet1_phi, "jet1_phi/F");
  MVA->Branch("jet2_phi", &jet2_phi, "jet2_phi/F");
  MVA->Branch("jet3_phi", &jet3_phi, "jet3_phi/F");
  MVA->Branch("jet4_phi", &jet4_phi, "jet4_phi/F");
  MVA->Branch("jet1_energy", &jet1_energy, "jet1_energy/F");
  MVA->Branch("jet2_energy", &jet2_energy, "jet2_energy/F");
  MVA->Branch("jet3_energy", &jet3_energy, "jet3_energy/F");
  MVA->Branch("jet4_energy", &jet4_energy, "jet4_energy/F");

  MVA->Branch("jet5_pt", &jet5_pt, "jet5_pt/F");
  MVA->Branch("jet5_eta", &jet5_eta, "jet5_eta/F");
  MVA->Branch("jet5_phi", &jet5_phi, "jet5_phi/F");
  MVA->Branch("jet5_energy", &jet5_energy, "jet5_energy/F");
  MVA->Branch("jet6_pt", &jet6_pt, "jet6_pt/F");
  MVA->Branch("jet6_eta", &jet6_eta, "jet6_eta/F");
  MVA->Branch("jet6_phi", &jet6_phi, "jet6_phi/F");
  MVA->Branch("jet6_energy", &jet6_energy, "jet6_energy/F");
  MVA->Branch("jet7_pt", &jet7_pt, "jet7_pt/F");
  MVA->Branch("jet7_eta", &jet7_eta, "jet7_eta/F");
  MVA->Branch("jet7_phi", &jet7_phi, "jet7_phi/F");
  MVA->Branch("jet7_energy", &jet7_energy, "jet7_energy/F");
  MVA->Branch("jet8_pt", &jet8_pt, "jet8_pt/F");
  MVA->Branch("jet8_eta", &jet8_eta, "jet8_eta/F");
  MVA->Branch("jet8_phi", &jet8_phi, "jet8_phi/F");
  MVA->Branch("jet8_energy", &jet8_energy, "jet8_energy/F");
  MVA->Branch("jets_btag", &jets_btag, "jets_btag/I");

  MVA->Branch("jet1_nChHad", jet_nChHad, "jet1_nChHad/I");
  MVA->Branch("jet1_nNeHad", jet_nNeHad, "jet1_nNeHad/I");
  MVA->Branch("jet1_RCN",    jet_RCN ,    "jet1_RCN/F");
  MVA->Branch("jet2_nChHad", jet_nChHad+1 , "jet2_nChHad/I");
  MVA->Branch("jet2_nNeHad", jet_nNeHad+1 , "jet2_nNeHad/I");
  MVA->Branch("jet2_RCN",    jet_RCN   +1 ,    "jet2_RCN/F");
  MVA->Branch("jet3_nChHad", jet_nChHad+2 , "jet3_nChHad/I");
  MVA->Branch("jet3_nNeHad", jet_nNeHad+2 , "jet3_nNeHad/I");
  MVA->Branch("jet3_RCN",    jet_RCN   +2 ,    "jet3_RCN/F");
  MVA->Branch("jet4_nChHad", jet_nChHad+3 , "jet4_nChHad/I");
  MVA->Branch("jet4_nNeHad", jet_nNeHad+3 , "jet4_nNeHad/I");
  MVA->Branch("jet4_RCN",    jet_RCN   +3 ,    "jet4_RCN/F");
  MVA->Branch("jet5_nChHad", jet_nChHad+4 , "jet5_nChHad/I");
  MVA->Branch("jet5_nNeHad", jet_nNeHad+4 , "jet5_nNeHad/I");
  MVA->Branch("jet5_RCN",    jet_RCN   +4 ,    "jet5_RCN/F");
  MVA->Branch("jet6_nChHad", jet_nChHad+5 , "jet6_nChHad/I");
  MVA->Branch("jet6_nNeHad", jet_nNeHad+5 , "jet6_nNeHad/I");
  MVA->Branch("jet6_RCN",    jet_RCN   +5 ,    "jet6_RCN/F");
  MVA->Branch("jet7_nChHad", jet_nChHad+6 , "jet7_nChHad/I");
  MVA->Branch("jet7_nNeHad", jet_nNeHad+6 , "jet7_nNeHad/I");
  MVA->Branch("jet7_RCN",    jet_RCN   +6 ,    "jet7_RCN/F");
  MVA->Branch("jet8_nChHad", jet_nChHad+7 , "jet8_nChHad/I");
  MVA->Branch("jet8_nNeHad", jet_nNeHad+7 , "jet8_nNeHad/I");
  MVA->Branch("jet8_RCN",    jet_RCN   +7 ,    "jet8_RCN/F");

  MVA->Branch("minDRJetLepton", &minDRJetLepton, "minDRJetLepton/F");
  MVA->Branch("minDRBJetLepton", &minDRBJetLepton, "minDRBJetLepton/F");
  MVA->Branch("minMLB", &minMLB, "minMLB/F");
  MVA->Branch("maxMLB", &maxMLB, "maxMLB/F");
  MVA->Branch("HT", &HT, "HT/F");
  MVA->Branch("leptonSumMass", &leptonSumMass, "leptonSumMass/F");
  MVA->Branch("CATopJetsBtag_disc1", &CATopJetsBtag_disc1, "CATopJetsBtag_disc1/F"); 
  MVA->Branch("CATopJetsBtag_disc2", &CATopJetsBtag_disc2, "CATopJetsBtag_disc2/F"); 
  MVA->Branch("CATopJets1_pt", &CATopJets1_pt, "CATopJets1_pt/F"); 
  MVA->Branch("CATopJets2_pt", &CATopJets2_pt, "CATopJets2_pt/F");
  MVA->Branch("CATopJets1_eta", &CATopJets1_eta, "CATopJets1_eta/F"); 
  MVA->Branch("CATopJets2_eta", &CATopJets2_eta, "CATopJets2_eta/F");
  MVA->Branch("CATopJets1_phi", &CATopJets1_phi, "CATopJets1_phi/F"); 
  MVA->Branch("CATopJets2_phi", &CATopJets2_phi, "CATopJets2_phi/F");
  MVA->Branch("CATopJets1_energy", &CATopJets1_energy, "CATopJets1_energy/F"); 
  MVA->Branch("CATopJets2_energy", &CATopJets2_energy, "CATopJets2_energy/F"); 
  MVA->Branch("CATopJets1_jetMass", &CATopJets1_jetMass, "CATopJets1_jetMass/F"); 
  MVA->Branch("CATopJets2_jetMass", &CATopJets2_jetMass, "CATopJets2_jetMass/F"); 
  MVA->Branch("CATopJets1_massDrop", &CATopJets1_massDrop, "CATopJets1_massDrop/F"); 
  MVA->Branch("CATopJets2_massDrop", &CATopJets2_massDrop, "CATopJets2_massDrop/F"); 
  MVA->Branch("CATopJets1_nbrConst", &CATopJets1_nbrConst, "CATopJets1_nbrConst/F"); 
  MVA->Branch("CATopJets2_nbrConst", &CATopJets2_nbrConst, "CATopJets2_nbrConst/F"); 
  MVA->Branch("CATopJets1_minMass", &CATopJets1_minMass, "CATopJets1_minMass/F"); 
  MVA->Branch("CATopJets2_minMass", &CATopJets2_minMass, "CATopJets2_minMass/F"); 

  MVA->Branch("CAWBtag_disc1", &CAWBtag_disc1, "CAWBtag_disc1/F"); 
  MVA->Branch("CAWBtag_disc2", &CAWBtag_disc2, "CAWBtag_disc2/F");
  MVA->Branch("CAWJets1_pt", &CAWJets1_pt, "CAWJets1_pt/F"); 
  MVA->Branch("CAWJets2_pt", &CAWJets2_pt, "CAWJets2_pt/F");
  MVA->Branch("CAWJets1_eta", &CAWJets1_eta, "CAWJets1_eta/F"); 
  MVA->Branch("CAWJets2_eta", &CAWJets2_eta, "CAWJets2_eta/F");
  MVA->Branch("CAWJets1_phi", &CAWJets1_phi, "CAWJets1_phi/F"); 
  MVA->Branch("CAWJets2_phi", &CAWJets2_phi, "CAWJets2_phi/F");
  MVA->Branch("CAWJets1_energy", &CAWJets1_energy, "CAWJets1_energy/F"); 
  MVA->Branch("CAWJets2_energy", &CAWJets2_energy, "CAWJets2_energy/F"); 
  MVA->Branch("CAWJets1_jetMass", &CAWJets1_jetMass, "CAWJets1_jetMass/F"); 
  MVA->Branch("CAWJets2_jetMass", &CAWJets2_jetMass, "CAWJets2_jetMass/F"); 
  MVA->Branch("CAWJets1_massDrop", &CAWJets1_massDrop, "CAWJets1_massDrop/F"); 
  MVA->Branch("CAWJets2_massDrop", &CAWJets2_massDrop, "CAWJets2_massDrop/F"); 
  MVA->Branch("CAWJets1_nbrConst", &CAWJets1_nbrConst, "CAWJets1_nbrConst/F"); 
  MVA->Branch("CAWJets2_nbrConst", &CAWJets2_nbrConst, "CAWJets2_nbrConst/F"); 
//   MVA->Branch("CAWJets1_minMass", &CAWJets1_minMass, "CAWJets1_minMass/F"); 
//   MVA->Branch("CAWJets2_minMass", &CAWJets2_minMass, "CAWJets2_minMass/F"); 

  MVA->Branch("nCAWJets", &nCAWJets, "nCAWJets/I");
  MVA->Branch("nCATopJets", &nCATopJets, "nCATopJets/I");
  MVA->Branch("amwt", &amwt, "amwt/F");
  MVA->Branch("amwtPeakWeight", &amwtPeakWeight, "amwtPeakWeight/F");
  MVA->Branch("leptonSum",&leptonSum	,"leptonSum/F");
  MVA->Branch("jetSum",&jetSum	,"jetSum/F");
  MVA->Branch("leptonMETSum",&leptonMETSum,"leptonMETSum/F");
  MVA->Branch("leptonJetsSum",&leptonJetsSum,"leptonJetsSum/F");
  MVA->Branch("jetsMETSum",&jetsMETSum	,"jetsMETSum/F");
  MVA->Branch("leptonJetsMETSum",&leptonJetsMETSum,"leptonJetsMETSum/F");
  MVA->Branch("isMuon1", &isMuon1, "isMuon1/I");
  MVA->Branch("isMuon2", &isMuon2, "isMuon2/I");
  MVA->Branch("isMuon3", &isMuon3, "isMuon3/I");
  MVA->Branch("isMuon4", &isMuon4, "isMuon4/I");
  MVA->Branch("isolation1", &isolation1, "isolation1/F");
  MVA->Branch("isolation2", &isolation2, "isolation2/F");
  MVA->Branch("isolation3", &isolation3, "isolation3/F");
  MVA->Branch("isolation4", &isolation4, "isolation4/F");
  MVA->Branch("eventType", &eventType, "eventType/I");
  MVA->Branch("Q", &Q);
  MVA->Branch("id1", &id1);
  MVA->Branch("id2", &id2);
  MVA->Branch("x1", &x1);
  MVA->Branch("x2", &x2);
  MVA->Branch("pdf1", &pdf1);
  MVA->Branch("pdf2", &pdf2);
  
}

void BDTVar::FillBranches(TprimeEvent* teve, treetop* theTree,
	TString type, bool doFakes){
  //return;
// cout << "start\n";
  run  = theTree->run  ;
  lumi = theTree->lumi ;
  event= theTree->event;


  eventType = -1;
  if (doFakes) {
    if (type == "Nt00") eventType=0;
    if (type == "Nt01") eventType=1;
    if (type == "Nt10") eventType=2;
    if (type == "Nt11") eventType=3;
    if (type == "Nt0") eventType=0;
    if (type == "Nt1") eventType=1;
    if (type == "Nt2") eventType=2;
    if (type == "Nt3") eventType=3;
  } else eventType=3;

  weight = teve->weight * scale_;
  nLeptons = teve->theGoodLeptons.size();
  isMuon1 = 0;
  isMuon2 = 0;
  isMuon3 = 0;
  isMuon4 = 0;
  if(nLeptons > 0){
    if (teve->theGoodLeptons[0]->isMuon()==1){
      if(teve->theGoodLeptons[0]->charge>0){
      isMuon1 = 13;  
    }
      else if(teve->theGoodLeptons[0]->charge<0){
      isMuon1 = -13; 
    }
  }
    else if (teve->theGoodLeptons[0]->isMuon()==0){
      if(teve->theGoodLeptons[0]->charge>0){
      isMuon1 = 11; 
    }
      else if(teve->theGoodLeptons[0]->charge<0){
      isMuon1 = -11;
    }
   }
  }
  if(nLeptons > 1){
    if (teve->theGoodLeptons[1]->isMuon()==1){
      if(teve->theGoodLeptons[1]->charge>0){
      isMuon2 = 13; 
    }
      else if(teve->theGoodLeptons[1]->charge<0){
      isMuon2 = -13;
    }
  }
    else if (teve->theGoodLeptons[1]->isMuon()==0){
      if(teve->theGoodLeptons[1]->charge>0){
      isMuon2 = 11;
      }
      else if(teve->theGoodLeptons[1]->charge<0){
      isMuon2 = -11;
    }
   } 
  } 
 if(nLeptons > 2 ){
  if (teve->theGoodLeptons[2]->isMuon()==1){
    if(teve->theGoodLeptons[2]->charge>0){
    isMuon3 = 13;
    }
    else if(teve->theGoodLeptons[2]->charge<0){
    isMuon3 = -13;
    }
  }
  else if (teve->theGoodLeptons[2]->isMuon()==0){
   if(teve->theGoodLeptons[2]->charge>0){
    isMuon3 = 11;
    }
   else if(teve->theGoodLeptons[2]->charge<0){
    isMuon3 = -11;
    }
   }
  }
  if(nLeptons > 3 ){
   if (teve->theGoodLeptons[3]->isMuon()==1){
     if(teve->theGoodLeptons[3]->charge>0){
     isMuon4 = 13;
     }
   else if(teve->theGoodLeptons[3]->charge<0){
     isMuon4 = -13;
    }
   }
  else if (teve->theGoodLeptons[3]->isMuon()==0){
   if(teve->theGoodLeptons[3]->charge>0){
    isMuon4 = 11;  
    }
   else if(teve->theGoodLeptons[3]->charge<0){
    isMuon4 = -11;
    }
   }
  }
  PF_met_pt = teve->met;
  nBTag = teve->nBTags();
  nJets = teve->vGoodJets.size();
  //nLeptons = teve->theGoodLeptons.size();
  minDRJetLepton   = teve->mindrJetLepton(teve->vGoodJets);
  minDRBJetLepton  = teve->mindrBJetLepton();
  minMLB           = teve->minMlb();
  maxMLB           = teve->maxMlb();

  HT               = teve->ht();
  leptonSumMass    = teve->leptonSum().M();
  nCAWJets         = teve->vGoodCAWJets.size();
  nCATopJets       = teve->vGoodCATopJets.size();
  CATopJetsBtag_disc1 = 0.0;
  CATopJetsBtag_disc2 = 0.0;
  CATopJets1_pt = 0.0;
  CATopJets2_pt = 0.0;
  CATopJets1_phi = 0.0;
  CATopJets2_phi = 0.0;
  CATopJets1_eta = 0.0;
  CATopJets2_eta = 0.0;
  CATopJets1_energy = 0.0;
  CATopJets2_energy = 0.0;
  CATopJets1_jetMass=0.;  CATopJets2_jetMass=0.;  CATopJets1_massDrop=0.;  CATopJets2_massDrop=0.;  
  CATopJets1_minMass=0.;  CATopJets2_minMass=0.;
  CATopJets1_nbrConst=0;  CATopJets2_nbrConst=0;
  if (nCATopJets > 0) CATopJetsBtag_disc1 = teve->vGoodCATopJets[0]->csvDiscr;
  if (nCATopJets > 1) CATopJetsBtag_disc2 = teve->vGoodCATopJets[1]->csvDiscr; 
  if (nCATopJets > 0) CATopJets1_pt = teve->vGoodCATopJets[0]->lv.Pt();
  if (nCATopJets > 1) CATopJets2_pt = teve->vGoodCATopJets[1]->lv.Pt();
  if (nCATopJets > 0) CATopJets1_phi = teve->vGoodCATopJets[0]->lv.Phi();
  if (nCATopJets > 1) CATopJets2_phi = teve->vGoodCATopJets[1]->lv.Phi();
  if (nCATopJets > 0) CATopJets1_eta = teve->vGoodCATopJets[0]->lv.Eta();
  if (nCATopJets > 1) CATopJets2_eta = teve->vGoodCATopJets[1]->lv.Eta();
  if (nCATopJets > 0) {
    CATopJets1_energy = teve->vGoodCATopJets[0]->lv.Energy();
    CATopJets1_jetMass  = teve->vGoodCATopJets[0]->mass;
//     CATopJets1_massDrop = teve->vGoodCATopJets[0]->Energy();
    CATopJets1_nbrConst = teve->vGoodCATopJets[0]->nDaughters;
    CATopJets1_minMass = teve->vGoodCATopJets[0]->minPairMass;
  }
  if (nCATopJets > 1) {
    CATopJets2_energy = teve->vGoodCATopJets[1]->lv.Energy();
    CATopJets2_jetMass  = teve->vGoodCATopJets[0]->mass;
//     CATopJets2_massDrop = teve->vGoodCATopJets[0]->Energy();
    CATopJets2_nbrConst = teve->vGoodCATopJets[0]->nDaughters;
    CATopJets2_minMass = teve->vGoodCATopJets[0]->minPairMass;
  }
  CAWBtag_disc1 = 0.0;
  CAWBtag_disc2 = 0.0;
  CAWJets1_pt = 0.0;
  CAWJets2_pt = 0.0;
  CAWJets1_phi = 0.0;
  CAWJets2_phi = 0.0;
  CAWJets1_eta = 0.0;
  CAWJets2_eta = 0.0;
  CAWJets1_energy = 0.0;
  CAWJets2_energy = 0.0;
  amwt = 0.0; 
  amwtPeakWeight = 0.0;
    lepton1_eta	=0.;
  lepton2_eta	=0.;
  lepton3_eta	=0.;
  lepton4_eta	=0.;
  lepton1_pt	=0.;
  lepton2_pt	=0.;
  lepton3_pt	=0.;
  lepton4_pt	=0.;
  lepton1_phi	=0.;
  lepton2_phi	=0.;
  lepton3_phi	=0.;
  lepton4_phi	=0.;
  lepton1_energy  =0.;
  lepton2_energy  =0.;
  lepton3_energy  =0.;
  lepton4_energy  =0.;
  isolation1 =  isolation2 =  isolation3 =  isolation4 = 0;
  jets_btag=0;
  jet1_pt 	=0.;
  jet2_pt 	=0.;
  jet3_pt 	=0.;
  jet4_pt 	=0.;
  jet1_eta	=0.;
  jet2_eta	=0.;
  jet3_eta	=0.;
  jet4_eta	=0.;
  jet1_phi	=0.;
  jet2_phi	=0.;
  jet3_phi	=0.;
  jet4_phi	=0.;
  jet1_energy	=0.;
  jet2_energy	=0.;
  jet3_energy	=0.;
  jet4_energy	=0.;
  jet5_pt 	=0.;
  jet6_pt 	=0.;
  jet7_pt 	=0.;
  jet8_pt 	=0.;
  jet5_eta	=0.;
  jet6_eta	=0.;
  jet7_eta	=0.;
  jet8_eta	=0.;
  jet5_phi	=0.;
  jet6_phi	=0.;
  jet7_phi	=0.;
  jet8_phi	=0.;
  jet5_energy	=0.;
  jet6_energy	=0.;
  jet7_energy	=0.;
  jet8_energy	=0.;
  CAWJets1_jetMass=0.;  CAWJets2_jetMass=0.;  CAWJets1_massDrop=0.;  CAWJets2_massDrop=0.;  
//   CAWJets1_minMass=0.;  CAWJets2_minMass=0.;
  CAWJets1_nbrConst=0;  CAWJets2_nbrConst=0;

  if (nCAWJets > 0) CAWBtag_disc1 = teve->vGoodCAWJets[0]->csvDiscr;
  if (nCAWJets > 1) CAWBtag_disc2 = teve->vGoodCAWJets[1]->csvDiscr;
  if (nCAWJets > 0) CAWJets1_pt   = teve->vGoodCAWJets[0]->lv.Pt();
  if (nCAWJets > 1) CAWJets2_pt   = teve->vGoodCAWJets[1]->lv.Pt();
  if (nCAWJets > 0) CAWJets1_phi  = teve->vGoodCAWJets[0]->lv.Phi();
  if (nCAWJets > 1) CAWJets2_phi  = teve->vGoodCAWJets[1]->lv.Phi();
  if (nCAWJets > 0) CAWJets1_eta  = teve->vGoodCAWJets[0]->lv.Eta();
  if (nCAWJets > 1) CAWJets2_eta  = teve->vGoodCAWJets[1]->lv.Eta();
  if (nCAWJets > 0) {
    CAWJets1_energy = teve->vGoodCAWJets[0]->lv.Energy();
    CAWJets1_jetMass  = teve->vGoodCAWJets[0]->mass;
//     CAWJets1_massDrop = teve->vGoodCAWJets[0]->Energy();
    CAWJets1_nbrConst = teve->vGoodCAWJets[0]->nDaughters;
//     CAWJets1_minMass = teve->vGoodCAWJets[0]->minPairMass;
  }

  if (nCAWJets > 1) {
    CAWJets2_energy = teve->vGoodCAWJets[1]->lv.Energy(); 
    CAWJets2_jetMass  = teve->vGoodCAWJets[0]->mass;
//     CAWJets2_massDrop = teve->vGoodCAWJets[0]->Energy();
    CAWJets2_nbrConst = teve->vGoodCAWJets[0]->nDaughters;
//     CAWJets2_minMass = teve->vGoodCAWJets[0]->minPairMass;
  }

  amwtPeakWeight = teve->amwtPeakWeight;
  if (teve->amwtPeakWeight>0.) amwt = teve->amwtPeakMass;   
  if (nLeptons > 0) lepton1_eta     = teve->theGoodLeptons[0]->lv.Eta();
  if (nLeptons > 1) lepton2_eta     = teve->theGoodLeptons[1]->lv.Eta();
  if (nLeptons > 2) lepton3_eta     = teve->theGoodLeptons[2]->lv.Eta();
  if (nLeptons > 3) lepton4_eta     = teve->theGoodLeptons[3]->lv.Eta();
  if (nLeptons > 0) lepton1_pt      = teve->theGoodLeptons[0]->lv.Pt();
  if (nLeptons > 1) lepton2_pt      = teve->theGoodLeptons[1]->lv.Pt();
  if (nLeptons > 2) lepton3_pt      = teve->theGoodLeptons[2]->lv.Pt();
  if (nLeptons > 3) lepton4_pt      = teve->theGoodLeptons[3]->lv.Pt();
  if (nLeptons > 0) lepton1_phi     = teve->theGoodLeptons[0]->lv.Phi();
  if (nLeptons > 1) lepton2_phi     = teve->theGoodLeptons[1]->lv.Phi();
  if (nLeptons > 2) lepton3_phi     = teve->theGoodLeptons[2]->lv.Phi();
  if (nLeptons > 3) lepton4_phi     = teve->theGoodLeptons[3]->lv.Phi();
  if (nLeptons > 0) lepton1_energy  = teve->theGoodLeptons[0]->lv.Energy();
  if (nLeptons > 1) lepton2_energy  = teve->theGoodLeptons[1]->lv.Energy();
  if (nLeptons > 2) lepton3_energy  = teve->theGoodLeptons[2]->lv.Energy();
  if (nLeptons > 3) lepton4_energy  = teve->theGoodLeptons[3]->lv.Energy();

  if (nLeptons > 0) isolation1  = teve->theGoodLeptons[0]->relIso;
  if (nLeptons > 1) isolation2  = teve->theGoodLeptons[1]->relIso;
  if (nLeptons > 2) isolation3  = teve->theGoodLeptons[2]->relIso;
  if (nLeptons > 3) isolation4  = teve->theGoodLeptons[3]->relIso;

  if (nJets    > 0) jet1_pt         = teve->vGoodJets[0]->lv.Pt();
  if (nJets    > 1) jet2_pt         = teve->vGoodJets[1]->lv.Pt();
  if (nJets    > 2) jet3_pt         = teve->vGoodJets[2]->lv.Pt();
  if (nJets    > 3) jet4_pt         = teve->vGoodJets[3]->lv.Pt();
  if (nJets    > 0) jet1_eta        = teve->vGoodJets[0]->lv.Eta();
  if (nJets    > 1) jet2_eta        = teve->vGoodJets[1]->lv.Eta();
  if (nJets    > 2) jet3_eta        = teve->vGoodJets[2]->lv.Eta();
  if (nJets    > 3) jet4_eta        = teve->vGoodJets[3]->lv.Eta();
  if (nJets    > 0) jet1_phi        = teve->vGoodJets[0]->lv.Phi();
  if (nJets    > 1) jet2_phi        = teve->vGoodJets[1]->lv.Phi();
  if (nJets    > 2) jet3_phi        = teve->vGoodJets[2]->lv.Phi();
  if (nJets    > 3) jet4_phi        = teve->vGoodJets[3]->lv.Phi();
  if (nJets    > 0) jet1_energy     = teve->vGoodJets[0]->lv.Energy();
  if (nJets    > 1) jet2_energy     = teve->vGoodJets[1]->lv.Energy();
  if (nJets    > 2) jet3_energy     = teve->vGoodJets[2]->lv.Energy();
  if (nJets    > 3) jet4_energy     = teve->vGoodJets[3]->lv.Energy();

  if (nJets    > 4) jet5_pt	    = teve->vGoodJets[4]->lv.Pt();
  if (nJets    > 5) jet6_pt	    = teve->vGoodJets[5]->lv.Pt();
  if (nJets    > 6) jet7_pt	    = teve->vGoodJets[6]->lv.Pt();
  if (nJets    > 7) jet8_pt	    = teve->vGoodJets[7]->lv.Pt();
  if (nJets    > 4) jet5_eta	    = teve->vGoodJets[4]->lv.Eta();
  if (nJets    > 5) jet6_eta	    = teve->vGoodJets[5]->lv.Eta();
  if (nJets    > 6) jet7_eta	    = teve->vGoodJets[6]->lv.Eta();
  if (nJets    > 7) jet8_eta	    = teve->vGoodJets[7]->lv.Eta();
  if (nJets    > 4) jet5_phi	    = teve->vGoodJets[4]->lv.Phi();
  if (nJets    > 5) jet6_phi	    = teve->vGoodJets[5]->lv.Phi();
  if (nJets    > 6) jet7_phi	    = teve->vGoodJets[6]->lv.Phi();
  if (nJets    > 7) jet8_phi	    = teve->vGoodJets[7]->lv.Phi();
  if (nJets    > 4) jet5_energy     = teve->vGoodJets[4]->lv.Energy();
  if (nJets    > 5) jet6_energy     = teve->vGoodJets[5]->lv.Energy();
  if (nJets    > 6) jet7_energy     = teve->vGoodJets[6]->lv.Energy();
  if (nJets    > 7) jet8_energy     = teve->vGoodJets[7]->lv.Energy();

  for (unsigned int j = 0;j<min((size_t)8,teve->vGoodJets.size());++j ) {
    jet_RCN[j] = teve->vGoodJets.at(j)->jetRCN;
    jet_nChHad[j] = teve->vGoodJets.at(j)->jetnChHad;
    jet_nNeHad[j] = teve->vGoodJets.at(j)->jetnNeuHad;
//     cout << jet_RCN[j]<< " "<<
// jet_nChHad[j]<< " "<<
// jet_nNeHad[j]<<endl;
  }
  
  for (int i=0;i<nJets;++i){
    if (teve->vGoodJets[i]->csvMedium != 0) jets_btag+= (1<<i);
  }
  
  leptonSum 		= teve->leptonScalSum();	    
  leptonMETSum 		= teve->leptonMETScalSum();
  leptonJetsSum 	= teve->leptonJetsScalSum();
  jetsMETSum 		= teve->jetsMETScalSum();
  leptonJetsMETSum	= teve->leptonJetsMETScalSum();

  if (theTree->withPDF){
    Q    = (theTree->Q	  );
    id1  = (theTree->id1  );
    id2  = (theTree->id2  );
    x1   = (theTree->x1	 );
    x2   = (theTree->x2	 );
    pdf1 = (theTree->pdf1);
    pdf2 = (theTree->pdf2);
  }

// cout << "fill\n";
  MVA->Fill();
// cout << "done\n";
 }

 void BDTVar::WriteBranches(){
  MVA->Write();
 } 
