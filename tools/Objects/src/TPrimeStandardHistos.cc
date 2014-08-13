#include "../interface/TPrimeStandardHistos.h"


void StandardHistos::initAllHistos(bool doFakes){
  vector <HistoData*> vhd(NHISTOS);

  vhd.at(hnEvents)      = new HistoData("nEvents",         1,  0.5,  1.5);

  vhd.at(hnPV)          = new HistoData("nPV",            50, -0.5,  49.5);
  vhd.at(hPUwgt)        = new HistoData("PUwght",         50,  0.,  3);

  vhd.at(Lep1Pt)    = new HistoData("Lep1Pt",   500,  0.,  500.);
  vhd.at(Lep2Pt)    = new HistoData("Lep2Pt",   500,  0.,  500.);
  vhd.at(Lep3Pt)    = new HistoData("Lep3Pt",   500,  0.,  500.);

  vhd.at(ElPt)      = new HistoData("ElPt",     500,  0.,  500.);
  vhd.at(MuPt)      = new HistoData("MuPt",     500,  0.,  500.);
  vhd.at(ElEta)     = new HistoData("ElEta",     60, -2.4, 2.4);
  vhd.at(MuEta)     = new HistoData("MuEta",     60, -2.4, 2.4);

  vhd.at(nElec)     = new HistoData("nElec",   9,  -0.5, 8.5);
  vhd.at(nMuon)     = new HistoData("nMuon",   9,  -0.5, 8.5);
  vhd.at(nLept)     = new HistoData("nLept",   9,  -0.5, 8.5);

  vhd.at(Jet1Pt)    = new HistoData("Jet1Pt",   200,  0.,  1000.);
  vhd.at(Jet2Pt)    = new HistoData("Jet2Pt",   200,  0.,  1000.);
  vhd.at(nJets)     = new HistoData("nJets",     15,  -0.5, 14.5);
  vhd.at(nBJets)    = new HistoData("nBJets",   10,  -0.5, 9.5);

  vhd.at(AK5Jet1Pt)     = new HistoData("AK5Jet1Pt",     480,  0.,  480.);
  vhd.at(AK5Jet2Pt)     = new HistoData("AK5Jet2Pt",     480,  0.,  480.);
  vhd.at(nAK5Jets)      = new HistoData("nAK5Jets",       9,  -0.5, 8.5);
  vhd.at(nAK5BJets)     = new HistoData("nAK5BJets",   9,  -0.5, 8.5);

  vhd.at(CA8Jet1Pt)     = new HistoData("CA8Jet1Pt",     480,  0.,  480.);
  vhd.at(CA8Jet2Pt)     = new HistoData("CA8Jet2Pt",     480,  0.,  480.);
  vhd.at(nCA8Jets)      = new HistoData("nCA8Jets",       9,  -0.5, 8.5);

  vhd.at(CAWJet1Pt)     = new HistoData("CAWJet1Pt",     480,  0.,  960.);
  vhd.at(CAWJet2Pt)     = new HistoData("CAWJet2Pt",     480,  0.,  960.);
  vhd.at(nCAWJets)      = new HistoData("nCAWJets",       9,  -0.5, 8.5);
  vhd.at(CAWBtag)       = new HistoData("CAWBtag",       100,  0.0, 1.0);

  vhd.at(CATopJet1Pt)   = new HistoData("CATopJet1Pt",   480,  0.,  960.);
  vhd.at(CATopJet2Pt)   = new HistoData("CATopJet2Pt",   480,  0.,  960.);
  vhd.at(nCATopJets)    = new HistoData("nCATopJets",     9,  -0.5, 8.5);
  vhd.at(CATopBtag)     = new HistoData("CATopBtag",       100,  0.0, 1.0);

  vhd.at(HT)        = new HistoData("HT", 500,  0.,  2500.);
  vhd.at(MET)       = new HistoData("MET", 200,  0., 1000.);
  vhd.at(sumPtL)    = new HistoData("sumPtL", 300,  0., 1500.);

  vhd.at(LepInvM)   = new HistoData("LepInvM",   500,  0.,  500.);
  vhd.at(AllLepPt)    = new HistoData("AllLepPt", 300,  0., 600.);
  vhd.at(dRlept)	= new HistoData("dRlept" , 100, 0., 5.);

  vhd.at(MinLepCATopdR) = new HistoData("MinLepCATopdR", 100,  0.,  5.);
  vhd.at(MinLepCAWdR)   = new HistoData("MinLepCAWdR",   100,  0.,  5.);
  vhd.at(MinLepCA8dR)   = new HistoData("MinLepCA8dR",   100,  0.,  5.);
  vhd.at(MinLepAK5dR)   = new HistoData("MinLepAK5dR",   100,  0.,  5.);
  vhd.at(MinLepJetdR)	= new HistoData("MinLepJetdR" , 100, 0., 5.);
  vhd.at(MinLepBJetdR)	= new HistoData("MinLepBJetdR" , 100, 0., 5.);

  vhd.at(leptonMETSum 	 )    = new HistoData("leptonMETSum", 500,  0., 2000.);
  vhd.at(leptonJetsSum   )    = new HistoData("leptonJetsSum", 500,  0., 2000.);
  vhd.at(jetsMETSum	 )    = new HistoData("jetsMETSum", 500,  0., 2000.);
  vhd.at(leptonJetsMETSum)    = new HistoData("leptonJetsMETSum", 500,  0., 2000.);

  vhd.at(minMlb)	= new HistoData("minMlb" , 300,  0.,  1500.);
  vhd.at(maxMlb)	= new HistoData("maxMlb" , 300,  0.,  1500.);
  vhd.at(Mlb)		= new HistoData("Mlb" , 300,  0.,  1500.);
  vhd.at(amwt)	        = new HistoData("amwt" , 1900, 100., 2000.);

  vhd.at(nCATopDau)       = new HistoData("nCATopDau",         8,  -0.5,  7.5);
  vhd.at(CATopMass)       = new HistoData("CATopMass",       240,   0.,  240.);
  vhd.at(CATopMinMass)    = new HistoData("CATopMinMass",    240,   0.,  240.);
  vhd.at(CAWMass)         = new HistoData("CAWMass",         240,   0.,  240.);

  vhd.at(nSelCATopDau)    = new HistoData("nSelCATopDau",      8,  -0.5,  7.5);
  vhd.at(SelCATopMass)    = new HistoData("SelCATopMass",    240,   0.,  240.);
  vhd.at(SelCATopMinMass) = new HistoData("SelCATopMinMass", 240,   0.,  240.);
  vhd.at(SelCAWMass)      = new HistoData("SelCAWMass",      240,   0.,  240.);

  vhd.at(CA8Mass)         = new HistoData("CA8Mass",         240,   0.,  240.);
  vhd.at(SelCA8Mass)      = new HistoData("SelCA8Mass",      240,   0.,  240.);

  vhd.at(AK5JetRCN)      = new HistoData("AK5JetRCN",      500,   0.,  100.);
  vhd.at(AK5JetnChHad)   = new HistoData("AK5JetnChHad",   100,   0.,  100.);
  vhd.at(AK5JetnNeuHad)  = new HistoData("AK5JetnNeuHad",  100,   0.,  100.);



  if (hname != "") hname = "_"+hname;

  goodHistos.resize(NHISTOS);

  for (unsigned int ih = 0; ih < NHISTOS; ih++){
//      cout << ih<<" "<<vhd.at(ih)->histoName<<" / "<<hname<<" / "<<sampleName_<<" / "<<suffix<<endl;
    goodHistos.at(ih) = new TH1F(vhd.at(ih)->histoName+hname+"_"+sampleName_+"_"+suffix, vhd.at(ih)->histoName+hname+"_"+sampleName_+"_"+suffix, 
				 vhd.at(ih)->nBins, vhd.at(ih)->min, vhd.at(ih)->max);
    goodHistos.at(ih)->Sumw2();
  }

  if (doFakes){
    Nt00Histos.resize(NHISTOS);
    Nt10Histos.resize(NHISTOS);
    Nt11Histos.resize(NHISTOS);
    Nt01Histos.resize(NHISTOS);

    for (unsigned int ih = 0; ih < NHISTOS; ih++){
      Nt00Histos.at(ih) = new TH1F(vhd.at(ih)->histoName+"Nt00"+hname+"_"+sampleName_+"_"+suffix, vhd.at(ih)->histoName+"Nt00"+hname+"_"+sampleName_+"_"+suffix, 
				   vhd.at(ih)->nBins, vhd.at(ih)->min, vhd.at(ih)->max);
      Nt00Histos.at(ih)->Sumw2();

      Nt01Histos.at(ih) = new TH1F(vhd.at(ih)->histoName+"Nt01"+hname+"_"+sampleName_+"_"+suffix, vhd.at(ih)->histoName+"Nt01"+hname+"_"+sampleName_+"_"+suffix, 
				   vhd.at(ih)->nBins, vhd.at(ih)->min, vhd.at(ih)->max);
      Nt01Histos.at(ih)->Sumw2();

      Nt10Histos.at(ih) = new TH1F(vhd.at(ih)->histoName+"Nt10"+hname+"_"+sampleName_+"_"+suffix, vhd.at(ih)->histoName+"Nt10"+hname+"_"+sampleName_+"_"+suffix, 
				   vhd.at(ih)->nBins, vhd.at(ih)->min, vhd.at(ih)->max);
      Nt10Histos.at(ih)->Sumw2();

      Nt11Histos.at(ih) = new TH1F(vhd.at(ih)->histoName+"Nt11"+hname+"_"+sampleName_+"_"+suffix, vhd.at(ih)->histoName+"Nt11"+hname+"_"+sampleName_+"_"+suffix, 
				   vhd.at(ih)->nBins, vhd.at(ih)->min, vhd.at(ih)->max);
      Nt11Histos.at(ih)->Sumw2();
    }

    }//closing the doFakes
}//closing function definition

void StandardHistos::fillAllHistos(TprimeEvent* teve, treetop* theTree, TString type, bool doFakes)
{
  if (not doFakes) fillVHD(teve, theTree, goodHistos);
  else{
   if(teve->theGoodLeptons.size()<3){
    if (type == "Nt00") fillVHD(teve, theTree, Nt00Histos);
    if (type == "Nt01") fillVHD(teve, theTree, Nt01Histos);
    if (type == "Nt10") fillVHD(teve, theTree, Nt10Histos);
    if (type == "Nt11") fillVHD(teve, theTree, Nt11Histos);
    }
  else if(teve->theGoodLeptons.size()>=3){
    if (type == "Nt0") fillVHD(teve, theTree, Nt00Histos);
    if (type == "Nt1") fillVHD(teve, theTree, Nt01Histos);
    if (type == "Nt2") fillVHD(teve, theTree, Nt10Histos);
    if (type == "Nt3") fillVHD(teve, theTree, Nt11Histos);
    }  
  }
}


void StandardHistos::fillVHD(TprimeEvent* teve, treetop* theTree, vector <TH1F*> myHistos){

  myHistos.at(hnEvents)->Fill(1, teve->weight);

  myHistos.at(hnPV)    ->Fill(theTree->nPV,      teve->weight);
  myHistos.at(hPUwgt)  ->Fill(teve->PUweight);

  myHistos.at(Lep1Pt) ->Fill(teve->theGoodLeptons[0]->lv.Pt(), teve->weight);
  myHistos.at(Lep2Pt) ->Fill(teve->theGoodLeptons[1]->lv.Pt(), teve->weight);
  if (theTree->goodLeptons.size()>2) myHistos.at(Lep3Pt)->Fill(theTree->goodLeptons[2]->lv.Pt(), teve->weight);

  //Fill histos split by lepton flavor
  if (teve->theGoodLeptons[0]->isMuon()){
    myHistos.at(MuPt) ->Fill(teve->theGoodLeptons[0]->lv.Pt(),  teve->weight);
    myHistos.at(MuEta)->Fill(teve->theGoodLeptons[0]->lv.Eta(), teve->weight);
  }
  else{
    myHistos.at(ElPt) ->Fill(teve->theGoodLeptons[0]->lv.Pt(),  teve->weight);
    myHistos.at(ElEta)->Fill(teve->theGoodLeptons[0]->lv.Eta(), teve->weight);
  }
  if (teve->theGoodLeptons[1]->isMuon()){
    myHistos.at(MuPt) ->Fill(teve->theGoodLeptons[1]->lv.Pt(),  teve->weight);
    myHistos.at(MuEta)->Fill(teve->theGoodLeptons[1]->lv.Eta(), teve->weight);
  }
  else{
    myHistos.at(ElPt) ->Fill(teve->theGoodLeptons[1]->lv.Pt(),  teve->weight);
    myHistos.at(ElEta)->Fill(teve->theGoodLeptons[1]->lv.Eta(), teve->weight);
  }

  int ie=0;
  int im=0;
  for (unsigned int i=0;i<teve->theGoodLeptons.size();++i){
    if (teve->theGoodLeptons[1]->isMuon()) ++im;
    else ++ie;
  }
  myHistos.at(nElec)  ->Fill(ie, teve->weight);
  myHistos.at(nMuon)  ->Fill(im, teve->weight);
  myHistos.at(nLept)  ->Fill(teve->theGoodLeptons.size(), teve->weight);

  if (teve->vGoodJets.size() > 0)  myHistos.at(Jet1Pt)->Fill(teve->vGoodJets.at(0)->lv.Pt(), teve->weight);
  if (teve->vGoodJets.size() > 1)  myHistos.at(Jet2Pt)->Fill(teve->vGoodJets.at(1)->lv.Pt(), teve->weight); 
  myHistos.at(nJets)  ->Fill(teve->vGoodJets.size(), teve->weight);
  myHistos.at(nBJets)  ->Fill(teve->nBTags(), teve->weight);
  if (teve->vGoodAK5Jets.size() > 0)  myHistos.at(AK5Jet1Pt)->Fill(teve->vGoodAK5Jets.at(0)->lv.Pt(), teve->weight);
  if (teve->vGoodAK5Jets.size() > 1)  myHistos.at(AK5Jet2Pt)->Fill(teve->vGoodAK5Jets.at(1)->lv.Pt(), teve->weight); 
  myHistos.at(nAK5Jets)->Fill(teve->vGoodAK5Jets.size(), teve->weight);
  myHistos.at(nAK5BJets)  ->Fill(teve->nAK5Tags(), teve->weight);
  for (unsigned int j = 0;j<teve->vGoodJets.size();++j ) {
    myHistos.at(AK5JetRCN)->    Fill(teve->vGoodJets.at(j)->jetRCN, teve->weight);
    myHistos.at(AK5JetnChHad)-> Fill(teve->vGoodJets.at(j)->jetnChHad, teve->weight);
    myHistos.at(AK5JetnNeuHad)->Fill(teve->vGoodJets.at(j)->jetnNeuHad, teve->weight);
  }

// //The selected CA8 jets
//   if (teve->vGoodCA8Jets.size() > 0)  myHistos.at(CA8Jet1Pt)->Fill(teve->vGoodCA8Jets.at(0)->lv.Pt(), teve->weight);
//   if (teve->vGoodCA8Jets.size() > 1)  myHistos.at(CA8Jet2Pt)->Fill(teve->vGoodCA8Jets.at(1)->lv.Pt(), teve->weight); 
//   myHistos.at(nCA8Jets)->Fill(teve->vGoodCA8Jets.size(), teve->weight);
//   for (unsigned int ui = 0; ui < teve->vGoodCA8Jets.size(); ui++){
//     myHistos.at(SelCA8Mass)->Fill(teve->vGoodCA8Jets.at(ui)->lv.M(), teve->weight);
//   }
// 
//   for (unsigned int ui = 0; ui < theTree->allCA8Jets.size(); ui++){
//     myHistos.at(CA8Mass)->Fill(theTree->allCA8Jets.at(ui)->lv.M(), teve->weight);
//   }
// 
// //The selected CAW jets
//   if (teve->vGoodCAWJets.size() > 0)  myHistos.at(CAWJet1Pt)->Fill(teve->vGoodCAWJets.at(0)->lv.Pt(), teve->weight);
//   if (teve->vGoodCAWJets.size() > 1)  myHistos.at(CAWJet2Pt)->Fill(teve->vGoodCAWJets.at(1)->lv.Pt(), teve->weight); 
//   myHistos.at(nCAWJets)->Fill(teve->vGoodCAWJets.size(), teve->weight);
//   for (unsigned int j = 0;j<teve->vGoodCAWJets.size();++j ) {
//     myHistos.at(CAWBtag)->Fill(teve->vGoodCAWJets.at(j)->csvDiscr, teve->weight);
//     myHistos.at(SelCAWMass)->Fill(teve->vGoodCAWJets.at(j)->lv.M(), teve->weight);
//   }
// 
// // All CAW jets
// 
//   for (unsigned int ui = 0; ui < theTree->allCAWJets.size(); ui++){
//     myHistos.at(CAWMass)->Fill(theTree->allCAWJets.at(ui)->lv.M(), teve->weight);
//   }
// 
// // Plots for the selected top jets:
// 
//   if (teve->vGoodCATopJets.size() > 0)  myHistos.at(CATopJet1Pt)->Fill(teve->vGoodCATopJets.at(0)->lv.Pt(), teve->weight);
//   if (teve->vGoodCATopJets.size() > 1)  myHistos.at(CATopJet2Pt)->Fill(teve->vGoodCATopJets.at(1)->lv.Pt(), teve->weight); 
//   myHistos.at(nCATopJets)->Fill(teve->vGoodCATopJets.size(), teve->weight);
// 
//   for (unsigned int j = 0;j<teve->vGoodCATopJets.size();++j ) {
//     myHistos.at(CATopBtag)->Fill(teve->vGoodCATopJets.at(j)->csvDiscr, teve->weight);
//     myHistos.at(nSelCATopDau)   ->Fill(teve->vGoodCATopJets.at(j)->nDaughters, teve->weight);
//     myHistos.at(SelCATopMass)   ->Fill(teve->vGoodCATopJets.at(j)->lv.M(), teve->weight);
//     myHistos.at(SelCATopMinMass)->Fill(teve->vGoodCATopJets.at(j)->minPairMass, teve->weight);
//   }
// 
// // Plots for all top jets:
// 
//   for (unsigned int ui = 0; ui < theTree->allCATopJets.size(); ui++){
//     myHistos.at(nCATopDau)   ->Fill(theTree->allCATopJets.at(ui)->nDaughters, teve->weight);
//     myHistos.at(CATopMass)   ->Fill(theTree->allCATopJets.at(ui)->lv.M(), teve->weight);
//     myHistos.at(CATopMinMass)->Fill(theTree->allCATopJets.at(ui)->minPairMass, teve->weight);     
//   }


//Other plots:

  myHistos.at(HT)   ->Fill(teve->ht(), teve->weight);
  myHistos.at(MET)  ->Fill(teve->met, teve->weight);

  myHistos.at(sumPtL)  ->Fill(teve->leptonScalSum(), teve->weight);

  myHistos.at(AllLepPt) ->Fill(teve->leptonSum().Pt(), teve->weight);
  if(theTree->goodLeptons.size()<3){
  myHistos.at(LepInvM)->Fill(teve->leptonSum().M(), teve->weight);
  }
  else if(theTree->goodLeptons.size()>=3){
  for(unsigned int il = 0; il<((theTree->goodLeptons.size())-1); il++){
    for(unsigned int ik = il+1;ik<theTree->goodLeptons.size();ik++) 
      myHistos.at(LepInvM)->Fill((teve->theGoodLeptons[il]->lv+teve->theGoodLeptons[ik]->lv).M(), teve->weight);
   }
  }
  myHistos.at(dRlept) ->Fill(teve->dRlepton(), teve->weight);

  myHistos.at(leptonMETSum    ) ->Fill(teve->leptonMETScalSum    (), teve->weight);
  myHistos.at(leptonJetsSum   ) ->Fill(teve->leptonJetsScalSum   (), teve->weight);
  myHistos.at(jetsMETSum      ) ->Fill(teve->jetsMETScalSum      (), teve->weight);
  myHistos.at(leptonJetsMETSum) ->Fill(teve->leptonJetsMETScalSum(), teve->weight);

  myHistos.at(MinLepAK5dR)->Fill(teve->mindrJetLepton(teve->vGoodAK5Jets), teve->weight);
  myHistos.at(MinLepCA8dR)->Fill(teve->mindrJetLepton(teve->vGoodCA8Jets), teve->weight);
  myHistos.at(MinLepCAWdR)->Fill(teve->mindrJetLepton(teve->vGoodCAWJets), teve->weight);
  myHistos.at(MinLepCATopdR)->Fill(teve->mindrJetLepton(teve->vGoodCATopJets), teve->weight);
  myHistos.at(MinLepJetdR)->Fill(teve->mindrJetLepton(teve->vGoodJets), teve->weight);

  vector<double> m(teve->mlb());
  if (m.size()>0) {
    myHistos.at(minMlb)  ->Fill(teve->minMlb(), teve->weight);
    for (unsigned int j = 0;j<m.size();++j )     myHistos.at(Mlb)  ->Fill( m[j], teve->weight);
    myHistos.at(maxMlb)  ->Fill(*(max_element(m.begin() ,m.end())), teve->weight);
  }
  myHistos.at(MinLepBJetdR)  ->Fill( teve->mindrBJetLepton(), teve->weight);
//   vector<double> dr(teve->drBJetLepton());
//   for (unsigned int j = 0;j<dr.size();++j )     myHistos.at(dR)  ->Fill( dr[j], teve->weight);

  if (teve->amwtPeakWeight>0.) myHistos.at(amwt)  ->Fill(teve->amwtPeakMass, teve->weight);

}


void StandardHistos::scaleAllHistos(bool doFakes){
  for (unsigned int ih = 0; ih < NHISTOS; ih++){
    goodHistos.at(ih)->Scale(scale_);

    if (doFakes){
       Nt00Histos.at(ih)->Scale(scale_);
       Nt01Histos.at(ih)->Scale(scale_);
       Nt10Histos.at(ih)->Scale(scale_);
       Nt11Histos.at(ih)->Scale(scale_);
    }
  }

}

void StandardHistos::getTotFakes(double f1, double f2, double p1, double p2){
  for (unsigned int ih = 0; ih < NHISTOS; ih++){

    double det = 1.0 / ((p1 - f1)*(p2 - f2)); 

    TH1F* hNff = (TH1F*)Nt00Histos.at(ih)->Clone("hNff");
    hNff->Scale(p1*p2);
    hNff->Add(Nt10Histos.at(ih), -1*(1-p1)*p2);
    hNff->Add(Nt01Histos.at(ih), -1*p1*(1-p2));
    hNff->Add(Nt11Histos.at(ih), (1-p1)*(1-p2));
    hNff->Scale(det*f1*f2);

    TH1F* hNpf = (TH1F*)Nt00Histos.at(ih)->Clone("hNpf");
    hNpf->Scale(-1*f1*p2);
    hNpf->Add(Nt10Histos.at(ih), (1-f1)*p2);
    hNpf->Add(Nt01Histos.at(ih), f1*(1-p2));
    hNpf->Add(Nt11Histos.at(ih), -1*(1-f1)*(1-p2));
    hNpf->Scale(det*p1*f2);

    TH1F* hNfp = (TH1F*)Nt00Histos.at(ih)->Clone("hNfp");
    hNfp->Scale(-1*p1*f2);
    hNfp->Add(Nt10Histos.at(ih), (1-p1)*f2);
    hNfp->Add(Nt01Histos.at(ih), p1*(1-f2));
    hNfp->Add(Nt11Histos.at(ih), -1*(1-p1)*(1-f2));
    hNfp->Scale(det*f1*p2);

    goodHistos.at(ih)->Add(hNff);
    goodHistos.at(ih)->Add(hNpf);
    goodHistos.at(ih)->Add(hNfp);
    
    delete hNff;
    delete hNpf;
    delete hNfp;
  }
}

void StandardHistos::getTotFakesTRI(double f, double p){
  for (unsigned int ih = 0; ih < NHISTOS; ih++){

    double det = 1.0/((p-f)*(p-f)*(p-f));

    TH1F* hNfff = (TH1F*)Nt00Histos.at(ih)->Clone("hNfff");
    hNfff->Scale(p*p*p);
    hNfff->Add(Nt01Histos.at(ih), -p*p*(1-p));
    hNfff->Add(Nt10Histos.at(ih), p*(1-p)*(1-p));
    hNfff->Add(Nt11Histos.at(ih), -(1-p)*(1-p)*(1-p));
    hNfff->Scale(det);

    
    TH1F* hNffp = (TH1F*)Nt00Histos.at(ih)->Clone("hNffp");
    hNffp->Scale(-3*p*p*f);
    hNffp->Add(Nt01Histos.at(ih), 2*p*f*(1-p)+p*p*(1-f));
    hNffp->Add(Nt10Histos.at(ih), -(f*(1-p)*(1-p) + 2*p*(1-p)*(1-f)));
    hNffp->Add(Nt11Histos.at(ih), 3*(1-p)*(1-p)*(1-f));
    hNffp->Scale(det);

    TH1F* hNfpp = (TH1F*)Nt00Histos.at(ih)->Clone("hNfpp");
    hNfpp->Scale(3*p*f*f);
    hNfpp->Add(Nt01Histos.at(ih), -(f*f*(1-p) + 2*p*f*(1-f)));
    hNfpp->Add(Nt10Histos.at(ih), (2*f*(1-p)*(1-f) + p*(1-f)*(1-f)));
    hNfpp->Add(Nt11Histos.at(ih), -3*(1-p)*(1-f)*(1-f));
    hNfpp->Scale(det);

    //cout << "fakes sum = " << -(1-p)*(1-p)*(1-p)*det*15 + 3*(1-p)*(1-p)*(1-f)*det*15 - 3*(1-p)*(1-f)*(1-f)*det*15 << endl;
    
    goodHistos.at(ih)->Add(hNfff, f*f*f);
    goodHistos.at(ih)->Add(hNffp, f*f*p);
    goodHistos.at(ih)->Add(hNfpp, f*p*p);

    delete hNfff;
    delete hNffp;
    delete hNfpp;
  }
}

