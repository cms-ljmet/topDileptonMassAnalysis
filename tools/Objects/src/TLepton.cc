#include "../interface/TLepton.h"

TLepton::TLepton(int chargeTemp, double ptTemp, double etaTemp, double phiTemp, double energyTemp, 
	double relIsoTemp, double ipXYTemp, double ipzTemp, int muonTypeTemp, 
	int innerTrackLayersTemp, double nChi2Temp, int MuonHitsTemp, int PixelHitsTemp, int nMatchTemp,
	int muQual) :
  charge(chargeTemp), pt(ptTemp), eta(etaTemp), phi(phiTemp), energy(energyTemp),
  relIso(relIsoTemp), ipXY(ipXYTemp), ipZ(ipzTemp), 
  isElectron_(false), isMuon_(true),
  innerTrackLayers(innerTrackLayersTemp), 
  nChi2(nChi2Temp), MuonHits(MuonHitsTemp), PixelHits(PixelHitsTemp), 
  nMatch(nMatchTemp)
{

  defaults();
  if (muonTypeTemp==8) {
    GlobalMuon = true;
    TrackerMuon= true;
  } else {
    GlobalMuon = (muonTypeTemp&0x4)>>1;
    TrackerMuon= muonTypeTemp&0x1;
  }
  qualLJ_ = (muQual & 0x2) >> 1; 
  qualDL_ = (muQual & 0x1);
  setLV();
}

TLepton::TLepton(int chargeTemp, double ptTemp, double etaTemp, double phiTemp, double energyTemp, 
	double relIsoTemp, double ipXYTemp, double ipzTemp,
	int eEEBtypes, int elQuality, int chargeConsistencyTemp, int notConvTemp,
	double elDeta, double elDphi,
	double elSihih, double elHoE, double elD0, double elOoemoop,
	int elMHits, int elVtxFitConv):
  charge(chargeTemp), pt(ptTemp), eta(etaTemp), phi(phiTemp), energy(energyTemp),
  relIso(relIsoTemp), ipXY(ipXYTemp), ipZ(ipzTemp), 
  isElectron_(true), isMuon_(false),
  chargeConsistency(chargeConsistencyTemp), notConv(notConvTemp),
  deltaEtaSuperClusterTrackAtVtx(elDeta), deltaPhiSuperClusterTrackAtVtx(elDphi),
  sigmaIetaIeta(elSihih), hadronicOverEm(elHoE), dB(elD0), ooemoop(elOoemoop),
  trackerExpectedHitsInner(elMHits), passConversionVeto(elVtxFitConv)
{
  defaults();
  notEEEBGap_ = (eEEBtypes & 0x4) >> 2;
  isEE_       = (eEEBtypes & 0x2) >> 1;
  isEB_       =  eEEBtypes & 0x1;

  elQualityT_ = (elQuality & 0x4) >> 2;
  elQualityM_ = (elQuality & 0x2) >> 1;
  elQualityL_ =  elQuality & 0x1;
  qualLJ_ = (elQuality & 0x10) >> 4; 
  qualDL_ = (elQuality & 0x8) >> 3;
  setLV();
}

bool TLepton::isDileptonLepton(double minPt)
{
  if (lv.Pt() < minPt)                            return false;
  return (qualDL_);
}
bool TLepton::isLeptJetLepton(double minPt)
{
  if (lv.Pt() < minPt)                            return false;
  return (qualLJ_);
}

int TLepton::origin() const
{
  if (!matched) return 0;
  int i = 0;
  while (i<10) {
    if (abs(mother[i].pdgID)==24) return 24; // The mother is a W
    if ((abs(mother[i].pdgID)!=11) && (abs(mother[i].pdgID)!=13) && (abs(mother[i].pdgID)!=15)) break; 
    ++i; // The mother is a lepton, continue checking whether the next mother is a W.
  }
  // So it is not a W
  for (i = 0;i<10;++i ) {
    if ((thousand(fabs(mother[i].pdgID))== 5) || (hundred(fabs(mother[i].pdgID))== 5)) return 5;
  }
  for (i = 0;i<10;++i ) {
    if ((thousand(fabs(mother[i].pdgID))== 4) || (hundred(fabs(mother[i].pdgID))== 4)) return 4;
  }
  for (i = 0;i<10;++i ) {
    int t = thousand(fabs(mother[i].pdgID));
    int h = hundred(fabs(mother[i].pdgID));
    if ((t==1)||(t==2)|| (t==3) || (h==1)||(h==2)|| (h==3)) return 1;
  }
  return -1;
}

void TLepton::setLV() {
  lv.SetPxPyPzE(pt*cos(phi), pt*sin(phi), pt*TMath::SinH(eta), energy);
//     if (matched) {
//       matchedLepton.setLV();
//       for (int i = 0;i<10;++i ) if (mother[i].pdgID!=0) mother[i].setLV();
//     }

}

bool TLepton::looseElectron(double minPt){
  //Common to good and loose
  if (isMuon_)                                    return false;
  if (lv.Pt() < minPt)                            return false;
  if (fabs(lv.Eta()) > MAX_EL_ETA)                return false;
  if (trackerExpectedHitsInner)                   return false; 
  if (not chargeConsistency)                      return false; 
  //Some are common, some are not
  if (isEB_){
    if (deltaEtaSuperClusterTrackAtVtx > 7.0e-03) return false;
    if (deltaPhiSuperClusterTrackAtVtx > 1.5e-01) return false; 
    if (sigmaIetaIeta                  > 1.0e-02) return false; //Common
    if (hadronicOverEm                 > 1.2e-01) return false; //Common
  }
  if (isEE_){
    if (deltaEtaSuperClusterTrackAtVtx > 9.0e-03) return false;
    if (deltaPhiSuperClusterTrackAtVtx > 1.0e-01) return false;
    if (sigmaIetaIeta                  > 3.0e-02) return false; //Common
    if (hadronicOverEm                 > 1.0e-01) return false; //Common
  }
  if (ooemoop                          > 5.0e-02) return false; //Common
  if (relIso > MAX_EL_ISO_LOOSE)                  return false;

  return true;
}

bool TLepton::goodElectron(double minPt){
  if (not looseElectron(minPt))                   return false;

  if (isEB_){
    if (deltaEtaSuperClusterTrackAtVtx > 4.0e-03) return false;
    if (deltaPhiSuperClusterTrackAtVtx > 3.0e-02) return false; 
  }
  if (isEE_){
    if (deltaEtaSuperClusterTrackAtVtx > 5.0e-03) return false;
    if (deltaPhiSuperClusterTrackAtVtx > 2.0e-02) return false;
  }
  if (fabs(ipZ) > 0.1 )                           return false;
  if (fabs(dB) > 0.02 )                           return false;
  if (not passConversionVeto)                     return false;
  if (relIso > MAX_EL_ISO)                        return false;  

  return true;
}

bool TLepton::looseMuon(double minPt){
  //Common to good and loose
  if (isElectron_)                 return false;
  if (lv.Pt() < minPt)             return false;
  if (fabs(lv.Eta()) > MAX_MU_ETA) return false;
  if (not GlobalMuon)              return false;
  if (not TrackerMuon)             return false;
  if (MuonHits < 1)                return false;
  if (PixelHits < 1)               return false;
  //Not common to good and loose
  if (nChi2 >= 50)                 return false;
  if (fabs(ipXY) > 2)              return false;
  if (relIso > MAX_MU_ISO_LOOSE)   return false;

  return true;
}

bool TLepton::goodMuon(double minPt){
  if (not looseMuon(minPt))         return false;

  if (nChi2 >= 10)                  return false;
  if (fabs(ipXY) > 0.2)             return false;
  if (fabs(ipZ)  > 0.5)             return false;
  if (innerTrackLayers <= 5)        return false;
  if (nMatch < 2)                   return false;
  if (relIso > MAX_MU_ISO)          return false;

  return true;
}

//For measurement of fake contribution
bool TLepton::lntLepton(double minPt){
  if (isElectron_ and looseElectron(minPt) and not goodElectron(minPt)) return true;
  if (isMuon_     and looseMuon(minPt)     and not goodMuon(minPt))     return true;
  return false;
}

bool TLepton::goodLepton(double minPt){
  if (goodElectron(minPt)) return true;
  if (goodMuon(minPt))     return true;
  return false;
}

void TLepton::printLepton(bool fullMatching){

  if (lv.Pt() < 1){
    cout<<"Lepton with pT < 1"<<endl;
//      return;
  }

  cout<<"pt: "<<lv.Pt()<<" eta: "<<lv.Eta()<<" phi: "<<lv.Phi()<<endl;
  cout<<"charge:            "<<charge<<endl;
  cout<<"iso:               "<<relIso<<endl;
  cout<<"ipXY / ipZ:        "<<ipXY<<" / "<< ipZ<<  endl;

  if (isElectron_){
    //Electrons only
    cout<<"EE/EB/notGap:        "<<isEE_<<" / "<< isEB_<<" / "<<  notEEEBGap_<<endl;
    cout<<"Cut based-id L/M/T   "<<elQualityL_<<" / "<< elQualityM_<<" / "<<  elQualityT_<<endl;

//     cout<<"hyperTightID:      "<<hyperTightID<<endl;
    cout<<"chargeConsistency:   "<<chargeConsistency<<endl;
    cout<<"Conversion veto:     "<<passConversionVeto<<endl;
    cout<<"trackerExpectedHitsInner:   "<<trackerExpectedHitsInner<<endl;
    cout<<"chargeConsistency:   "<<chargeConsistency<<endl;
    cout<<"delta Eta/Phi SuperClusterTrackAtVtx: "<<
	deltaEtaSuperClusterTrackAtVtx<<" / "<< deltaPhiSuperClusterTrackAtVtx<<endl;
    cout<<"sigmaIetaIeta:   "<<sigmaIetaIeta<<endl;
    cout<<"hadronicOverEm:   "<<hadronicOverEm<<endl;
    cout<<"dB:   "<<dB<<endl;
    cout<<"ooemoop:   "<<ooemoop<<endl;

  }
  else{
    //Muons only
    cout<<"Global / Tracker   "<<GlobalMuon<<" / "<< TrackerMuon<<endl;
    cout<<"innerTrackLayers:    "<<innerTrackLayers<<endl;
    cout<<"NormalizedChi2:    "<<nChi2<<endl;
    cout<<"MuonHits:          "<<MuonHits<<endl;
    cout<<"PixelHits:         "<<PixelHits<<endl;
    cout<<"nMatch:            "<<nMatch<<endl;

    if (goodMuon(30.)) cout << "Good muon (at 30 GeV...)\n";
    else if (looseMuon(30.)) cout << "Loose muon (at 30 GeV...)\n";
  }

  cout<<"MC match:          "<<matched<<" "<< origin()<<endl;
  if (fullMatching && matched) {
    cout<<"Matched lepton pt: "<<matchedLepton.lv.Pt()<<" eta: "<<matchedLepton.lv.Eta()
        <<" phi: "<<matchedLepton.lv.Phi()<<endl;
    cout <<"Mothers\n";
    for (int i = 0;i<10;++i )
      if (mother[i].pdgID!=0) cout<<i<<" : pdgID: "<< mother[i].pdgID<<" pt: "
	<< mother[i].lv.Pt()<<" eta: "<<mother[i].lv.Eta()
        <<" phi: "<<mother[i].lv.Phi()<<endl;
  }
  cout<<endl;
}

void TLepton::defaults(){
  MIN_LEP_PT_FR  = 25.;

//   MAX_EL_ISO  = 0.15;
//   MAX_MU_ISO  = 0.20;

  MAX_EL_ISO  = 0.10;
  MAX_MU_ISO  = 0.12;

  MAX_EL_ISO_LOOSE  = 0.60;
  MAX_MU_ISO_LOOSE  = 0.40;

  MAX_EL_ETA  = 2.4;
  MAX_MU_ETA  = 2.4;
}

// void TLepton::initLepton(){
//   lv.SetPxPyPzE(0,0,0,0);
//   relIso = SEN_HIGH;
//   ipXY = SEN_HIGH;
//   ipZ = SEN_HIGH;
//   charge = 0;
// 
//   //Electrons only
//   notEEEBGap_ = true;
//   chargeConsistency = false;
//   notConv = false;
// 
//   //Muons only
//   GlobalMuon = false;
//   TrackerMuon = false;
//   innerTrackLayers = SEN_LOW;
//   nChi2 = SEN_HIGH;
//   MuonHits = SEN_LOW;
//   PixelHits = SEN_LOW;
//   nMatch = SEN_LOW;
// //     relPtUnc = SEN_HIGH;
// }

TJet::TJet(double ptTemp, double etaTemp, double phiTemp, double energyTemp, bool bTagTemp,
  	int jetnChHadTemp, int jetnNeuHadTemp, double jetRCNTemp) :
  pt(ptTemp), eta(etaTemp), phi(phiTemp), energy(energyTemp), csvMedium(bTagTemp),
  jetnChHad(jetnChHadTemp), jetnNeuHad(jetnNeuHadTemp), jetRCN(jetRCNTemp)
{
  setLV();
  defaults();
}

TJet::TJet(double ptTemp, double etaTemp, double phiTemp, double energyTemp, bool CA8DauTemp, 
	   int motherIndexTemp) :
  pt(ptTemp), eta(etaTemp), phi(phiTemp), energy(energyTemp), CA8Dau(CA8DauTemp), 
  motherIndex(motherIndexTemp)
{
  setLV();
  defaults();
}

TJet::TJet(double ptTemp, double etaTemp, double phiTemp, double energyTemp, int indexTemp, 
	   int nDaughtersTemp, double massTemp, double btagTemp) :
  pt(ptTemp), eta(etaTemp), phi(phiTemp), energy(energyTemp), index(indexTemp), 
  nDaughters(nDaughtersTemp), mass(massTemp), csvDiscr(btagTemp)
{
  setLV();
  defaults();
}

TJet::TJet(double ptTemp, double etaTemp, double phiTemp, double energyTemp, int indexTemp, 
	   int nDaughtersTemp, double topMassTemp, double minPairMassTemp, double btagTemp) :
  pt(ptTemp), eta(etaTemp), phi(phiTemp), energy(energyTemp), index(indexTemp), 
  nDaughters(nDaughtersTemp), topMass(topMassTemp), minPairMass(minPairMassTemp),
  csvDiscr(btagTemp)
{
  setLV();
  defaults();
}


TJet::TJet()
{
  defaults();
}

void TJet::setLV() {
  lv.SetPxPyPzE(pt*cos(phi), pt*sin(phi), pt*TMath::SinH(eta), energy);
}

void TJet::printJet(){

  if (lv.Pt() < 1){
    cout<<"Jet with pT < 1"<<endl;
//      return;
  }

  cout<<"pt: "<<lv.Pt()<<" eta: "<<lv.Eta()<<" phi: "<<lv.Phi()<<" energy: "<<lv.E()<<endl;
  cout<<"CSV:       "<<csvMedium  <<endl;
  cout << jetRCN<< " "<< jetnChHad<< " "<< jetnNeuHad<<endl;
  cout<<endl;
}


bool TJet::checkJetLeptondR(const TLepton & lep1, const TLepton & lep2){

  bool dRcut = true; 

  if (lv.Pt() < MIN_JET_PT) return false; 

  if (lep1.lv.DeltaR(lv) < 0.3) dRcut = false;

  if (lep2.lv.DeltaR(lv) < 0.3) dRcut = false;

  return dRcut;
}

bool TJet::checkJetLeptondR(const TLepton & lep1){

  bool dRcut = true; 

  if ((lv.Pt() < MIN_JET_PT) || fabs(lv.Eta()) > MAX_JET_ETA) return false; 

  if (lep1.lv.DeltaR(lv) < 0.3) dRcut = false;

  return dRcut;
}

bool TJet::checkJetJetdR(const TJet & jet1){

  bool dRcut = true; 

  if ((lv.Pt() < MIN_JET_PT) || fabs(lv.Eta()) > MAX_JET_ETA) return false; 

  if (jet1.lv.Pt() > MIN_JET_PT){
    if (jet1.lv.DeltaR(lv) < 0.8) dRcut = false;
  }

  return dRcut;
}

void TJet::defaults(){
  MIN_JET_PT  = 30.;
  MAX_JET_ETA = 2.4;
}

