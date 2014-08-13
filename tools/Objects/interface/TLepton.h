#ifndef ObjectID_cxx
#define ObjectID_cxx

#include <TLorentzVector.h>
#include <cmath>
#include <iostream>



using namespace std;

//Particle masses and acceptable deltas
//Do they belong here or somewhere else?
const double MTOP  = 173.5;
const double DMTOP = 30.0;
const double MW    = 80.385;
const double DMW   = 20.0;

class MCLepton{

public:
  
  TLorentzVector lv;

  double pt;
  double eta;
  double phi;
  double energy;
  int pdgID;

  void setLV() {
    lv.SetPxPyPzE(pt*cos(phi), pt*sin(phi), pt*TMath::SinH(eta), energy);
  }

};

class TLepton{

public:
  
  TLepton(int chargeTemp, double ptTemp, double etaTemp, double phiTemp, double energyTemp, 
	double relIsoTemp, double ipXYTemp, double ipzTemp,
	int muonTypeTempm,
	int innerTrackLayersTemp, double nChi2Temp, int MuonHitsTemp, int PixelHitsTemp, int nMatchTemp, int muQual);

  TLepton(int chargeTemp, double ptTemp, double etaTemp, double phiTemp, double energyTemp, 
	double relIsoTemp, double ipXYTemp, double ipzTemp,
	int eEEBtypes, int elQuality, int chargeConsistencyTemp, int notConvTemp,
	double elDeta, double elDphi, 
	double elSihih, double elHoE, double elD0, double elOoemoop,
	int elMHits, int elVtxFitConv);

  bool isElectron() { return isElectron_;}
  bool isMuon(){ return isMuon_;}
  
  TLorentzVector lv;
  int charge;
  double pt;
  double eta;
  double phi;
  double energy;
  double relIso;
  double ipXY, ipZ;
  bool   isElectron_, isMuon_;

  bool qualLJ_, qualDL_;
  //Electrons only
  bool notEEEBGap_, isEE_, isEB_;
  bool elQualityT_, elQualityM_, elQualityL_;
  int    chargeConsistency;
  int    notConv;
  double deltaEtaSuperClusterTrackAtVtx, deltaPhiSuperClusterTrackAtVtx;
  double sigmaIetaIeta, hadronicOverEm, dB, ooemoop; 
  int trackerExpectedHitsInner;
  int passConversionVeto;


  //Muons only
  int    GlobalMuon;
  int    TrackerMuon;
  int    innerTrackLayers;
  double nChi2;
  int    MuonHits; // number of valid muon hits
  int    PixelHits;
  int    nMatch;
//   double relPtUnc;
  //mc
  MCLepton mother[10];
  MCLepton matchedLepton;
  bool matched;
  double gen_reco_DR;

  int origin() const;
  bool goodElectron(double minPt);
  bool looseElectron(double minPt);
  bool goodMuon(double minPt);
  bool looseMuon(double minPt);
  bool lntLepton(double minPt);
  bool goodLepton(double minPt);
  void printLepton(bool fullMatching = false);
  bool isDileptonLepton(double minPt);
  bool isLeptJetLepton(double minPt);



  inline int hundred (int i) const
  {
    int j = (int) i % 1000;
    return (int) j/100;
  }
  inline int thousand (int i) const
  {
    int j = (int) i % 10000;
    return (int) j/1000;
  }

private :

  void setLV() ;
  
  double MIN_LEP_PT_FR ;

  double MAX_EL_ISO;
  double MAX_MU_ISO;

  double MAX_EL_ISO_LOOSE;
  double MAX_MU_ISO_LOOSE;

  double MAX_EL_ETA;
  double MAX_MU_ETA;

  TLepton(){defaults();}
  void defaults();
//   void initLepton();
};

class TJet{

public:
  
  TJet(double ptTemp, double etaTemp, double phiTemp, double energyTemp, bool bTagTemp=-1.,
  	int jetnChHadTemp=0, int jetnNeuHadTemp=0, double jetRCNTemp=0.);     //AK5 jets
  TJet(double ptTemp, double etaTemp, double phiTemp, double energyTemp, bool CA8DauTemp, 
       int motherIndexTemp);                                                                 //For daughters
  TJet(double ptTemp, double etaTemp, double phiTemp, double energyTemp, int indexTemp, 
       int nDaughtersTemp, double massTemp, double btagTemp);                                //CAW jets
  TJet(double ptTemp, double etaTemp, double phiTemp, double energyTemp, int indexTemp, 
       int nDaughtersTemp, double topMassTemp, double minPairMassTemp, double btagTemp);     //CATop jets


  TLorentzVector lv;

  double pt;
  double eta;
  double phi;
  double energy;
  
  int csvMedium;
  int partonFlavour;
  int jetnChHad, jetnNeuHad;
  double jetRCN;

  //For jet substructure
  bool   CA8Dau;
  int    index;
  int    nDaughters;
  int    motherIndex;
  double mass;
  double topMass;
  double minPairMass;

  double csvDiscr;

  void setLV();
  void printJet();
  bool checkJetLeptondR(const TLepton & lep1, const TLepton & lep2);
  bool checkJetLeptondR(const TLepton & lep1);
  bool checkJetJetdR(const TJet & jet1);
private:
  double MIN_JET_PT;
  double MAX_JET_ETA;

  TJet();
  void defaults();
};




#endif
