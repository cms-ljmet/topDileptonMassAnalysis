#ifndef TprimeEvent_h
#define TprimeEvent_h
#include <TString.h>
#include "TLepton.h"

class TprimeEvent{

public:

  TString suffix;

  unsigned int isam;
  int channel;

  double weight, PUweight;
  int    nPV;
 
//   TLepton* lepton1;
//   TLepton* lepton2;
  vector <TLepton *> theGoodLeptons;
  vector <TLepton *> additionalLeptons;

  double leptonScalSum();
  double leptonMETScalSum();
  double leptonJetsScalSum();
  double jetsMETScalSum();
  double leptonJetsMETScalSum();

  TLorentzVector leptonSum();
  TLorentzVector jetVectSum();
  TLorentzVector leptonMETVectSum();
  TLorentzVector leptonJetsVectSum();
  TLorentzVector jetsMETVectSum();
  TLorentzVector leptonJetsMETVectSum();

  vector <TJet*>        vGoodCATopJets;
  vector <TJet*>        vGoodCA8Jets;
  vector <TJet*>        vGoodCAWJets;
  vector <TJet*>        vGoodAK5Jets;
  vector <TJet*>        vGoodJets;
  vector <unsigned int>   vGoodCsvs;

  double ht();
  double met;
  double met_phi;
  double mt;
  
  int nBTags();
  int nAK5Tags();
  
  TprimeEvent()  { clearEvent(); }

  void clearEvent();
  void printEvent();

  double dRlepton() const;

  vector<double> drJetLepton (const vector <TJet*> & theJets) const;
  double mindrJetLepton(const vector <TJet*> & theJets) const;

  vector<double> drBJetLepton () const;
  double mindrBJetLepton() const;

  vector<double> mlb () const;
  double minMlb() const;
  double maxMlb() const;
  
  vector <TJet*> selectJetsForMass();
  float amwtPeakMass, amwtPeakWeight;

};

#endif
