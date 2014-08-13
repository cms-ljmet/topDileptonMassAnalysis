#ifndef StandardHistos_cxx
#define StandardHistos_cxx

#include "TLorentzVector.h"
#include "TVector2.h"
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TString.h>

#include <vector>
#include <algorithm>

#include "../interface/TprimeEvent.h"
#include "../interface/treetop.h"

// const unsigned int     = 63;

enum Histos_t {Lep1Pt,  Lep2Pt,  Lep3Pt,			// 3
	       ElPt,    MuPt,    ElEta, MuEta, 			// 4	7
	       nElec, nMuon, nLept,				// 3	10
	       Jet1Pt,  Jet2Pt,  nJets, nBJets,			// 4	14
	       AK5Jet1Pt,   AK5Jet2Pt,   nAK5Jets, nAK5BJets,   // 4	18
	       CAWJet1Pt,   CAWJet2Pt,   nCAWJets, CAWBtag,     // 4	22
	       CATopJet1Pt, CATopJet2Pt, nCATopJets, CATopBtag, // 4	26
	       HT, MET, AllLepPt,     				// 3	29
	       LepInvM, sumPtL,  dRlept, 			// 3	32
	       MinLepAK5dR, MinLepCAWdR, MinLepCATopdR, MinLepJetdR,//4 36
	       MinLepBJetdR, minMlb, Mlb,			//3	39
	       amwt, hnPV, hnEvents, hPUwgt,			//4	43
	       nSelCATopDau, SelCATopMass, SelCATopMinMass,	//3	46
	       nCATopDau, CATopMass, CATopMinMass,		//3	49
	       CA8Jet1Pt,   CA8Jet2Pt,   nCA8Jets, MinLepCA8dR, //4	53
	       CAWMass, SelCAWMass, CA8Mass, SelCA8Mass,	//4	57
	       leptonMETSum, leptonJetsSum, jetsMETSum,		//3	60
	       leptonJetsMETSum, maxMlb,			//2	62
	       AK5JetRCN, AK5JetnChHad, AK5JetnNeuHad,		//3	65
	       NHISTOS


};

const unsigned int NHISTOS_FR = 9;

enum Histos_FR_t {Lep1Pt_FR, Jet1Pt_FR, nJets_FR, nEvents_FR,
		  nPV_FR,    mJL_FR,    dRJL_FR,  MET_FR,
		  MT_FR
};



class HistoData{

public:

  TString histoName;
  int     nBins;
  double  min;
  double  max;

  HistoData();

  HistoData(TString hn, int nb, double mn, double mx){
    histoName = hn;
    nBins     = nb;
    min       = mn;
    max       = mx;
  }

};

class StandardHistos{

public:

  StandardHistos(const TString & aName, const TString & aSuffix, const TString & sampleName,
  	const double scale, bool doFakes) :
	hname(aName), suffix(aSuffix), sampleName_(sampleName), scale_(scale)
  {
    initAllHistos(doFakes);
  }

  void fillVHD(TprimeEvent* teve, treetop* theTree, vector <TH1F*> myHistos);

  void fillAllHistos(TprimeEvent* teve, treetop* theTree, TString type="", bool doFakes=false);

  void scaleAllHistos(bool doFakes=false);

  void getTotFakes(double f1, double f2, double p1, double p2);

  void getTotFakesTRI(double f, double p);

  double scale() {return scale_;}

private:
  TString hname;
  TString suffix;
  TString sampleName_;
  double scale_;
  vector <TH1F*> goodHistos;
  vector <TH1F*> Nt00Histos;
  vector <TH1F*> Nt01Histos;
  vector <TH1F*> Nt10Histos;
  vector <TH1F*> Nt11Histos;
  void initAllHistos(bool doFakes=false);

};


// class StandardHistosFR{
// 
// private:
//   StandardHistosFR(){}
// public:
// 
//   StandardHistosFR(const TString & aName, const TString & aSuffix, const TString & aLabel, const double scale) :
// 	hname(aName), suffix(aSuffix), outLabel(aLabel), scale_(scale)
//   {
//     initAllHistos();
//   }
// 
//   TString hname;
// 
//   TString suffix;
//   TString outLabel;
//   double scale_;
// 
//   vector <TH1F*> looseHistos;
//   vector <TH1F*> tightHistos;
//   vector <TH1F*> ratioHistos;
// 
//   void initAllHistos(){
// 
//     vector <HistoData*> vhd(NHISTOS_FR);
//     vhd.at(Lep1Pt_FR)   = new HistoData("Lep1Pt",   60,  0.,  240.);
//     vhd.at(Jet1Pt_FR)   = new HistoData("Jet1Pt",   60,  0.,  420.);
//     vhd.at(nJets_FR)    = new HistoData("nJets",    9,  -0.5, 8.5);
//     vhd.at(nEvents_FR)  = new HistoData("nEvents",  1,   0.5, 1.5);
//     vhd.at(nPV_FR)      = new HistoData("nPV",      30, -0.5, 29.5);
//     vhd.at(mJL_FR)      = new HistoData("mJL",      60,  0.,  240.);
//     vhd.at(dRJL_FR)     = new HistoData("dRJL",     50,  0.,  10);
//     vhd.at(MET_FR)      = new HistoData("MET",      50,  0.,  250);
//     vhd.at(MT_FR)       = new HistoData("MT",       50,  0.,  250);
// 
//     looseHistos.resize(NHISTOS_FR);
//     tightHistos.resize(NHISTOS_FR);
//     ratioHistos.resize(NHISTOS_FR);
// 
//     for (unsigned int ih = 0; ih < NHISTOS_FR; ih++){
//       looseHistos.at(ih) = new TH1F(vhd.at(ih)->histoName+"_loose_"+outLabel+"_"+suffix, vhd.at(ih)->histoName+"_loose_"+outLabel+"_"+suffix,
// 				    vhd.at(ih)->nBins, vhd.at(ih)->min, vhd.at(ih)->max);
//       looseHistos.at(ih)->Sumw2();
// 
//       tightHistos.at(ih) = new TH1F(vhd.at(ih)->histoName+"_tight_"+outLabel+"_"+suffix, vhd.at(ih)->histoName+"_tight_"+outLabel+"_"+suffix,
// 				    vhd.at(ih)->nBins, vhd.at(ih)->min, vhd.at(ih)->max);
//       tightHistos.at(ih)->Sumw2();
// 
//       ratioHistos.at(ih) = new TH1F(vhd.at(ih)->histoName+"_ratio_"+outLabel+"_"+suffix, vhd.at(ih)->histoName+"_ratio_"+outLabel+"_"+suffix,
// 				    vhd.at(ih)->nBins, vhd.at(ih)->min, vhd.at(ih)->max);
//       ratioHistos.at(ih)->Sumw2();
//     }
//   }
// 
//   void fillLooseHistos(TprimeEvent* teve){
//     looseHistos.at(Lep1Pt_FR) ->Fill(teve->lepton1->lv.Pt(), teve->weight);
//     looseHistos.at(nJets_FR)  ->Fill(teve->vGoodJets.size(), teve->weight);
//     looseHistos.at(nEvents_FR)->Fill(1, teve->weight);
//     looseHistos.at(nPV_FR)    ->Fill(teve->nPV, teve->weight);
//     looseHistos.at(MET_FR)    ->Fill(teve->met, teve->weight);
//     looseHistos.at(MT_FR)     ->Fill(teve->mt, teve->weight);
// 
//     if (teve->vGoodJets.size() > 0){
//       looseHistos.at(Jet1Pt_FR)->Fill(teve->vGoodJets.at(0)->lv.Pt(), teve->weight);
//       for (unsigned int ij = 0; ij < teve->vGoodJets.size(); ij++){
// 	looseHistos.at(mJL_FR)->Fill((teve->lepton1->lv + teve->vGoodJets.at(ij)->lv).M(), teve->weight);
// 	looseHistos.at(dRJL_FR)->Fill(teve->lepton1->lv.DeltaR(teve->vGoodJets.at(ij)->lv), teve->weight);
//       }
//     }
//   }
// 
//   void fillTightHistos(TprimeEvent* teve){
//     tightHistos.at(Lep1Pt_FR) ->Fill(teve->lepton1->lv.Pt(), teve->weight);
//     tightHistos.at(nJets_FR)  ->Fill(teve->vGoodJets.size(), teve->weight);
//     tightHistos.at(nEvents_FR)->Fill(1, teve->weight);
//     tightHistos.at(nPV_FR)    ->Fill(teve->nPV, teve->weight);
//     tightHistos.at(MET_FR)    ->Fill(teve->met, teve->weight);
//     tightHistos.at(MT_FR)     ->Fill(teve->mt, teve->weight);
// 
//     if (teve->vGoodJets.size() > 0){
//       tightHistos.at(Jet1Pt_FR)->Fill(teve->vGoodJets.at(0)->lv.Pt(), teve->weight);
//       for (unsigned int ij = 0; ij < teve->vGoodJets.size(); ij++){
// 	tightHistos.at(mJL_FR)->Fill((teve->lepton1->lv + teve->vGoodJets.at(ij)->lv).M(), teve->weight);
// 	tightHistos.at(dRJL_FR)->Fill(teve->lepton1->lv.DeltaR(teve->vGoodJets.at(ij)->lv), teve->weight);
//       }
//     }
//   }
// 
//   void scaleTLHistos(){
//     for (unsigned int ih = 0; ih < NHISTOS_FR; ih++){
//       looseHistos.at(ih)->Scale(scale_);
//       tightHistos.at(ih)->Scale(scale_);
//     }
//   }
// 
//   void getRatioHistos(){
//     for (unsigned int ih = 0; ih < NHISTOS_FR; ih++){
//       ratioHistos.at(ih)->Add(tightHistos.at(ih));
//       ratioHistos.at(ih)->Divide(looseHistos.at(ih));
//     }
//   }
// 
//   void fillTLHistos(TprimeEvent* teve, bool b_Tight){
//     fillLooseHistos(teve);
//     if (b_Tight) fillTightHistos(teve);
//   }
// };



#endif
