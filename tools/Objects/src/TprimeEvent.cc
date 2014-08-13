#include "../interface/TprimeEvent.h"

double TprimeEvent::leptonMETScalSum()
{
  return (leptonScalSum()+ met);
}

double TprimeEvent::leptonJetsScalSum()
{
 return ( leptonScalSum() + ht());
}

double TprimeEvent::jetsMETScalSum()
{
 return ( ht() + met);
}

double TprimeEvent::leptonJetsMETScalSum()
{
 return ( leptonScalSum() + ht() + met);
}



TLorentzVector TprimeEvent::leptonMETVectSum()
{
  TLorentzVector sumlepMET = leptonSum();
  sumlepMET.SetPx(sumlepMET.Px()+met*cos(met_phi));
  sumlepMET.SetPy(sumlepMET.Py()+met*sin(met_phi));
  return sumlepMET;
}

TLorentzVector TprimeEvent::leptonJetsVectSum()
{
 TLorentzVector sumlepJets = leptonSum() + jetVectSum();
  return sumlepJets;

}

TLorentzVector TprimeEvent::jetsMETVectSum()
{
 TLorentzVector sumJetsMET = jetVectSum();
  sumJetsMET.SetPx(sumJetsMET.Px()+met*cos(met_phi));
  sumJetsMET.SetPy(sumJetsMET.Py()+met*sin(met_phi));
  return sumJetsMET;

}

TLorentzVector TprimeEvent::leptonJetsMETVectSum()
{
 TLorentzVector sumLeptonJetsMET = leptonSum() + jetVectSum();
  sumLeptonJetsMET.SetPx(sumLeptonJetsMET.Px()+met*cos(met_phi));
  sumLeptonJetsMET.SetPy(sumLeptonJetsMET.Py()+met*sin(met_phi));
  return sumLeptonJetsMET;
}


int TprimeEvent::nBTags()
{
  int n = 0;
  for (unsigned int ij = 0; ij < vGoodJets.size(); ij++){
    if (vGoodJets.at(ij)->csvMedium != 0) ++n;
  }
  return n;
}
int TprimeEvent::nAK5Tags()
{
  int n = 0;
  for (unsigned int ij = 0; ij < vGoodAK5Jets.size(); ij++){
    if (vGoodAK5Jets.at(ij)->csvMedium != 0) ++n;
  }
  return n;
}

TLorentzVector TprimeEvent::leptonSum()
{
  TLorentzVector sum;
  for (unsigned int j = 0;j<theGoodLeptons.size();++j )
    sum+= theGoodLeptons[j]->lv;
  return sum;
}

double TprimeEvent::leptonScalSum()
{
  double sum = 0.;
  for (unsigned int j = 0;j<theGoodLeptons.size();++j )
    sum+= theGoodLeptons[j]->lv.Pt();
  return sum;
}

TLorentzVector TprimeEvent::jetVectSum()
{
  TLorentzVector sum;
  for (unsigned int j = 0;j<vGoodJets.size();++j )
    sum+= vGoodJets[j]->lv;
  return sum;
}

double TprimeEvent::ht()
{
  double ht_=0.;
  for (unsigned int ij = 0; ij < vGoodJets.size(); ij++)
    ht_ += vGoodJets.at(ij)->lv.Pt();
  return ht_;
}

double TprimeEvent::dRlepton() const {
  return theGoodLeptons[0]->lv.DeltaR(theGoodLeptons[1]->lv);
}

vector<double> TprimeEvent::drBJetLepton () const {
  vector<double> dr;
  for (unsigned int i=0;i<vGoodJets.size();++i) {
    if (vGoodJets[i]->csvMedium != 0) {
      for (unsigned int j = 0;j<theGoodLeptons.size();++j ) {
	dr.push_back(theGoodLeptons[j]->lv.DeltaR(vGoodJets[i]->lv));
      }
    }
  }
  return dr;
}

double TprimeEvent::mindrBJetLepton() const {
     if (vGoodJets.size()==0) return -1;
     vector<double> dr(drBJetLepton());
     if (dr.size()==0) return -1;
     return *(min_element(dr.begin() ,dr.end()));

}

vector<double> TprimeEvent::drJetLepton (const vector <TJet*> & theJets) const {
  vector<double> dr;
  for (unsigned int i=0;i<theJets.size();++i) {
    for (unsigned int j = 0;j<theGoodLeptons.size();++j ) {
      dr.push_back(theGoodLeptons[j]->lv.DeltaR(theJets[i]->lv));
    }
  }
  return dr;
}

double TprimeEvent::mindrJetLepton(const vector <TJet*> & theJets) const {
   if (theJets.size()==0) return -1;
   vector<double> dr(drJetLepton(theJets));
   if (dr.size()==0) return -1;
   return *(min_element(dr.begin() ,dr.end()));
}


vector<double> TprimeEvent::mlb () const {
  vector<double> m;
  for (unsigned int i=0;i<vGoodJets.size();++i) {
    if (vGoodJets[i]->csvMedium != 0) {
      for (unsigned int j = 0;j<theGoodLeptons.size();++j )
        m.push_back((theGoodLeptons[j]->lv+vGoodJets[i]->lv).M());
    }
  }
  return m;
}
double TprimeEvent::minMlb() const {
     if (vGoodJets.size()==0) return -1;
     vector<double> m(mlb());
     if (m.size()==0) return -1;
     return *(min_element(m.begin() ,m.end()));
}
double TprimeEvent::maxMlb() const {
     if (vGoodJets.size()==0) return -1;
     vector<double> m(mlb());
     if (m.size()==0) return -1;
     return *(max_element(m.begin() ,m.end()));
}

void TprimeEvent::clearEvent(){
  weight = -1000.0;
  theGoodLeptons.clear();
  vGoodCATopJets.clear();
  vGoodCA8Jets.clear();
  vGoodAK5Jets.clear();
  vGoodJets.clear();
  vGoodCsvs.clear();
  met = 0;
  mt  = 0;
  amwtPeakMass = 0.0;
  amwtPeakWeight = 0.0;

}

void TprimeEvent::printEvent(){
  for (unsigned int ij = 0; ij < vGoodJets.size(); ij++){
    cout<<"Jet"<<ij<<endl;
    cout<<"Pt: "<<vGoodJets.at(ij)->lv.Pt()<<'\t'<<"Eta: "<<vGoodJets.at(ij)->lv.Eta()<<'\t'<<"Phi: "<<vGoodJets.at(ij)->lv.Phi()<<endl;
  }
  cout<<"HT: "<<ht()<<endl<<endl;
  cout<<"InvM: "<<leptonSum().M()<<endl;
}

vector<TJet*> TprimeEvent::selectJetsForMass()
{
  vector <TJet*> jets;
  if (vGoodJets.size()>2) return jets;
  unsigned int j[2];
  
  unsigned int nj=0;
  for (unsigned int ij = 0; ij < vGoodJets.size(); ij++){
    if ((vGoodJets.at(ij)->csvMedium != 0)&&(nj<2)) j[nj++]=ij;
  }
  if (nj<2) {
    for (unsigned int ij = 0; ij < vGoodJets.size(); ij++){
      if ((nj==1)&&(j[0]!=ij)) j[nj++]=ij;
      if (nj==0) j[nj++]=ij;
    }
  }
  jets.push_back(vGoodJets.at(j[0]));
  jets.push_back(vGoodJets.at(j[1]));
  return jets;
}
