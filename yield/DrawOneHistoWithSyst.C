#ifndef DrawOneHistoWithSyst_cxx
#define DrawOneHistoWithSyst_cxx

#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
//include "LoadData_2012.C"
#include "tdrstyle.C"

// #include "tPrimeStandardHistos.C"

using namespace std;
#include "CMS_lumi.C"

enum Samples_t {Data,    ChargeMisID, FakeRate,
		TTbarHad,   TTbarSL, TTbarLept, Top172, ZJets,   DY1050,WJets, SingleTop, DrellYan,
		WW,      WZ,      ZZ,      
		WGEl,    WGMu,
		T_tW,    Tbar_tW, T_t,     Tbar_t, T_s,   Tbar_s,
		WWm,     WWp,     WWW,  
		TTW,     TTZ, 	  TTWW,
		ZZZNoGs, WWZNoGs, WZZNoGs, WWds,
		WGstarToLNu2E, WGstarToLNu2Mu, WGstarToLNu2Tau, ZGToLLG,
		Top166, Top169, Top171, Top173, Top175, Top178,
 		last
};

const unsigned int totalSamples = last;

Samples_t allSamples[totalSamples] = {Data,    ChargeMisID, FakeRate,
		TTbarHad, TTbarSL, TTbarLept, Top172, ZJets,   DY1050,WJets, SingleTop, DrellYan,
		WW,      WZ,      ZZ,      
		WGEl,    WGMu,
		T_tW,    Tbar_tW, T_t,     Tbar_t, T_s,   Tbar_s,
		WWm,     WWp,     WWW,  
		TTW,     TTZ, 	  TTWW, 
		ZZZNoGs, WWZNoGs, WZZNoGs, WWds,
		WGstarToLNu2E, WGstarToLNu2Mu, WGstarToLNu2Tau, ZGToLLG,
		Top166, Top169, Top171, Top173, Top175, Top178
};


TString allNames[totalSamples] = {"Data",    "ChargeMisID", "FakeRate",
			    "TTbarHad",   "TTbarSL", "TTbarLept", "Top172", "ZJets",   "DY1050","WJets", "tW", "DrellYan",
			    "WW",      "WZ",      "ZZ", 
			    "WGEl",    "WGMu",
			    "T_tW",    "Tbar_tW", "T_t",     "Tbar_t", "T_s", "Tbar_s",
			    "WWm",     "WWp",     "WWW",  
			    "TTW",     "TTZ", "TTWW",
			    "ZZZNoGs", "WWZNoGs", "WZZNoGs", "WWds",
			    "WGstarToLNu2E", "WGstarToLNu2Mu", "WGstarToLNu2Tau", "ZGToLLG",
			    "Top166", "Top169", "Top171", "Top173", "Top175", "Top178"
};

TString legendNames[totalSamples] = {"Data",    "Charge MisID", "Non-prompt",
			    "TTbarHad",   "TTbarSL","t#bar{t}", "t#bar{t}",   "ZJets",   "DY1050","W+Jets", "tW", "Drell Yan",
			    "WW",      "WZ",      "ZZ",
			    "T_tW",    "Tbar_tW", "T_t",     "Tbar_t", "T_s", "Tbar_s",
			    "WWm",     "WWp",     "WWW",
			    "t#bar{t}W",     "t#bar{t}Z", 	  "t#bar{t}WW",
			     "ZZZNoGs", "WWZNoGs", "WZZNoGs", "Tri-bosons",
			     "WGstarToLNu2E", "WGstarToLNu2Mu", "WGstarToLNu2Tau", "ZGToLLG",
			    "Top166", "Top169", "Top171", "Top173", "Top175", "Top178"
};

Color_t color[totalSamples] = {
	kBlack, kCyan, kGray,				// Data,    ChargeMisID, FakeRate,
	kRed+1,kRed+1,kRed+1,kRed+1,   kAzure-2,kAzure-2,  kGreen-3, 	kMagenta, kAzure-2, // TTbar,   ZJets,   DY1050,Single Top, DY,WJets, SingleTop, DrellYan
	kOrange,kOrange,kOrange,  //WW,      WZ,      ZZ,
	//kGreen-9,      kViolet+1,      kAzure-2,
	kMagenta,    kMagenta, kMagenta,  kMagenta, kMagenta,   kMagenta, //Single Top
	kGreen-3,     kGreen-3,     kGreen+2,  	// WWm,     WWp,     WWW,
	kOrange,     kOrange,kOrange+2,		// TTW,     TTZ,     TTWW,
	};

vector <unsigned int> indices(totalSamples, -1);

//This is for OS:


#define REDUCED
// const unsigned int NSAMPLES = 4;
// Samples_t samples[NSAMPLES] = {
// Data, TTbar, SingleTop, DrellYan,
// //  T_tW,    Tbar_tW, T_t,     Tbar_t, T_s,   Tbar_s,
// //ZJets, DY1050
// };

const unsigned int NSAMPLES = 10;
Samples_t samples[NSAMPLES] = {
Data,
// TTbarHad, TTbarSL, TTbarLept, 
Top172,
ZJets, DY1050, T_tW,    Tbar_tW, WJets,
WW,      WZ,      ZZ
};


// float dySF[4] = {1.13825,1.,1.13171,1.}; was for pass3_test2;
float dySF[4] = {1.2524,1.0,1.32265,1.};


vector <TString> sNames;

vector <double> xSection;
vector <TH1D*>  Histos;
vector <double> sScale;
vector <double> systUnc;


/**
#DrawOneHisto parameter list:
#1: Input directory in root file
#2: Channel (ElEl, ElMu or MuMu
#3: Prefix for histograms
#4: Name of histogram
#5: Log or not log
#6: x-Axis name
#7: Optional: rebin by this number
*/

float getMyEntries(TH1* histo)
{
  float tot=0.;
  for (int i=0;i<histo->GetNbinsX()+1;++i) {
    if (histo->GetBinContent(i)>0.) tot+=histo->GetBinContent(i);
  }
  return tot;
}

float getMinEntries(TH1* histo)
{
  float min=999999.;
  for (int i=0;i<histo->GetNbinsX()+1;++i) {
    if ((histo->GetBinContent(i)>0.) && (histo->GetBinContent(i)<min)) min=histo->GetBinContent(i);
  }
  return min;
}
float getMinEntries(TH1* histo, float min)
{
  for (int i=0;i<histo->GetNbinsX()+1;++i) {
    if ((histo->GetBinContent(i)>0.) && (histo->GetBinContent(i)<min)) min=histo->GetBinContent(i);
  }
  return min;
}

double uncertainty (Samples_t sampleIndex );

int main( int argc, const char* argv[] ){
  if (argc < 8){
    cout<<"Need at least 7 arguments. Only "<<argc<<" found."<<endl;
    return 1;
  }

  setTDRStyle();

  //Root file information
  TString rootFile  = argv[1];
  TString inDir  = argv[2];
  TString inSuff = argv[3];
  TString inPref = argv[4];
  TString inName = argv[5];

  //Drawing information
  int yLog = atoi(argv[6]);
  TString xTitle = argv[7];

  int rebin = 1;
  if (argc > 8) rebin = atoi(argv[8]);
  float minX = -999, maxX;
  if (argc > 9) {
    minX = atof(argv[9]);
    maxX = atof(argv[10]);
  }

  cout<<inDir<<" ";
  cout<<inSuff<<" ";
  cout<<"inPref " <<inPref<<" ";
  cout<<"inName "<<inName<<" ";
  cout<<yLog<<" ";
  cout<<xTitle<<" ";
  cout<<rebin<<" ";
  cout<<minX<<" ";
  cout<<maxX<<endl;

  const unsigned int NCHAN = 4;

  TFile *file0 = TFile::Open(rootFile);
  if (file0==0) {
    cout << "\n\nERROR: "<<rootFile<< " does not exist\n\n";
    return 1;
  }
  file0->cd();

  TH1F* histArray[NSAMPLES][NCHAN];
//   if (inPref != "Top") 
  inPref = "_"+inPref;
//   else                 inPref = "";
  if (inDir != "Top") file0->cd(inDir);
  for (unsigned int isam = 0; isam < NSAMPLES; isam++){
//   cout << "ab "<<samples[isam]<< indices.size()<< endl;
    indices[samples[isam]] = isam;
    sNames.push_back(allNames[samples[isam]]);

    TString sampleName;
    if (isam < NSAMPLES) sampleName = sNames[isam];
    else                 sampleName = "FakeRate";

    cout << isam<<" "<<allNames[samples[isam]]<<" "<<color[samples[isam]]<<" "<<inName<<inPref<<"_"<<sampleName<<endl;;

    histArray[isam][0] = (TH1F*)gDirectory->Get(inName+inPref+"_"+sampleName+"_"+"MuMu");
    histArray[isam][1] = (TH1F*)gDirectory->Get(inName+inPref+"_"+sampleName+"_"+"ElMu");
    histArray[isam][2] = (TH1F*)gDirectory->Get(inName+inPref+"_"+sampleName+"_"+"ElEl");
    if ((histArray[isam][0]==0) || (histArray[isam][0]==0) || (histArray[isam][0]==0)) {
      cout << "\n\nERROR: One of the histograms do not exist: "<< inName+inPref+"_"+sampleName+"_"+"ElEl" <<" "<<inName<<inPref<<"_"<<sampleName<<"\n\n";
      gDirectory->ls();
      return 1;
    }
// cout <<  histArray[isam][0]->GetName()<<endl;
// cout <<  histArray[isam][0]->GetNbinsX()<<endl;
    if (rebin != 1){
      histArray[isam][0]->Rebin(rebin);
      histArray[isam][1]->Rebin(rebin);
      histArray[isam][2]->Rebin(rebin);
    }
    if (minX >-990){
      histArray[isam][0]->SetAxisRange(minX,maxX);
      histArray[isam][1]->SetAxisRange(minX,maxX);
      histArray[isam][2]->SetAxisRange(minX,maxX);
    }
    if ((samples[isam]== DY1050) || (samples[isam]== ZJets)) {
//       double initial =  histArray[isam][0]->Integral()+histArray[isam][1]->Integral()+histArray[isam][2]->Integral();
      histArray[isam][0]->Scale(dySF[0]);
      histArray[isam][1]->Scale(dySF[1]);
      histArray[isam][2]->Scale(dySF[2]);
//       double final =  histArray[isam][0]->Integral()+histArray[isam][1]->Integral()+histArray[isam][2]->Integral();
//       dySF[3] = final/inital;
    }
//     if (samples[isam]== Top172) {
//       histArray[isam][0]->Scale(1.020528046);
//       histArray[isam][1]->Scale(1.013066315);
//       histArray[isam][2]->Scale(1.02435166 );
//     }

    if ( (samples[isam]!=Data) && (samples[isam]!=FakeRate) && (samples[isam]!=ChargeMisID)){
//       histArray[isam][0]->Scale(19.62/16.52);
//       histArray[isam][1]->Scale(19.62/16.52);
//       histArray[isam][2]->Scale(19.62/16.52);
    }

    histArray[isam][3] = (TH1F*)histArray[isam][0]->Clone("histArray_"+sampleName);
    histArray[isam][3]->Add(histArray[isam][1]);
    histArray[isam][3]->Add(histArray[isam][2]);
  }

  for (int iChan = 0;iChan<4;++iChan) {
    switch( iChan ) {
      case 0: inSuff = "MuMu"; break;
      case 1: inSuff = "ElMu"; break;
      case 2: inSuff = "ElEl"; break;
      case 3: inSuff = "All" ; break;
    }
    vector <TH1F*> histos;
    for (unsigned int isam = 0; isam < NSAMPLES; isam++) {
      histos.push_back(histArray[isam][iChan]);
    }

    //Format data histogram
    histos.at(indices[Data])->SetMarkerStyle(20);
    for (int ibin = 0; ibin < histos.at(indices[Data])->GetNbinsX() + 1; ibin++){
      if (histos.at(indices[Data])->GetBinContent(ibin) == 0) histos.at(indices[Data])->SetBinContent(ibin, -10);
    }

    THStack* mystack = new THStack("mystack","mystack");
    TH1F * histo1D_mc = (TH1F *) histos.at(Data)->Clone("MC");
    histo1D_mc->Reset();
    TH1F * histo1D_multiBoson = (TH1F *) histos.at(Data)->Clone("MB");
    histo1D_multiBoson->Reset();
    TH1F * histo1D_ttBoson = (TH1F *) histos.at(Data)->Clone("ttv");
    histo1D_ttBoson->Reset();
    TH1F * histo1D_signal = (TH1F *) histos.at(Data)->Clone("signal");
    histo1D_signal->Reset();
    histo1D_signal->SetMarkerStyle(22);
    histo1D_signal->SetMarkerSize(1.5);
    histo1D_signal->SetMarkerColor(kBlue);
    TH1F * ratio = (TH1F *) histos.at(Data)->Clone("ratio");

// Legend:

  TLegend * leg = new TLegend(0.7, 0.8, 0.88, 0.94, "","brNDC");
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetLineStyle(0);

  leg->AddEntry(histos.at(Data), "Data", "lp");
//   leg->AddEntry(histo1D_signal, "T (800 GeV)", "lp");

  cout << "Yield : "<< inSuff<<endl;
  cout << "Yield: Data  "<<getMyEntries(histos.at(Data))<<endl;

  float singeTopYield = 0.;
  float ttYield = 0.;
  float dyYield = 0.;
  float totalBackground = 0.;
  float minimum=999999;
  double signalMax = 0.;
  bool mb=false;
  for (unsigned int isam = 0; isam < NSAMPLES; isam++){
    if (samples[isam]!=Data) {
      histos.at(isam)->SetFillColor(color[samples[isam]]);
cout << samples[isam]<< endl;
        minimum=getMinEntries(histos.at(isam),minimum);
        totalBackground+=histos.at(isam)->Integral();
	if (!((samples[isam]>=TTbarHad)&& (samples[isam]<=TTbarLept))) {
	  if (!((samples[isam]>=T_tW)&& (samples[isam]<=Tbar_s))) {
	    if ((samples[isam]!= DY1050) && (samples[isam]!= ZJets))
              cout << "Yield: "<< std::fixed<<std::setprecision(1)<<legendNames[samples[isam]] << " "<< histos.at(isam)->Integral()<<" "<<
	      endl;
	      else {dyYield+=histos.at(isam)->Integral();
	      }
	  } else singeTopYield+=histos.at(isam)->Integral();
	  } else ttYield+=histos.at(isam)->Integral();


	if (samples[isam]==T_tW)
	  leg->AddEntry(histos.at(isam), "Single Top", "f");
	if (samples[isam]==ZJets)
	  leg->AddEntry(histos.at(isam), "Drell Yann", "f");
	if (samples[isam]==TTbarHad)
	  leg->AddEntry(histos.at(isam), "t#bar{t}", "f");
// 	if (samples[isam]==ZZ)
// 	  leg->AddEntry(histos.at(isam), "WW", "f");

	if ((samples[isam]<=WJets) || 
	  ((samples[isam]>=T_tW)&& (samples[isam]<=Tbar_s)))
	  {
            cout << "Adding: "<< std::fixed<<std::setprecision(1)<<legendNames[samples[isam]] << " "<< histos.at(isam)->Integral()<<" "<<
	    histos.at(isam)->GetBinContent(1)<<
	    endl;
	  	mystack->Add(histos.at(isam));
		 if (!((samples[isam]>=TTbarHad)&& (samples[isam]<=TTbarLept)) &&
		 !((samples[isam]>=T_tW)&& (samples[isam]<=Tbar_s)) &&
		 (samples[isam]!= DY1050) && (samples[isam]!= ZJets) &&
		 (samples[isam]!= WWp) && (samples[isam]!= WWm))
		   leg->AddEntry(histos.at(isam), legendNames[samples[isam]], "f");
		}
	else if ((samples[isam]==WW) || (samples[isam]==WZ) || (samples[isam]==ZZ) ||
	(samples[isam]==WWm) || (samples[isam]==WWp) || (samples[isam]==WWW) ||
	(samples[isam]>=ZZZNoGs)){
	  if (!mb) leg->AddEntry(histo1D_multiBoson, "Multi-bosons", "f");
	  mb=true;
	  histo1D_multiBoson->Add(histos.at(isam));
	} else if ((samples[isam]==TTW) || (samples[isam]==TTZ) || (samples[isam]==TTWW)){
	if (histo1D_ttBoson->Integral()==0.) leg->AddEntry(histo1D_ttBoson, "t#bar{t}+bosons", "f");
		histo1D_ttBoson->Add(histos.at(isam));
    
        } else {
	cout << " Unknown : " << legendNames[samples[isam]]<<endl;
	}
	histo1D_mc->Add(histos.at(isam));

    }
  }

  cout << histo1D_multiBoson->Integral() <<" "<< histo1D_multiBoson->GetBinContent(1)<<
  endl;
  if (histo1D_multiBoson->Integral()>0.) {
    histo1D_multiBoson->SetFillColor(color[WW]);
    mystack->Add(histo1D_multiBoson);
  }
  cout << "ttV int "<<histo1D_ttBoson->Integral()<<" "<< histo1D_ttBoson->GetBinContent(1)<<endl;
  if (histo1D_ttBoson->Integral()>0.)    {
    histo1D_ttBoson->SetFillColor(color[TTW]);
    mystack->Add(histo1D_ttBoson);
  }


  cout << "Yield: tt: "<< std::fixed<<std::setprecision(1)<<ttYield<<endl;
  cout << "Yield: DY: "<< std::fixed<<std::setprecision(1)<<dyYield<<endl;
  cout << "Yield: Single top: "<< singeTopYield<<endl;
  cout << "Yield: total background: "<< std::fixed<<std::setprecision(1)<<totalBackground<<endl;
  cout << "Yield: Signal Tprime 500: "<< histo1D_signal->Integral()/500.<<endl;

    //Calculate systematics
    vector <TH1F*> selectedHistos;
    vector <double> selectedUnc;

    for (unsigned int isam = 0; isam < NSAMPLES; isam++){
      if (samples[isam]!=Data  ) {
	selectedHistos.push_back(histos.at(isam));
	selectedUnc.push_back(uncertainty(samples[isam]));
	//cout <<  "Unc: "<< allNames[samples[isam]]<<" "<<uncertainty(samples[isam])<<endl;
      }
    }

    TH1F* h_lumiBand = (TH1F*) histo1D_mc->Clone("h_lumiBand");

    for (int ibin = 1; ibin < h_lumiBand->GetNbinsX()+1; ibin++){

      double uncStat = 0;
      double uncSyst = 0;
      double uncTot  = 0;

      for (unsigned int ih = 0; ih < selectedHistos.size(); ih++){
	uncStat += selectedHistos.at(ih)->GetBinError(ibin)*selectedHistos.at(ih)->GetBinError(ibin);
	uncSyst += selectedHistos.at(ih)->GetBinContent(ibin)*selectedHistos.at(ih)->GetBinContent(ibin)*selectedUnc.at(ih)*selectedUnc.at(ih);
        //cout << allNames[samples[ih+1]]<< " "<<selectedHistos.at(ih)->GetBinContent(ibin) << " "<< selectedHistos.at(ih)->GetBinError(ibin) << " "<< selectedUnc.at(ih)<< endl;
      }

      uncTot = sqrt(uncStat + uncSyst);
      h_lumiBand->SetBinError(ibin, uncTot);
      //cout << "Bin "<< ibin << " "<< h_lumiBand->GetBinContent(ibin) << " "<< sqrt(uncStat) << " "<<sqrt(uncSyst)<< " "<<uncTot<<endl;
    }

  if (mystack->GetMaximum() > histos.at(indices[Data])->GetMaximum() + sqrt(histos.at(indices[Data])->GetMaximum())){
    if (yLog == 0) histos.at(indices[Data])->SetMaximum((mystack->GetMaximum() + sqrt(histos.at(indices[Data])->GetMaximum())) * 1.05);
    else           histos.at(indices[Data])->SetMaximum((mystack->GetMaximum() + sqrt(histos.at(indices[Data])->GetMaximum())) * 1.40);
  }


  //Larger description of channel
  TString schan = "";
  if (inSuff == "ElEl") schan = " - ee";
  if (inSuff == "ElMu") schan = " - e#mu";
  if (inSuff == "MuMu") schan = " - #mu#mu";
  if (inSuff == "All")  schan = "";

  if (yLog == 0) histos.at(indices[Data])->SetMinimum(0);
  else           histos.at(indices[Data])->SetMinimum(getMinEntries(histo1D_mc)/5);
//histos.at(indices[Data])->SetMinimum(10);
  double maximum = max(histos.at(indices[Data])->GetMaximum(),signalMax);
  maximum = max(maximum,histo1D_mc->GetMaximum());
  cout << "max: "<<maximum<<endl;
  if (yLog == 0) histos.at(indices[Data])->SetMaximum(maximum * 1.2);
  else           histos.at(indices[Data])->SetMaximum(maximum * 5);

//  histos.at(indices[Data])->GetXaxis()->SetTitle(xTitle+schan);
  histos.at(indices[Data])->GetYaxis()->SetTitle("Events");
  histos.at(indices[Data])->GetYaxis()->SetTitleOffset(1.3);

  TCanvas* canv = new TCanvas("canv", "canv", 600, 800);
  canv->SetLeftMargin(0.1);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.05);
  canv->SetBottomMargin(0.3);
  canv->cd();

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.7 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  int iPeriod = 2;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV 

  histos.at(indices[Data])->SetMarkerSize(1.2);
  histos.at(indices[Data])->GetXaxis()->SetLabelSize(0);
  histos.at(indices[Data])->GetXaxis()->SetTitleSize(0);
  histos.at(indices[Data])->Draw("P E x0");
  mystack->Draw("hist same");
//   mystack->GetHistogram()->Draw("same");
  histos.at(indices[Data])->Draw("P E same");

#ifndef REDUCED
//     for (unsigned int isam = 0; isam < NSAMPLES; isam++){
//       if (samples[isam]>=Tprime400_bWbW) {
// // 	histos.at(isam)->Draw("P E same x0");
//  	histos.at(isam)->Draw("hist same");
//       }
//     }

//   histo1D_signal->Draw("P E same x0");
#endif

  gStyle->SetHatchesLineWidth(4);
  gStyle->SetErrorX(0.5);

  h_lumiBand->SetFillStyle(3005);
  h_lumiBand->SetFillColor(1);
  h_lumiBand->SetMarkerStyle(1);
  h_lumiBand->Draw("same e2");

  leg->Draw();

  //Mandatory text boxes
  CMS_lumi( canv, iPeriod, 0 );

  if (getMyEntries(histos.at(Data))>0. && (yLog == 1) ) gPad->SetLogy(yLog);

    TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
    pad->SetTopMargin(0.7);
    pad->SetLeftMargin(0.1);
    pad->SetRightMargin(0.03);
    pad->SetFillColor(0);
    pad->SetFillStyle(0);
//     pad->SetTickx(1);
//     pad->SetTicky(1);
    pad->Draw();
    pad->cd(0);
    pad->SetGridy(true);
    ratio->Divide(h_lumiBand);
//    ratio->Divide(histo1D_mc);
    ratio->GetXaxis()->SetTitle(xTitle+schan);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->Draw("p e x0");
    ratio->SetMaximum( 1.5);
    ratio->SetMinimum( 0.5);


//    histos.at(Data)->Draw("same");


  TString outName = inName+inPref;
  if (inSuff != "All") outName += "_"+inSuff;


  canv->SaveAs(outName+".png");
  canv->SaveAs(outName+".pdf");
//  canv->SaveAs(outName+".C");
  }
}
// int main( int argc, const char* argv[] ){
//
// TString channel[4] = {"ElEl","ElMu","MuMu", "All"}
//
//   TString rootFile  = argv[1];
//   TString inDir  = argv[2];
//   TString inSuff = argv[3];
//   TString inPref = argv[4];
//   TString inName = argv[5];
//
//   //Drawing information
//   int yLog = atoi(argv[6]);
//   TString xTitle = argv[7];
//
//   vector<TString> legend(NHISTOS);
//   legend.at(Lep1Pt)    = TString("Leading lepton p_{T} [GeV/c]");
// legend.at(Lep2Pt)    = TString("Second lepton p_{T} [GeV/c]");
// legend.at(Lep3Pt)    = TString("Third lepton p_{T} [GeV/c]");
//
// legend.at(ElPt)	  = TString("Electron p_{T} [GeV/c]");
// legend.at(MuPt)	  = TString("Muon p_{T} [GeV/c]");
//
// legend.at(ElEta)	  = TString("");
// legend.at(MuEta)	  = TString("");
//
// legend.at(Jet1Pt)    = TString("Leading jet p_{T} [GeV/c]");
// legend.at(Jet2Pt)    = TString("Second jet p_{T} [GeV/c]");
//
// legend.at(nJets)	  = TString("Number of Jets"	    );
// legend.at(nBJets)    = TString("Number of b-tagged jets");
// legend.at(nElec)	  = TString("Number of Electrons");
// legend.at(nMuon)	  = TString("Number of Muons");
// legend.at(nLept)	  = TString("Number of leptons");
// legend.at(sumPtL)    = TString("Sum of lepton p_{T} [GeV/c]");
// legend.at(sumPtJ)    = TString("Sum of jet p_{T} [GeV/c]");
// legend.at(MET)	  = TString("MET [GeV]");
//
// legend.at(HT)	  = TString("HT [GeV]");
//
// legend.at(LepInvM)   = TString("Invariant Mass of Lepton Pair [GeV/c^{2}]");
// legend.at(LepJInvM)  = TString("");
// legend.at(dRlept)      = TString("#delta R (l1,l2)");
// legend.at(mindR)	    = TString("Min (#delta R (l,b-jet))");
// legend.at(dR)	    = TString("#delta R (l,b-jet)");
// legend.at(minMlb)      = TString("Min(m_{lb}})");
// legend.at(Mlb)	    = TString("m_{lb}}");
//
// };
//
//
//   for (int=0,i<4;++i)
//   {
//     for (unsigned int ih = 0; ih < NHISTOS; ih++){
//       drawSinglePlot(rootFile, inDir, channel[i] , "Lep1Pt"  1 legend[ih], 10);
//
// }
//
// }

double uncertainty (Samples_t sampleIndex )
{
  double mcUnc = sqrt (2 * 0.006 * 0.006 + //Trigger
		       2 * 0.03  * 0.03);  //Lepton Efficiency
  double systUnc;

      switch( sampleIndex )
      {
        case Data:	   systUnc = 0.  ; break;
	case ChargeMisID:  systUnc = 0.20; break;
	case FakeRate:	   systUnc = 0.50; break;

	case TTbarHad:	   systUnc = 0.12; break;
	case TTbarSL:	   systUnc = 0.12; break;
	case TTbarLept:	   systUnc = 0.12; break;
	case ZJets:	   systUnc = 0.30; break;
	case DY1050:	   systUnc = 0.30; break;
	case WJets:	   systUnc = 0.30; break;

	case T_tW:	   systUnc = 0.20;  break;
	case Tbar_tW:	   systUnc = 0.20;  break;
	case T_t:	   systUnc = 0.20;  break;
	case Tbar_t:	   systUnc = 0.20;  break;
	case T_s:	   systUnc = 0.20;  break;
	case Tbar_s:	   systUnc = 0.20;  break;

	case WZ:	   systUnc = 0.17; break;
	case ZZ:	   systUnc = 0.07; break;

	case WWm:	   systUnc = 0.50; break;
	case WWp:	   systUnc = 0.50; break;
	case WWW:	   systUnc = 0.50; break;

	case TTW:	   systUnc = 0.32; break;
	case TTWW:	   systUnc = 0.50; break;
	case TTZ:	   systUnc = 0.50; break;

	default:systUnc = mcUnc;
      }
  return systUnc;
}
#endif
