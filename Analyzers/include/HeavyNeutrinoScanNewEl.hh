// ---------------------------------------------------------------                 
// History:                                                                
//                                                                              
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018        
//                                                                                
// --------------------------------------------------------------- 

#ifndef HEAVYNEUTRINOSCANNEWEL_HH
#define HEAVYNEUTRINOSCANNEWEL_HH

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include "Analyzer.hh"
#include "TwoLinesCDA.hh"
#include "PointLineDistance.hh"

class TH1I;
class TH2F;
class TGraph;
class TTree;
class TGraphAsymmErrors;

class HeavyNeutrinoScanNewEl : public NA62Analysis::Analyzer {

public:

  HeavyNeutrinoScanNewEl(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrinoScanNewEl() {}
  void InitHist();
  void InitOutput();
  void DefineMCSimple() {}
  void Process(Int_t);
  void StartOfBurstUser() {}
  void EndOfBurstUser() {}
  void StartOfRunUser() {}
  void EndOfRunUser() {}
  void EndOfJobUser();
  void PostProcess() {}
  void DrawPlot() {}
  void CosmeticsGraph(TGraph*, const char*, const char*, Int_t);
  void CosmeticsGraph(TGraphAsymmErrors*, const char*, const char*, Int_t);
  void CosmeticsHisto(TH1*, const char*, const char*, Int_t);
  void EvaluateUL(TH2*, TGraph*);
  std::vector<TGraph*> ExtractContours(TH2*);
  void SumGraphs(TGraphAsymmErrors*, TGraphAsymmErrors*, TGraphAsymmErrors*);
  void Normalize(TH1D*, TH1D*, TString, TString, TString);


protected:
  
  TRandom3* r;
  TwoLinesCDA *fCDAcomp;
  PointLineDistance *fDistcomp;
  Bool_t fReadingData;
  Double_t fInitialFV;
  Double_t fLFV;
  Double_t fMomStart;
  Double_t fMomStop;
  Double_t fMomStep;
  Double_t fMassForSingleValue;
  Double_t fCouplingForSingleValue;
  Double_t fCouplingStart;
  Double_t fCouplingStop;
  Double_t fCouplingStep;
  Double_t fMassStart;
  Double_t fMassStop;
  Double_t fMassStep;
  Double_t fTMass;
  Double_t fTCoupling;
  Double_t fTNEvents;
  Double_t fTSumGood;
  Int_t fMode;
  Int_t fN;
  Int_t fNMom;
  Int_t fNMass;
  Int_t fNContours;

  std::map<Double_t, std::map<Double_t, Double_t>> fNEvents;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGood;
  std::map<Double_t, std::map<Double_t, Double_t>> fYield;
  std::map<Double_t, Double_t> fCouplings;
  std::map<Double_t, Double_t> fMasses;
};

#endif
