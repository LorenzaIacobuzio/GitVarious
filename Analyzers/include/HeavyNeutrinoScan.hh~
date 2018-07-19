// ---------------------------------------------------------------                 
// History:                                                                
//                                                                              
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018        
//                                                                                
// --------------------------------------------------------------- 

#ifndef HEAVYNEUTRINOSCAN_HH
#define HEAVYNEUTRINOSCAN_HH

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include "Analyzer.hh"

class TH1I;
class TH2F;
class TGraph;
class TTree;
class TGraphAsymmErrors;

class HeavyNeutrinoScan : public NA62Analysis::Analyzer {

public:

  HeavyNeutrinoScan(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrinoScan() {}
  void InitHist();
  void InitOutput() {}
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
  
protected:
  
  // Scan variables                                                                 

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
  Int_t fMode;
  Int_t fN;
  Int_t fNMom;

  std::map<Double_t, std::map<Double_t, Double_t>> fGammaTot;
  std::map<Double_t, std::map<Double_t, Double_t>> fTau;
  std::map<Double_t, std::map<Double_t, Int_t>>    fNevents;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGood;
  std::map<Double_t, std::map<Double_t, Double_t>> fYield;
  std::map<Double_t, Double_t>                     fCouplings;
  std::map<Double_t, Double_t>                     fMasses;
};

#endif
