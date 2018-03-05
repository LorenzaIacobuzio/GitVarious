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

class HeavyNeutrinoScan : public NA62Analysis::Analyzer {

public:

  HeavyNeutrinoScan(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrinoScan();
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
  std::vector<Double_t> ComputeRMS(std::vector<Double_t>);
  void PlotErrorBars(TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, std::string);

protected:

  // Scan variables                                                                                   

  Bool_t   fReadingData;
  fstream  fErrorFile;
  fstream  fErrorFileTarget;
  fstream  fErrorFileTAX;
  Double_t fInitialFV;
  Double_t fLFV;
  Double_t fMassForSingleValue;
  Double_t fCouplingStart;
  Double_t fCouplingStop;
  Double_t fCouplingStep;
  Bool_t   fEnableCouplingScan;
  Int_t    fN;
  Int_t    fErrorCounter;
  Int_t    fErrorCounterTarget;
  Int_t    fErrorCounterTAX;
  Int_t    fErrorStep;

  std::map<Double_t, std::map<Double_t, Int_t>>    fNevents;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumAll;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGood;
  std::map<Double_t, std::map<Double_t, Int_t>>    fNeventsTarget;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGoodTarget;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumAllTarget;
  std::map<Double_t, std::map<Double_t, Int_t>>    fNeventsTAX;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGoodTAX;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumAllTAX;
  std::map<Double_t, std::map<Double_t, Double_t>> fAcc;
  std::map<Double_t, std::map<Double_t, Double_t>> fAccTarget;
  std::map<Double_t, std::map<Double_t, Double_t>> fAccTAX;
  std::map<Double_t, std::map<Double_t, Double_t>> fYield;
  std::map<Double_t, std::map<Double_t, Double_t>> fYieldTarget;
  std::map<Double_t, std::map<Double_t, Double_t>> fYieldTAX;
  std::map<Double_t, std::map<Double_t, Double_t>> fProb;
  std::map<Double_t, std::map<Double_t, Double_t>> fGammaTot;
  std::map<Double_t, std::map<Double_t, Double_t>> fTau;
  std::map<Double_t, Double_t>                     fCouplings;
  std::map<Double_t, Double_t>                     fMasses;

  // One value histos

  TH1D *fhZDProd;
  TH1D *fhZDDecay;
  TH1D *fhDTheta;
  TH1D *fhDLambda;
  TH1D *fhDPath;
  TH1D *fhDMom;
  TH1D *fhZHNLDecay;
  TH1D *fhHNLGamma;
  TH1D *fhHNLDecayProb;
  TH1D *fhHNLReachProb;
  TH1D *fhHNLTheta;
  TH1D *fhHNLMom;
  TH1D *fhWeight;
  TH1D *fhCoupling;
  TH1D *fhMass;
  TH1D *fhAcc;
  TH1D *fhAccTarget;
  TH1D *fhAccTAX;
  TH1D *fhYield;
  TH1D *fhYieldTarget;
  TH1D *fhYieldTAX;

  // Scan histos                                                                  

  TH2D *fhReachCoupling;
  TH2D *fhDecayCoupling;
  TH2D *fhWeightCoupling;
  TH2D *fhReachMass;
  TH2D *fhDecayMass;
  TH2D *fhWeightMass;

  TGraphErrors *fgAccCoupling;
  TGraphErrors *fgAccCouplingTarget;
  TGraphErrors *fgAccCouplingTAX;
  TGraphErrors *fgYieldCoupling;
  TGraphErrors *fgYieldCouplingTarget;
  TGraphErrors *fgYieldCouplingTAX;
  TGraphErrors *fgGammaTotCoupling;
  TGraphErrors *fgTauCoupling;
  TGraphErrors *fgAccMass;
  TGraphErrors *fgAccMassTarget;
  TGraphErrors *fgAccMassTAX;
  TGraphErrors *fgYieldMass;
  TGraphErrors *fgYieldMassTarget;
  TGraphErrors *fgYieldMassTAX;
  TGraphErrors *fgGammaTotMass;
  TGraphErrors *fgTauMass;
  TGraphErrors *fgExclusion;

  // Error bar histos

  TGraphErrors *fgErrorAccCoupling;
  TGraphErrors *fgErrorAccCouplingTarget;
  TGraphErrors *fgErrorAccCouplingTAX;
  TGraphErrors *fgErrorYieldCoupling;
  TGraphErrors *fgErrorYieldCouplingTarget;
  TGraphErrors *fgErrorYieldCouplingTAX;
  TGraphErrors *fgErrorAccMass;
  TGraphErrors *fgErrorAccMassTarget;
  TGraphErrors *fgErrorAccMassTAX;
  TGraphErrors *fgErrorYieldMass;
  TGraphErrors *fgErrorYieldMassTarget;
  TGraphErrors *fgErrorYieldMassTAX;
  TGraphErrors *fgErrorExclusion;  
};

#endif
