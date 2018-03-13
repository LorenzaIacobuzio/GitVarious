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
  void PlotErrorBars(TGraphErrors*, TGraphErrors*, TGraphErrors*, TGraphErrors*, std::string);
  void PlotErrorBarsMom(TGraphErrors*, std::string);
  void CosmeticsGraph(TGraphErrors*, const char*, const char*, Int_t);
  std::vector<Double_t> ComputeRMS(std::vector<Double_t>);
  std::vector<std::map<std::string, std::vector<Double_t>>> SplitVector(std::vector<Double_t>, std::string);
  
protected:
  
  // Scan variables                                                                                   

  Bool_t   fReadingData;
  Bool_t   fEnableCouplingScan;
  fstream  fErrorFile;
  fstream  fErrorFileTarget;
  fstream  fErrorFileTAX;
  fstream  fErrorFileMom;
  Double_t fInitialFV;
  Double_t fLFV;
  Double_t fMomStart;
  Double_t fMomStop;
  Double_t fMomStep;
  Double_t fMassForSingleValue;
  Double_t fCouplingStart;
  Double_t fCouplingStop;
  Double_t fCouplingStep;
  Int_t    fN;
  Int_t    fNMom;
  Int_t    fSplitStep;

  std::map<Double_t, std::map<Double_t, Int_t>>    fNevents;
  std::map<Double_t, std::map<Double_t, Int_t>>    fNeventsTarget;
  std::map<Double_t, std::map<Double_t, Int_t>>    fNeventsTAX;
  std::map<Double_t, Int_t>                        fNeventsMom;
  std::map<Double_t, Double_t>                     fCouplings;
  std::map<Double_t, Double_t>                     fMasses;
  std::map<Double_t, Double_t>                     fMomenta;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumAll;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumAllTarget;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumAllTAX;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGood;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGoodTarget;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGoodTAX;
  std::map<Double_t, Double_t>                     fSumGoodMom;
  std::map<Double_t, std::map<Double_t, Double_t>> fAcc;
  std::map<Double_t, std::map<Double_t, Double_t>> fAccTarget;
  std::map<Double_t, std::map<Double_t, Double_t>> fAccTAX;
  std::map<Double_t, std::map<Double_t, Double_t>> fYield;
  std::map<Double_t, std::map<Double_t, Double_t>> fYieldTarget;
  std::map<Double_t, std::map<Double_t, Double_t>> fYieldTAX;
  std::map<Double_t, Double_t>                     fYieldMom;
  std::map<Double_t, std::map<Double_t, Double_t>> fProb;
  std::map<Double_t, std::map<Double_t, Double_t>> fGammaTot;
  std::map<Double_t, std::map<Double_t, Double_t>> fTau;

  // One value histos

  TH1D *fhZMotherProd;
  TH1D *fhZDProd;
  TH1D *fhZTauProd;
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
  TH2D *fhProbCoupling;
  TH2D *fhWeightCoupling;
  TH2D *fhReachMass;
  TH2D *fhDecayMass;
  TH2D *fhProbMass;
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
  TGraphErrors *fgYieldMom;
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
  TGraphErrors *fgErrorYieldMom;
  TGraphErrors *fgErrorExclusion;  
};

#endif
