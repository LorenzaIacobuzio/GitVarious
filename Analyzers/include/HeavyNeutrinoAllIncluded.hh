#ifndef HEAVYNEUTRINOALLINCLUDED_HH
#define HEAVYNEUTRINOALLINCLUDED_HH

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <TCanvas.h>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include "DetectorAcceptance.hh"
#include "SpectrometerTrackVertex.hh"
#include "TwoLinesCDA.hh"
#include "PointLineDistance.hh"
#include "LAVMatching.hh"
#include "SAVMatching.hh"
#include "MCInfo.hh"

class TH1I;
class TH2F;
class TGraph;
class TTree;

class HeavyNeutrinoAllIncluded : public NA62Analysis::Analyzer {

public:

  HeavyNeutrinoAllIncluded(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrinoAllIncluded();
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

protected:

  // Scan variables                                                                                   

  Double_t fCouplingStart;
  Double_t fCouplingStop;
  Double_t fCouplingStep;
  Bool_t fEnableCouplingScan;
  Int_t fN;
  std::map<Double_t, std::map<Double_t, Int_t>>    fNevents;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumAll;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGood;
  std::map<Double_t, std::map<Double_t, Int_t>>    fNeventsTarget;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumAllTarget;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGoodTarget;
  std::map<Double_t, std::map<Double_t, Int_t>>    fNeventsTAX;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumAllTAX;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGoodTAX;
  std::map<Double_t, std::map<Double_t, Double_t>> fAcc;
  std::map<Double_t, std::map<Double_t, Double_t>> fYield;
  std::map<Double_t, std::map<Double_t, Double_t>> fYieldTarget;
  std::map<Double_t, std::map<Double_t, Double_t>> fYieldTAX;
  std::map<Double_t, std::map<Double_t, Double_t>> fProb;
  std::map<Double_t, std::map<Double_t, Double_t>> fGammaTot;
  std::map<Double_t, std::map<Double_t, Double_t>> fTau;
  std::map<Double_t, Double_t>                     fCouplings;
  std::map<Double_t, Double_t>                     fMasses;

  // Scan histos                                                                  

  TH1D *fhAcc;
  TH1D *fhYield;
  TH1D *fhYieldTarget;
  TH1D *fhYieldTAX;
  TH2D *fhReachCoupling;
  TH2D *fhDecayCoupling;
  TH2D *fhWeightCoupling;
  TH2D *fhReachMass;
  TH2D *fhDecayMass;
  TH2D *fhWeightMass;

  TGraph *fgAccCoupling;
  TGraph *fgYieldCoupling;
  TGraph *fgGammaTotCoupling;
  TGraph *fgTauCoupling;
  TGraph *fgAccMass;
  TGraph *fgYieldMass;
  TGraph *fgGammaTotMass;
  TGraph *fgTauMass;
  TGraph *fgExclusion;
};

#endif
