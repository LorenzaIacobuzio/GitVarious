#ifndef HEAVYNEUTRINOMASSSCAN_HH
#define HEAVYNEUTRINOMASSSCAN_HH

#include <stdlib.h>
#include <fstream>
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

class HeavyNeutrinoMassScan : public NA62Analysis::Analyzer {

public:

  HeavyNeutrinoMassScan(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrinoMassScan();
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

  // Other parameters                                                                                   

  TwoLinesCDA *fCDAcomp;
  PointLineDistance *fDistcomp;
  LAVMatching *fLAVMatching;
  SAVMatching *fSAVMatching;

  // Scan variables

  static constexpr Double_t fInitialMass = 100.;
  static constexpr Double_t fFinalMass = 2000.;
  static constexpr Double_t fMassStep = 100.;
  static constexpr Int_t fN = round((std::abs(fFinalMass-fInitialMass))/fMassStep);
  Int_t fNevents[fN];
  Double_t fSumAll[fN];
  Double_t fSumGood[fN];
  Double_t fAcc[fN];
  Double_t fGammaTot[fN];
  Double_t fTau[fN];
  Double_t fProb[fN];
  Double_t fYield[fN];
  Double_t fMasses[fN];

  // Histos

  TH2D *fhReach;
  TH2D *fhDecay;
  TH2D *fhWeight;

  TGraph *fgAcc;
  TGraph *fgYield;
  TGraph *fgGammaTot;
  TGraph *fgTau;
};

#endif
