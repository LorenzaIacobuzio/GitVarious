#ifndef HEAVYNEUTRINOTOTALSCAN_HH
#define HEAVYNEUTRINOTOTALSCAN_HH

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
#include "TH3.h"
#include "TGraph2D.h"

class TH1I;
class TH2F;
class TGraph;
class TTree;

class HeavyNeutrinoTotalScan : public NA62Analysis::Analyzer {

public:

  HeavyNeutrinoTotalScan(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrinoTotalScan();
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

  TwoLinesCDA *fCDAcomp;
  PointLineDistance *fDistcomp;
  LAVMatching *fLAVMatching;
  SAVMatching *fSAVMatching;

  // Scan variables

  static constexpr Double_t fCouplingStart = -10.;
  static constexpr Double_t fCouplingStop = 0.;
  static constexpr Double_t fCouplingStep = 0.1;
  static constexpr Int_t fN = round((std::abs(fCouplingStop-fCouplingStart))/fCouplingStep);
  std::map<Double_t, std::map<Double_t, Int_t>>    fNevents;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumAll;
  std::map<Double_t, std::map<Double_t, Double_t>> fSumGood;
  std::map<Double_t, std::map<Double_t, Double_t>> fAcc;
  std::map<Double_t, std::map<Double_t, Double_t>> fYield;
  std::map<Double_t, std::map<Double_t, Double_t>> fProb;
  std::map<Double_t, Double_t>                     fCouplings;
  std::map<Double_t, Double_t>                     fMasses;

  // Histos

  TGraph   *fgExclusion;
};

#endif
