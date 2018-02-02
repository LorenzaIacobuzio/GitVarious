#ifndef HEAVYNEUTRINOMASSSCAN_HH
#define HEAVYNEUTRINOMASSSCAN_HH

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <iterator>
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

  TwoLinesCDA *fCDAcomp;
  PointLineDistance *fDistcomp;
  LAVMatching *fLAVMatching;
  SAVMatching *fSAVMatching;

  // Scan variables

  std::map<Double_t, Int_t>    fNevents;
  std::map<Double_t, Double_t> fSumAll;
  std::map<Double_t, Double_t> fSumGood;
  std::map<Double_t, Double_t> fAcc;
  std::map<Double_t, Double_t> fGammaTot;
  std::map<Double_t, Double_t> fTau;
  std::map<Double_t, Double_t> fProb;
  std::map<Double_t, Double_t> fYield;
  std::map<Double_t, Double_t> fMasses;

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
