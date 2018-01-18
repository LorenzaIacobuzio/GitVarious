#ifndef HEAVYNEUTRINOPIMUSELECTIONSCAN_HH
#define HEAVYNEUTRINOPIMUSELECTIONSCAN_HH

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

class HeavyNeutrinoPiMuSelectionScan : public NA62Analysis::Analyzer {

public:

  HeavyNeutrinoPiMuSelectionScan(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrinoPiMuSelectionScan();
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

  Double_t ComputeHNLMass(KinePart*);
  Double_t ComputeL(TVector3, TVector3, TVector3);
  Double_t ComputeNDecayProb(KinePart*, Double_t, Double_t);
  Double_t ComputeNReachProb(KinePart*, Double_t, Double_t);
  Double_t GammaPiE(Double_t);
  Double_t GammaPiMu(Double_t);
  Double_t GammaENuNu(Double_t);
  Double_t GammaTot(Double_t, Double_t, Double_t, Double_t);
  Double_t Tau(Double_t);
  Double_t NIntoPiMuBR(Double_t, Double_t, Double_t);
  Double_t PhaseSpace(Double_t, Double_t, Double_t);
  Double_t PhaseSpaceFactor(Double_t, Double_t, Double_t, Double_t);
  Double_t ComputeBR(KinePart*, Double_t);
  Double_t ComputeTotalBR(KinePart*, Double_t, Double_t); 

protected:

  TwoLinesCDA *fCDAcomp;
  PointLineDistance *fDistcomp;
  LAVMatching *fLAVMatching;
  SAVMatching *fSAVMatching;
  Int_t fIndex;
  Double_t fMN;
  static constexpr Double_t fCouplingStart = -10.;
  static constexpr Double_t fCouplingStop = 0.;
  static constexpr Double_t fCouplingStep = 0.1;
  static constexpr Int_t fN = round((std::abs(fCouplingStop-fCouplingStart))/fCouplingStep);
  Int_t fNevents[fN];
  Double_t fSumAll[fN];
  Double_t fSumGood[fN];
  Double_t fAcc[fN];
  Double_t fGammaTot[fN];
  Double_t fTau[fN];
  Double_t fProb[fN];
  Double_t fYield[fN];
  Double_t fCouplings[fN];
  ofstream ExclusionFile;

  TH2D *fhReach;
  TH2D *fhDecay;
  TH2D *fhWeight;

  TGraph *fgAcc;
  TGraph *fgYield;
  TGraph *fgGammaTot;
  TGraph *fgTau;
};

#endif
