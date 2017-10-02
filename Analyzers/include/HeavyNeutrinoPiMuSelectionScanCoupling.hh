#ifndef HEAVYNEUTRINOPIMUSELECTIONSCANCOUPLING_HH
#define HEAVYNEUTRINOPIMUSELECTIONSCANCOUPLING_HH

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

class HeavyNeutrinoPiMuSelectionScanCoupling : public NA62Analysis::Analyzer {

public:

  HeavyNeutrinoPiMuSelectionScanCoupling(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrinoPiMuSelectionScanCoupling();
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
  Double_t fDBeProdProb;
  Double_t fDCuProdProb;
  Double_t fDDecayProb;
  Double_t fNPiMu;
  Double_t fNevents;
  Double_t fCoupling;
  Double_t fCouplingStart = 1.E-9;
  Double_t fCouplingStop = 1.E-5;
  Double_t fCouplingStep = 0.1;
  Int_t fN = TMath::Abs((fCouplingStop-fCouplingStart)/fCouplingStep);
  Double_t fSumAll[fN];
  Double_t fSumGood[fN];
  Double_t fAcc[fN];
  Double_t fProb[fN];
  Double_t fYield[fN];
};

#endif
