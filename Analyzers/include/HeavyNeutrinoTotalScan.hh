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

  // NA62 parameters                                                                                    

  Double_t fpMom;
  Double_t fBeA;
  Double_t fBeDensity;
  Double_t fpBeLambda;
  Double_t ftargetLength;
  Double_t fCuA;
  Double_t fCuDensity;
  Double_t fpCuLambda;
  Double_t fTAXLength;
  Double_t fTAXDistance;
  Double_t fbeamLength;
  Double_t fzCHOD;
  Double_t fzMUV3;
  Double_t fLFV;
  Double_t fLInitialFV;
  Double_t frMinStraw;
  Double_t frMaxStraw;
  Double_t fzCHODPlane;
  Double_t frMinCHOD;
  Double_t frMaxCHOD;
  Double_t fzStraw[4];
  Double_t fxStrawChamberCentre[4];

  // Other parameters                                                                                   

  Double_t fDBeProdProb;
  Double_t fDCuProdProb;
  Double_t fDDecayProb;
  TwoLinesCDA *fCDAcomp;
  PointLineDistance *fDistcomp;
  LAVMatching *fLAVMatching;
  SAVMatching *fSAVMatching;

  // Scan variables

  Int_t fCouplingIndex;
  Int_t fMassIndex;
  Int_t fCounter;
  Double_t fMN;
  Double_t fTemp;
  static constexpr Double_t fCouplingStart = -10.;
  static constexpr Double_t fCouplingStop = 0.;
  static constexpr Double_t fCouplingStep = 0.1;
  static constexpr Int_t fNcoupling = round((std::abs(fCouplingStop-fCouplingStart))/fCouplingStep);
  static constexpr Double_t fInitialMass = 100.;
  static constexpr Double_t fFinalMass = 2000.;
  static constexpr Double_t fMassStep = 100.;
  static constexpr Int_t fNmass = round((std::abs(fFinalMass-fInitialMass))/fMassStep);
  Int_t fNevents[fNcoupling][fNmass];
  Double_t fSumAll[fNcoupling][fNmass];
  Double_t fSumGood[fNcoupling][fNmass];
  Double_t fAcc[fNcoupling][fNmass];
  Double_t fYield[fNcoupling][fNmass];
  Double_t fProb[fNcoupling][fNmass];
  Double_t fCouplings[fNcoupling];
  Double_t fMasses[fNmass];

  // Histos

  TGraph *fgExclusion;
};

#endif
