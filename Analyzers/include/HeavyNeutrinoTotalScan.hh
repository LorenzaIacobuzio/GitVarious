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

  Double_t ComputeHNLMass(KinePart*);
  Double_t ComputeL(TVector3, TVector3, TVector3);
  Double_t ComputeNDecayProb(KinePart*, Double_t, Double_t);
  Double_t ComputeNReachProb(KinePart*, Double_t, Double_t);
  Double_t PhaseSpace(Double_t, Double_t, Double_t);
  Double_t PhaseSpaceFactor(Double_t, Double_t, Double_t);
  Double_t TwoBodyBR(Double_t, Double_t, Double_t, Int_t, Bool_t);
  Double_t ThreeBodyBR(Double_t, Double_t, Double_t, Double_t, Int_t, Bool_t);
  std::string ThreeBodyFunction(Double_t, Double_t);
  Double_t Gamma2(Double_t, Double_t, Double_t, Double_t, Bool_t);                          
  Double_t GammaLeptonNu3(Double_t, Double_t, Double_t, Bool_t);                      
  Double_t GammaTot(Double_t);                                               
  Double_t tauN(Double_t);                                                                             
  Double_t lambda(Double_t, Double_t, Double_t);
  Double_t ComputeProd(KinePart*, Double_t);
  Double_t ComputeDecay(Double_t);

protected:

  // Histos

  TGraph *fgExclusion;

  // Masses                                                                                        

  Double_t fMe;
  Double_t fMmu;
  Double_t fMtau;
  Double_t fMpi;
  Double_t fMpi0;
  Double_t fMrho;
  Double_t fMrho0;
  Double_t fMeta;
  Double_t fMetaprime;
  Double_t fMD;
  Double_t fMDS;
  Double_t fMD0;
  Double_t fMK;
  Double_t fMK0;
  Double_t fMp;
  Double_t fMKStar;
  Double_t fMK0Star;

  // Lifetimes                                                                                          

  Double_t fDlife;
  Double_t fDSlife;
  Double_t fD0life;
  Double_t ftaulife;

  // Constants                                                                                          

  Double_t fhc;
  Double_t fcLight;
  Double_t fGF;
  Double_t fPi;
  Double_t fRho;
  Double_t fD;
  Double_t fDS;
  Double_t fK;
  Double_t fEta;
  Double_t fEtaprime;
  Double_t fsigmacc;

  // CKM                                                                                              

  Double_t fVcs;
  Double_t fVcd;
  Double_t fVud;
  Double_t fVus;

  // Form factors, pseudoscalar and vector mesons                                                       

  Double_t fDK0;
  Double_t fDpi0;
  Double_t fD0K;
  Double_t fD0pi;
  Double_t fgDK0;
  Double_t fgDpi0;
  Double_t fgD0K;
  Double_t fgD0pi;
  Double_t fA0D;
  Double_t fA1D;
  Double_t fA2D;
  Double_t fVD;
  Double_t fA0D0;
  Double_t fA1D0;
  Double_t fA2D0;
  Double_t fVD0;

  // Fragmentation fractions                                                                            

  Double_t ffD;
  Double_t ffD0;
  Double_t ffDS;

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
  Double_t* fzStraw;
  Double_t* fxStrawChamberCentre;

  // Other parameters                                                                                   

  Double_t fDBeProdProb;
  Double_t fDCuProdProb;
  Double_t fDDecayProb;
  Double_t fUSquared;
  Double_t fUeSquared;
  Double_t fUmuSquared;
  Double_t fUtauSquared;
  Double_t fUeSquaredRatio;
  Double_t fUmuSquaredRatio;
  Double_t fUtauSquaredRatio;

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
  Double_t fCouplings[fNcoupling];
  Double_t fMasses[fNmass];

  // Other variables

  TwoLinesCDA *fCDAcomp;
  PointLineDistance *fDistcomp;
  LAVMatching *fLAVMatching;
  SAVMatching *fSAVMatching;
};

#endif
