#ifndef HEAVYNEUTRINO_HH
#define HEAVYNEUTRINO_HH

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

class HeavyNeutrino : public NA62Analysis::Analyzer {

public:

  HeavyNeutrino(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrino();
  void InitHist();
  void InitOutput() {}
  void DefineMCSimple() {}
  void Process(Int_t);
  void StartOfBurstUser() {}
  void EndOfBurstUser();
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
  Double_t fNevents;
  Double_t fSumAll;
  Double_t fSumGood;
  TwoLinesCDA *fCDAcomp;
  PointLineDistance *fDistcomp;
  LAVMatching *fLAVMatching;
  SAVMatching *fSAVMatching;
  std::vector<Int_t> fID;
  std::vector<std::string> fStream;

  // Histos

  TH1D *fhNk3pi;
  TH1D *fhNbursts;
  TH1D *fhNEvents;
  TH1D *fhN2tracks;
  TH1D *fhNtracks;

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
  TH1D *fhMomPi;
  TH1D *fhMomMu;

  TH2D *fhXYSpec0Reco;
  TH2D *fhXYSpec1Reco;
  TH2D *fhXYSpec2Reco;
  TH2D *fhXYSpec3Reco;
  TH2D *fhXYCHODReco;
  TH2D *fhXYCHODTrue;
  TH2D *fhXYMUV3True;
  TH2D *fhP1vsP2;

  TH1D *fhPhysicsEventsVsCuts;

  TH2D *fhCDAvsZVertex_TotMomToBeamlineInitial;
  TH2D *fhCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack;
  TH2D *fhCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts;
  TH2D *fhCDAvsZVertex_TotMomToBeamlineAfterGeomCuts;
  TH2D *fhCDAvsZVertex_TotMomToBeamlineAfterVetoes;
  TH2D *fhCDAvsZVertex_TotMomToBeamlineFinal;

  TH2D *fhZvertexvsBeamlineDistInitial;
  TH2D *fhZvertexvsBeamlineDistAfterDownstreamTrack;
  TH2D *fhZvertexvsBeamlineDistAfterEnergyCuts;
  TH2D *fhZvertexvsBeamlineDistAfterGeomCuts;
  TH2D *fhZvertexvsBeamlineDistAfterVetoes;
  TH2D *fhZvertexvsBeamlineDistFinal;

  TH2D *fhCDAvsZVertex_TrackToBeamlineInitial;
  TH2D *fhCDAvsZVertex_TrackToTrackInitial;
  TH2D *fhCDAvsZVertex_TrackToBeamlineAfterCut;
  TH2D *fhCDAvsZVertex_TrackToTrackAfterCut;
  TH2D *fhCDAvsZVertex_TrackToBeamlineFinal;
  TH2D *fhCDAvsZVertex_TrackToTrackFinal;

  TH2D *fhBeamlineDistvsTargetDist_TotMomInitial;
  TH2D *fhBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack;
  TH2D *fhBeamlineDistvsTargetDist_TotMomAfterEnergyCuts;
  TH2D *fhBeamlineDistvsTargetDist_TotMomAfterGeomCuts;
  TH2D *fhBeamlineDistvsTargetDist_TotMomAfterVetoes;
  TH2D *fhBeamlineDistvsTargetDist_TotMomFinal;

  TH1D *fhDeltaTimeFromCHOD;
  TH1D *fhNMUV3CandAssocToTrack;
  TH1D *fhNCHODCandAssocToTrack;

  TH1D *fhEoP;
  TH2D *fhEoPMuVsPi;

  TH2D *fhSingleAddEnLKrHit;
  TH2D *fhSingleAddEnLKrCand;
  TH2D *fhAddEnLKrHit;
  TH2D *fhAddEnLKrCand;

  TH1D *fhInvMassMC;
  TH1D *fhInvMassReco;

  TH1D *fhAcc;
  TH1D *fhYield;
};

#endif
