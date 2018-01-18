#ifndef HEAVYNEUTRINOPIMUSELECTION_HH
#define HEAVYNEUTRINOPIMUSELECTION_HH

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

class HeavyNeutrinoPiMuSelection : public NA62Analysis::Analyzer {

public:

  HeavyNeutrinoPiMuSelection(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrinoPiMuSelection();
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
  Double_t fcos2ThetaC;
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

  // Other parameters                                                                                   

  Int_t fNSpecies;
  Double_t* fPTotal;
  Double_t* fmesonMass;
  Double_t* fmesonTau;
  Double_t fDBeProdProb;
  Double_t fDCuProdProb;
  Double_t fDDecayProb;
  Double_t fUSquared;
  Double_t fUeSquared;
  Double_t fUmuSquared;
  Double_t fUtauSquared;
  TwoLinesCDA *fCDAcomp;
  PointLineDistance *fDistcomp;
  LAVMatching *fLAVMatching;
  SAVMatching *fSAVMatching;
  Double_t fNevents;
  Double_t fSumAll;
  Double_t fSumGood;
};

#endif
