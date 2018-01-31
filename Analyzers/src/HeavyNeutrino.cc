#include <stdlib.h>
#include <iostream>
#include <string>
#include <cmath>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include "HeavyNeutrino.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "GeometricAcceptance.hh"
#include "TriggerConditions.hh"
#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/GaussLegendreIntegrator.h"
#include "TF1.h"
#include "TF2.h"
#include "HNLFunctions.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

#define TRIGGER_L0_PHYSICS_TYPE 1
#define MCTriggerMask 0xFF
#define LabelSize 0.05

/// \class HeavyNeutrino

HeavyNeutrino::HeavyNeutrino(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrino") {

  RequestAllMCTrees();
  RequestAllRecoTrees();
  RequestL0Data();
  RequestL1Data();

  AddParam("USquared", &fUSquared, 1.E-6);
  AddParam("UeSquaredRatio", &fUeSquaredRatio, 1.);
  AddParam("UmuSquaredRatio", &fUmuSquaredRatio, 16.);
  AddParam("UtauSquaredRatio", &fUtauSquaredRatio, 3.8);

  // Other parameters                                                                                   

  fSumGood     = 0.;
  fNevents     = 0.;
  fSumAll      = 0.;
  fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
  fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
  fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;
  fCDAcomp     = new TwoLinesCDA();
  fDistcomp    = new PointLineDistance();
  fLAVMatching = new LAVMatching();
  fSAVMatching = new SAVMatching();

  // Parameters for L0 trigger conditions

  fStream = {"RICH-Q2-MO1", "RICH-Q2-M1", "RICH-Q2-MO1-LKr10", "RICH-Q2-M1-LKr20", "RICH-Q2-MO2-nLKr20", "RICH-Q2-MO2", "RICH-Q2-M2", "RICH-QX-LKr20", "RICH-LKr20", "RICH-Q2-nMUV-LKr20", "RICH-Q2-MO1-LKr20",  "RICH-Q2-MO2-nLKr30"};

  for (UInt_t i = 0; i < fStream.size(); i++) {
    fID.push_back(TriggerConditions::GetInstance()->GetL0TriggerID(fStream[i]));
  }

  // Histos

  fhNk3pi    = nullptr;
  fhNbursts  = nullptr;
  fhNEvents  = nullptr;
  fhN2tracks = nullptr;
  fhNtracks  = nullptr;

  fhZDProd       = nullptr;
  fhZDDecay      = nullptr;
  fhDTheta       = nullptr;
  fhDLambda      = nullptr;
  fhDPath        = nullptr;
  fhDMom         = nullptr;
  fhZHNLDecay    = nullptr;
  fhHNLGamma     = nullptr;
  fhHNLDecayProb = nullptr;
  fhHNLReachProb = nullptr;
  fhHNLTheta     = nullptr;
  fhHNLMom       = nullptr;
  fhWeight       = nullptr;
  fhMomPi        = nullptr;
  fhMomMu        = nullptr;

  fhXYSpec0Reco = nullptr;
  fhXYSpec1Reco = nullptr;
  fhXYSpec2Reco = nullptr;
  fhXYSpec3Reco = nullptr;
  fhXYCHODReco  = nullptr;
  fhXYCHODTrue  = nullptr;
  fhXYMUV3True  = nullptr;
  fhP1vsP2      = nullptr;

  fhPhysicsEventsVsCuts = nullptr;

  fhCDAvsZVertex_TotMomToBeamlineInitial              = nullptr;
  fhCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack = nullptr;
  fhCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts      = nullptr;
  fhCDAvsZVertex_TotMomToBeamlineAfterGeomCuts        = nullptr;
  fhCDAvsZVertex_TotMomToBeamlineAfterVetoes          = nullptr;
  fhCDAvsZVertex_TotMomToBeamlineFinal                = nullptr;

  fhZvertexvsBeamlineDistInitial              = nullptr;
  fhZvertexvsBeamlineDistAfterDownstreamTrack = nullptr;
  fhZvertexvsBeamlineDistAfterEnergyCuts      = nullptr;
  fhZvertexvsBeamlineDistAfterGeomCuts        = nullptr;
  fhZvertexvsBeamlineDistAfterVetoes          = nullptr;
  fhZvertexvsBeamlineDistFinal                = nullptr;

  fhCDAvsZVertex_TrackToBeamlineInitial  = nullptr;
  fhCDAvsZVertex_TrackToTrackInitial     = nullptr;
  fhCDAvsZVertex_TrackToBeamlineAfterCut = nullptr;
  fhCDAvsZVertex_TrackToTrackAfterCut    = nullptr;
  fhCDAvsZVertex_TrackToBeamlineFinal    = nullptr;
  fhCDAvsZVertex_TrackToTrackFinal       = nullptr;

  fhBeamlineDistvsTargetDist_TotMomInitial              = nullptr;
  fhBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack = nullptr;
  fhBeamlineDistvsTargetDist_TotMomAfterEnergyCuts      = nullptr;
  fhBeamlineDistvsTargetDist_TotMomAfterGeomCuts        = nullptr;
  fhBeamlineDistvsTargetDist_TotMomAfterVetoes          = nullptr;
  fhBeamlineDistvsTargetDist_TotMomFinal                = nullptr;

  fhDeltaTimeFromCHOD     = nullptr;
  fhNMUV3CandAssocToTrack = nullptr;
  fhNCHODCandAssocToTrack = nullptr;

  fhEoP       = nullptr;
  fhEoPMuVsPi = nullptr;

  fhSingleAddEnLKrHit  = nullptr;
  fhSingleAddEnLKrCand = nullptr;
  fhAddEnLKrHit        = nullptr;
  fhAddEnLKrCand       = nullptr;

  fhInvMassMC   = nullptr;
  fhInvMassReco = nullptr;

  fhAcc   = nullptr;
  fhYield = nullptr;    
}

void HeavyNeutrino::InitHist() {

  BookHisto("hNk3pi",    new TH1D("Nk3pi",    "Total number of K3pi events",       1, 0., 1.));
  BookHisto("hNbursts",  new TH1D("Nbursts",  "Total number of processed bursts",  1, 0., 1.));
  BookHisto("hNEvents",  new TH1D("NEvents",  "Number of total processed events" , 1, 0., 1.));
  BookHisto("hNtracks",  new TH1D("Ntracks",  "Number of tracks",                  4, -0.5, 3.5));
  BookHisto("hN2tracks", new TH1D("N2tracks", "Number of two-tracks events",       1, 0., 1.));

  BookHisto("hZDProd",        new TH1D("ZDProd", "Z of D meson production point", 20000., -250., 33000.));
  BookHisto("hZDDecay",       new TH1D("ZDDecay", "Z of D meson decay point",     20000., -250., 33000.));
  BookHisto("hDTheta",        new TH1D("DTheta",     "D meson theta",              100,  0., 0.3));
  BookHisto("hDLambda",       new TH1D("DLambda",    "D meson decay length",       100, -1., 40.));
  BookHisto("hDPath",         new TH1D("DPath",      "D meson path in Z",          100, -1., 50.));
  BookHisto("hDMom",          new TH1D("DMom",       "D meson momentum",           100, -1., 170.));

  BookHisto("hZHNLDecay",     new TH1D("ZHNLDecay",    "Z of HNL decay point",           100., 90., 190.));
  BookHisto("hHNLGamma",      new TH1D("HNLGamma",     "Lorentz gamma of HNL",           50., 0., 170.));
  BookHisto("hHNLDecayProb",  new TH1D("HNLDecayProb", "HNL decay probability",          100., 0., 0.0065));
  BookHisto("hHNLReachProb",  new TH1D("HNLReachProb", "HNL probability of reaching FV", 100., 0.99, 1.001));
  BookHisto("hHNLTheta",      new TH1D("HNLTheta",     "HNL theta",                      100., 0., 0.5));
  BookHisto("hHNLMom",        new TH1D("HNLMom",       "HNL momentum",                   100., -0.5, 200.));
  BookHisto("hWeight",        new TH1D("Weight",       "Weight",                         1000, 0., 0.08E-12));
  BookHisto("hMomPi",         new TH1D("MomPi",        "Pion momentum",                  100, -0.5, 200.));
  BookHisto("hMomMu",         new TH1D("MomMu",        "Muon momentum",                  100, -0.5, 200.));

  BookHisto("hXYSpec0Reco", new TH2D("XYSpec0Reco",        "Two-track reconstructed events at CH1",  100, -1.5, 1.5, 100, -1.5, 1.5)); 
  BookHisto("hXYSpec1Reco", new TH2D("XYSpec1Reco",        "Two-track reconstructed events at CH2",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYSpec2Reco", new TH2D("XYSpec2Reco",        "Two-track reconstructed events at CH3",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYSpec3Reco", new TH2D("XYSpec3Reco",        "Two-track reconstructed events at CH4",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYCHODReco",  new TH2D("XYCHODReco",         "Two-track reconstructed events at CHOD", 100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYCHODTrue",  new TH2D("XYCHODTrue",         "X,Y of HNL daughters at CHOD, from MC",  100, -2., 2., 100, -2., 2.));
  BookHisto("hXYMUV3True",  new TH2D("XYMUV3True",         "X,Y of HNL daughters at MUV3, from MC",  100, -2., 2., 100, -2., 2.));
  BookHisto("hP1vsP2",      new TH2D("P1vsP2",      "Trimomentum of the two HNL daughters, from MC", 100, 0., 100., 100, 0., 100.));

  BookHisto("hPhysicsEventsVsCuts", new TH1D("PhysicsEventsVsCuts", "Physics events passing the selection cuts", 35, 0., 35.));

  BookHisto("hCDAvsZVertex_TotMomToBeamlineInitial",              new TH2D("CDAvsZVertex_TotMomToBeamlineInitial",              "Two-track total momentum wrt beam axis, before all cuts",            100, 0., 240., 100, 0., 0.5));
  BookHisto("hCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack", new TH2D("CDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack", "Two-track total momentum wrt beam axis, after downstream selection", 100, 0., 240., 100, 0., 0.5));
  BookHisto("hCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts",      new TH2D("CDAvsZVertex_TotMomToBeamlineAfterEnergyCuts",      "Two-track total momentum wrt beam axis, after energy cuts",          100, 0., 240., 100, 0., 0.5));
  BookHisto("hCDAvsZVertex_TotMomToBeamlineAfterVetoes",          new TH2D("CDAvsZVertex_TotMomToBeamlineAfterVetoes",          "Two-track total momentum wrt beam axis, after vetoes",               100, 0., 240., 100, 0., 0.5));
  BookHisto("hCDAvsZVertex_TotMomToBeamlineAfterGeomCuts",        new TH2D("CDAvsZVertex_TotMomToBeamlineAfterGeomCuts",        "Two-track total momentum wrt beam axis, after geometrical cuts",     100, 0., 240., 100, 0., 0.5));
  BookHisto("hCDAvsZVertex_TotMomToBeamlineFinal",                new TH2D("CDAvsZVertex_TotMomToBeamlineFinal",                "Two-track total momentum wrt beam axis, after all selections",       100, 0., 240., 100, 0., 0.5));
  
  BookHisto("hZvertexvsBeamlineDistInitial",              new TH2D("ZvertexvsBeamlineDistInitial",              "Two track vertex wrt beam axis, before all cuts",            50, 0., 240., 100, 0., 1.));
  BookHisto("hZvertexvsBeamlineDistAfterDownstreamTrack", new TH2D("ZvertexvsBeamlineDistAfterDownstreamTrack", "Two track vertex wrt beam axis, after downstream selection", 50, 0., 240., 100, 0., 1.));
  BookHisto("hZvertexvsBeamlineDistAfterEnergyCuts",      new TH2D("ZvertexvsBeamlineDistAfterEnergyCuts",      "Two track vertex wrt beam axis, after energy cuts",          50, 0., 240., 100, 0., 1.));
  BookHisto("hZvertexvsBeamlineDistAfterVetoes",          new TH2D("ZvertexvsBeamlineDistAfterVetoes",          "Two track vertex wrt beam axis, after vetoes",               50, 0., 240., 100, 0., 1.));
  BookHisto("hZvertexvsBeamlineDistAfterGeomCuts",        new TH2D("ZvertexvsBeamlineDistAfterGeomCuts",        "Two track vertex wrt beam axis, after geometrical cuts",     50, 0., 240., 100, 0., 1.));
  BookHisto("hZvertexvsBeamlineDistFinal",                new TH2D("ZvertexvsBeamlineDistFinal",                "Two track vertex wrt beam axis, after all selections",       50, 0., 240., 100, 0., 1.));

  BookHisto("hCDAvsZVertex_TrackToBeamlineInitial",  new TH2D("CDAvsZVertex_TrackToBeamlineInitial",  "Track wrt beam axis, before all cuts",      100, 0., 240., 100, 0., 0.2));
  BookHisto("hCDAvsZVertex_TrackToBeamlineAfterCut", new TH2D("CDAvsZVertex_TrackToBeamlineAfterCut", "Track wrt beam axis, after cut",            100, 0., 240., 100, 0., 1.));
  BookHisto("hCDAvsZVertex_TrackToBeamlineFinal",    new TH2D("CDAvsZVertex_TrackToBeamlineFinal",    "Track wrt beam axis, after all selections", 100, 0., 240., 100, 0., 1.));
  BookHisto("hCDAvsZVertex_TrackToTrackInitial",     new TH2D("CDAvsZVertex_TrackToTrackInitial",     "Track1 wrt track2, before all cuts",        100, 0., 240., 100, 0., 0.2));
  BookHisto("hCDAvsZVertex_TrackToTrackAfterCut",    new TH2D("CDAvsZVertex_TrackToTrackAfterCut",    "Track1 wrt track2, after cut",              100, 0., 240., 100, 0., 1.));
  BookHisto("hCDAvsZVertex_TrackToTrackFinal",       new TH2D("CDAvsZVertex_TrackToTrackFinal",       "Track1 wrt track2, after all selections",   100, 0., 240., 100, 0., 1.));
  
  BookHisto("hBeamlineDistvsTargetDist_TotMomInitial",              new TH2D("BeamlineDistvsTargetDist_TotMomInitial",              "Two-track total momentum wrt beam axis, before all cuts",            100, 0., 1., 100, 0., 1.));
  BookHisto("hBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack", new TH2D("BeamlineDistvsTargetDist_TotMomAfterDownstreamTrack", "Two-track total momentum wrt beam axis, after downstream selection", 100, 0., 1., 100, 0., 1.));
  BookHisto("hBeamlineDistvsTargetDist_TotMomAfterEnergyCuts",      new TH2D("BeamlineDistvsTargetDist_TotMomAfterEnergyCuts",      "Two-track total momentum wrt beam axis, after energy cuts",          100, 0., 1., 100, 0., 1.));
  BookHisto("hBeamlineDistvsTargetDist_TotMomAfterVetoes",          new TH2D("BeamlineDistvsTargetDist_TotMomAfterVetoes",          "Two-track total momentum wrt beam axis, after vetoes",               100, 0., 1., 100, 0., 1.));
  BookHisto("hBeamlineDistvsTargetDist_TotMomAfterGeomCuts",        new TH2D("BeamlineDistvsTargetDist_TotMomAfterGeomCuts",        "Two-track total momentum wrt beam axis, after geometrical cuts",     100, 0., 1., 100, 0., 1.));
  BookHisto("hBeamlineDistvsTargetDist_TotMomFinal",                new TH2D("BeamlineDistvsTargetDist_TotMomFinal",                "Two-track total momentum wrt beam axis, after all selections",       100, 0., 1., 100, 0., 1.));
  
  BookHisto("hDeltaTimeFromCHOD",     new TH1D("DeltaTimeFromCHOD",     "Time difference of two tracks (CHOD candidates)",    200, -20., 20.));  
  BookHisto("hNMUV3CandAssocToTrack", new TH1D("NMUV3CandAssocToTrack", "Number of MUV3 candidates associated to each track", 4, -0.5, 3.5));
  BookHisto("hNCHODCandAssocToTrack", new TH1D("NCHODCandAssocToTrack", "Number of CHOD candidates associated to each track", 10, 0., 10.));
  
  BookHisto("hEoP", new TH1D("EoP", "E/p in LKr", 100, 0., 1.2));
  BookHisto("hEoPMuVsPi", new TH2D("EoPMuVsPi", "Muon E/p vs pion E/p in LKr", 100, 0., 1.2, 100, 0., 0.1));
  
  BookHisto("hSingleAddEnLKrHit",  new TH2D("SingleAddEnLKrHit", "Single hit energy vs time, for additional LKr hits",              250, -100., 100, 250, 0., 100.));
  BookHisto("hSingleAddEnLKrCand", new TH2D("SingleAddEnLKrCand", "Single candidate energy vs time, for additional LKr candidates", 250, -100., 100, 250, 0., 100.));
  BookHisto("hAddEnLKrHit",        new TH2D("AddEnLKrHit", "Additional energy vs time, for LKr hits",                               250, -100., 100, 250, 0., 10.));
  BookHisto("hAddEnLKrCand",       new TH2D("AddEnLKrCand", "Additional energy vs time, for LKr candidates",                        250, -100., 100, 250, 0., 10.));
  
  BookHisto("hInvMassMC",   new TH1D("InvMassMC",   "Invariant mass MC", 100, 999., 1001.));
  BookHisto("hInvMassReco", new TH1D("InvMassReco", "Invariant mass Reco", 50, 960., 1040.));

  BookHisto("hAcc",   new TH1D("Acc",   "Acceptance", 50, 1.E-6, 5.E-5));
  BookHisto("hYield", new TH1D("Yield", "Yield",      50, 1.E-20, 1.E-18));
}

void HeavyNeutrino::Process(Int_t) {

  //TRecoLKrEvent*          LKrEvent  = (TRecoLKrEvent*)          GetEvent("LKr");
  TRecoLAVEvent*          LAVEvent  = (TRecoLAVEvent*)          GetEvent("LAV");
  TRecoIRCEvent*          IRCEvent  = (TRecoIRCEvent*)          GetEvent("IRC");
  TRecoSACEvent*          SACEvent  = (TRecoSACEvent*)          GetEvent("SAC");

  // Counter for cuts

  Int_t CutID = 0;
  
  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  Double_t MN             = 0.;
  Double_t HNLTau         = 0.;
  Double_t gammaTot       = 0.;
  Double_t NDecayProb     = 0.;
  Double_t NReachProb     = 0.;
  Double_t LReach         = 0.;
  Double_t LeptonUSquared = 0.;
  Double_t ProdFactor     = 0.;
  Double_t DecayFactor    = 0.;
  Double_t Weight         = 0.;
  Double_t DProdProb      = 0.;
  Int_t counter           = 0;
  TLorentzVector mom1;
  TLorentzVector mom2;
  TVector3 point1;
  TVector3 point2;
  TVector3 momentum1;
  Double_t p1,p2;

  // Some plots of KinePart quantities

  if (GetWithMC()) {
    Event *evt = GetMCEvent();
    counter = 0;
    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);
      if (p->GetParentID() == 0) {
	counter++;
	Double_t xCHOD = p->xAt(fzCHOD);                      
	Double_t yCHOD = p->yAt(fzCHOD); 
	Double_t xMUV3 = p->xAt(fzMUV3);
	Double_t yMUV3 = p->yAt(fzMUV3);
	FillHisto("hXYCHODTrue", xCHOD/1000., yCHOD/1000.);
	FillHisto("hXYMUV3True", xMUV3/1000., yMUV3/1000.);
	if (p->GetParticleName() == "pi+" || p->GetParticleName() == "pi-") {
	  mom1 = p->GetInitial4Momentum();
	  p1 = TMath::Sqrt(mom1.Px()*mom1.Px()+mom1.Py()*mom1.Py()+mom1.Pz()*mom1.Pz());
	}
	else if (p->GetParticleName() == "mu+" || p->GetParticleName() == "mu-") {
	  mom2 = p->GetInitial4Momentum();
	  p2 = TMath::Sqrt(mom2.Px()*mom2.Px()+mom2.Py()*mom2.Py()+mom2.Pz()*mom2.Pz());
	}
	if (counter == 2) {
	  FillHisto("hMomPi",  p1/1000.);
	  FillHisto("hMomMu",  p2/1000.);
	  FillHisto("hP1vsP2", p1/1000., p2/1000.);
	  FillHisto("hInvMassMC", (mom1+mom2).M2()/1000.);
	}
      }
    }

    // Computation of coupling-related quantities of all HNLs (good and bad)

    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);      
      if (p->GetParentID() == -1 && p->GetPDGcode() == 999) {
	fNevents++;
	point1.SetXYZ(p->GetProdPos().X(), p->GetProdPos().Y(), p->GetProdPos().Z());
	point2.SetXYZ(0., 0., fLInitialFV);
	momentum1.SetXYZ(p->GetInitial4Momentum().Px(), p->GetInitial4Momentum().Py(), p->GetInitial4Momentum().Pz());
	MN = ComputeHNLMass(p);
	gammaTot = GammaTot(MN);
	HNLTau = tauN(MN);
	LReach = ComputeL(point1, point2, momentum1);
	ProdFactor = ComputeProd(p, MN);
	DecayFactor = ComputeDecay(MN);
	NReachProb = ComputeNReachProb(p, HNLTau, LReach);
	NDecayProb = ComputeNDecayProb(p, HNLTau, fLFV);
	if (p->GetParticleName().Contains("e") && !p->GetParticleName().Contains("nu_tau"))
	  LeptonUSquared = fUeSquared;
	else if (p->GetParticleName().Contains("mu") && !p->GetParticleName().Contains("nu_tau"))
	  LeptonUSquared = fUmuSquared;
	else if (p->GetParticleName() == "DS->Ntau" || p->GetParticleName().Contains("nu_tau") || (p->GetParticleName().Contains("tau") && (p->GetParticleName().Contains("rho") || p->GetParticleName().Contains("pi")|| p->GetParticleName().Contains("K"))))
	  LeptonUSquared = fUtauSquared;

	if (p->GetProdPos().Z() < fTAXDistance)
	  DProdProb = fDBeProdProb;
	else if (p->GetProdPos().Z() >= fTAXDistance)
	  DProdProb = fDCuProdProb;

	// Weight to be associated to each HNL

	Weight = DProdProb*fDDecayProb*NReachProb*NDecayProb*DecayFactor*ProdFactor*LeptonUSquared;
	fSumAll += Weight;

	// Some more plots of KinePart quantities

	FillHisto("hZDProd", p->GetPosAtCheckPoint(0).z());
	FillHisto("hZDDecay", p->GetProdPos().Z());
	FillHisto("hDTheta", p->GetPosAtCheckPoint(0).x());
	FillHisto("hDLambda", p->GetPosAtCheckPoint(0).y());
	FillHisto("hDPath", p->GetMomAtCheckPoint(0).X());
	FillHisto("hDMom", p->GetMomAtCheckPoint(0).Y()/1000.);
	FillHisto("hZHNLDecay", p->GetEndPos().Z()/1000.);
	FillHisto("hHNLGamma", p->GetInitial4Momentum().Gamma());
	FillHisto("hHNLReachProb", NReachProb);
	FillHisto("hHNLDecayProb", NDecayProb);
	FillHisto("hHNLTheta", p->GetMomAtCheckPoint(0).Z());
	FillHisto("hHNLMom", p->GetMomAtCheckPoint(0).T()/1000.);
	FillHisto("hWeight",Weight);
      }
    }
  }

  // Compute number of processed events

  FillHisto("hNEvents", 0.5);
  
  // K3Pi
  
  Bool_t k3pi     = *(Bool_t*) GetOutput("K3piSelection.EventSelected");
  Int_t RunNumber = GetWithMC() ? 0 : GetEventHeader()->GetRunID();

  if (k3pi && 0x10)
    FillHisto("hNk3pi", 0.5);

  // L0 data
  
  L0TPData *L0TPData      = GetL0Data();
  UChar_t L0DataType      = GetWithMC() ? 0x11 : L0TPData->GetDataType();
  UInt_t L0TriggerFlags   = GetWithMC() ? 0xFF : L0TPData->GetTriggerFlags();
  Bool_t PhysicsTriggerOK = (L0DataType & TRIGGER_L0_PHYSICS_TYPE);
  Bool_t TriggerFlagsOK   = L0TriggerFlags & MCTriggerMask;
  Bool_t TriggerOK        = PhysicsTriggerOK && TriggerFlagsOK;

  if (!TriggerOK)
    return;
      
  // If real data
  // Select physics triggers and L0 trigger conditions

  if (!GetWithMC()) {  
    Bool_t L0OK = kFALSE;
    Bool_t L1OK = kFALSE;
    
    for (UInt_t i = 0; i < fID.size(); i++) {
      L0OK |= TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0TPData, fID[i]);
    }
    
    if (!L0OK) return;
    
    // Check for notKTAG o KTAG don't care at L1                                                      
    
    for (UInt_t i = 0; i < fID.size(); i++) {
      std::string L1Algo       = (std::string)TriggerConditions::GetInstance()->GetL1TriggerConditionName(RunNumber, fID[i]);
      std::size_t foundKTAG    = L1Algo.find("KTAG");
      std::size_t foundnotKTAG = L1Algo.find("notKTAG");
      if (L1Algo != "none") {
	if (foundKTAG == std::string::npos || foundnotKTAG != std::string::npos) {
	  L1OK = kTRUE;
	}
      }
    }
    
    if (!L1OK) return;
  }
  
  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  // Select two-track events
  
  std::vector<DownstreamTrack> Tracks = *(std::vector<DownstreamTrack>*) GetOutput("DownstreamTrackBuilder.Output");

  FillHisto("hNtracks", Tracks.size());
  
  if (Tracks.size() != 2)
    return;

  FillHisto("hN2tracks", 0.5);
  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  // Track features
    
  Int_t Charge1                                 = Tracks[0].GetCharge();
  Int_t Charge2                                 = Tracks[1].GetCharge();
  Double_t ChiSquare1                           = Tracks[0].GetChi2();
  Double_t ChiSquare2                           = Tracks[1].GetChi2();
  TRecoSpectrometerCandidate* SpectrometerCand1 = Tracks[0].GetSpectrometerCandidate();
  TRecoSpectrometerCandidate* SpectrometerCand2 = Tracks[1].GetSpectrometerCandidate();
  TVector3 Mom1                                 = SpectrometerCand1->GetThreeMomentumBeforeMagnet();
  TVector3 Mom2                                 = SpectrometerCand2->GetThreeMomentumBeforeMagnet();
  TVector3 TotMom                               = Mom1 + Mom2;
  
  // (X,Y) of reconstructed tracks for all Spectrometer chambers
    
  for (UInt_t i = 0; i < Tracks.size(); i++) {
    TRecoSpectrometerCandidate* Cand = Tracks[i].GetSpectrometerCandidate();
    for (Int_t j = 0; j < 4; j++) {
      Double_t x        = Cand->xAt(fzStraw[j]);
      Double_t y        = Cand->yAt(fzStraw[j]);
      Double_t r        = sqrt(x*x + y*y); 
      Double_t rShifted = sqrt(pow(x-fxStrawChamberCentre[j],2) + y*y); 
      if (rShifted > frMinStraw && r < frMaxStraw) {
	TString name = Form("hXYSpec%dReco",j);
	FillHisto(name, x/1000., y/1000.);
      }
    }
    Double_t x  = Cand->xAtAfterMagnet(fzCHODPlane);
    Double_t y  = Cand->yAtAfterMagnet(fzCHODPlane);
    Double_t r  = sqrt(x*x+y*y);
    Double_t r1 = frMinCHOD;
    Double_t r2 = frMaxCHOD;
    if (r > r1 && r < r2)
	FillHisto("hXYCHODReco", x/1000., y/1000.);
  }

  // Compute CDA of track1,2 wrt kaon axis, between each other and of two-track total momentum wrt beam axis
  
  fCDAcomp->SetLine1Point1(0.0, 0.0, 102000.0);     // kaon axis
  fCDAcomp->SetDir1(1.2E-3, 0., 1.);
  
  fCDAcomp->SetLine2Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);     // track1
  fCDAcomp->SetDir2(Mom1);
  fCDAcomp->ComputeVertexCDA();
  
  Double_t CDA1     = fCDAcomp->GetCDA();
  Double_t Zvertex1 = fCDAcomp->GetVertex().z();     // kaon axis-track1
  
  fCDAcomp->SetLine2Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand2->yAtBeforeMagnet(10.0), 10.0);     // track2
  fCDAcomp->SetDir2(Mom2);
  fCDAcomp->ComputeVertexCDA();
  
  Double_t CDA2     = fCDAcomp->GetCDA();
  Double_t Zvertex2 = fCDAcomp->GetVertex().z();     // kaon axis-track2
  
  fCDAcomp->SetLine1Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);     // track1
  fCDAcomp->SetDir1(Mom1);
  fCDAcomp->SetLine2Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand2->yAtBeforeMagnet(10.0), 10.0);     // track2
  fCDAcomp->SetDir2(Mom2);
  fCDAcomp->ComputeVertexCDA();
  
  Double_t CDA     = fCDAcomp->GetCDA();
  TVector3 Vertex  = fCDAcomp->GetVertex();
  Double_t Zvertex = fCDAcomp->GetVertex().z();     // track1-track2
  
  fCDAcomp->SetLine1Point1(0.0, 0.0, 102000.0);     // beam axis
  fCDAcomp->SetDir1(0., 0., 1.);
  
  fCDAcomp->SetLine2Point1(Vertex);     // total momentum of track1+track2
  fCDAcomp->SetDir2(TotMom);
  fCDAcomp->ComputeVertexCDA();
  
  Double_t CDAMom = fCDAcomp->GetCDA();
  
  // Compute distance of two-track momentum wrt target (extrapolation at target)
  
  fDistcomp->SetLineDir(TotMom);    
  fDistcomp->SetLinePoint1(Vertex);
  fDistcomp->SetPoint(0., 0., 0.);  
  fDistcomp->ComputeDistance();
  
  Double_t TargetDist = fDistcomp->GetDistance();
  
  // Compute distance of two-track vertex wrt beam axis
  
  fDistcomp->SetLinePoint1(0., 0., 102000.);
  fDistcomp->SetLineDir(0., 0., 1.);
  fDistcomp->SetPoint(Vertex);
  
  fDistcomp->ComputeDistance();
  
  Double_t BeamlineDist = fDistcomp->GetDistance();
  
  FillHisto("hCDAvsZVertex_TrackToBeamlineInitial",      Zvertex1 / 1000.,         CDA1 / 1000.);
  FillHisto("hCDAvsZVertex_TrackToBeamlineInitial",      Zvertex2 / 1000.,         CDA2 / 1000.);
  FillHisto("hCDAvsZVertex_TrackToTrackInitial",          Zvertex / 1000.,          CDA / 1000.);
  FillHisto("hCDAvsZVertex_TotMomToBeamlineInitial",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 1
  FillHisto("hZvertexvsBeamlineDistInitial",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 1
  FillHisto("hBeamlineDistvsTargetDist_TotMomInitial", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 1

  // Track selection, CUT 1: Two tracks in Spectrometer acceptance

  for (Int_t i = 0; i < 4; i++) {
      Double_t x1        = SpectrometerCand1->xAt(fzStraw[i]);
      Double_t y1        = SpectrometerCand1->yAt(fzStraw[i]);
      Double_t r1        = sqrt(x1*x1 + y1*y1);
      Double_t rShifted1 = sqrt(pow(x1-fxStrawChamberCentre[i],2) + y1*y1);
      Double_t x2        = SpectrometerCand2->xAt(fzStraw[i]);
      Double_t y2        = SpectrometerCand2->yAt(fzStraw[i]);
      Double_t r2        = sqrt(x2*x2 + y2*y2);
      Double_t rShifted2 = sqrt(pow(x2-fxStrawChamberCentre[i],2) + y2*y2);
      Bool_t inAcc       = false;

      if ((rShifted1 > frMinStraw && r1 < frMaxStraw) && (rShifted2 > frMinStraw && r2 < frMaxStraw))
	inAcc = true;
      if (!inAcc) 
	return;

      FillHisto("hPhysicsEventsVsCuts", CutID);
      CutID++;
  }

  // Track selection, CUT 2: Chi2 and momentum cuts

  if (ChiSquare1 >= 20. || ChiSquare2 >= 20.)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  if (SpectrometerCand1->GetNChambers() < 3 || SpectrometerCand2->GetNChambers() < 3)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  // Track selection, CUT 3: Opposite-charged tracks
  
  if (Charge1 + Charge2 != 0)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  // Plot the number of candidates associated to each track, for MUV3 and CHOD
  
  FillHisto("hNMUV3CandAssocToTrack", Tracks[0].GetNMUV3AssociationRecords());
  FillHisto("hNMUV3CandAssocToTrack", Tracks[1].GetNMUV3AssociationRecords());
  FillHisto("hNCHODCandAssocToTrack", Tracks[0].GetNCHODAssociationRecords());
  FillHisto("hNCHODCandAssocToTrack", Tracks[1].GetNCHODAssociationRecords());
  
  // Downstream track selection, CUT 4: Extrapolation and association to CHOD

  Bool_t CHODAssoc = (Tracks[0].CHODAssociationExists() && Tracks[1].CHODAssociationExists());
  Double_t x1      = SpectrometerCand1->xAtAfterMagnet(fzCHODPlane);
  Double_t y1      = SpectrometerCand1->yAtAfterMagnet(fzCHODPlane);
  Double_t r1      = sqrt(x1*x1+y1*y1);
  Double_t x2      = SpectrometerCand2->xAtAfterMagnet(fzCHODPlane);
  Double_t y2      = SpectrometerCand2->yAtAfterMagnet(fzCHODPlane);
  Double_t r2      = sqrt(x2*x2+y2*y2);
  Bool_t inAcc     = false;

  if ((r1 > frMinCHOD && r1 < frMaxCHOD) && (r2 > frMinCHOD && r2 < frMaxCHOD))
    inAcc = true;
  if (!inAcc)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  if (!CHODAssoc)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  // Downstream track selection, CUT 5: Extrapolation and association to LKr

  Bool_t LKrAssoc = (Tracks[0].LKrAssociationExists() && Tracks[1].LKrAssociationExists());

  if (!GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kLKr) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kLKr))
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  if (!LKrAssoc)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  // Downstream track selection, CUT 6: Extrapolation and association to MUV3
    
  if (!GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kMUV3) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kMUV3))
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  Bool_t Assoc1 = Tracks[0].MUV3AssociationExists();
  Bool_t Assoc2 = Tracks[1].MUV3AssociationExists();
  Int_t Assoc   = 0;
  
  if (Assoc1 && !Assoc2)
    Assoc = 1;
  else if (!Assoc1 && Assoc2)
    Assoc = 2;
  else
    Assoc = 0;

  if (Assoc != 1 && Assoc != 2)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  FillHisto("hCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 2
  FillHisto("hZvertexvsBeamlineDistAfterDownstreamTrack",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 2
  FillHisto("hBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 2
  
  // Compute time of MUV3 and CHOD candidates for better resolution wrt Spectrometer tracks
  
  Double_t MUV3Time;
  
  if (Assoc == 1)
    MUV3Time = Tracks[0].GetMUV3Time(0);
  else if (Assoc == 2)
    MUV3Time = Tracks[1].GetMUV3Time(0);
  
  Double_t CHODTime1 = Tracks[0].GetCHODTime();
  Double_t CHODTime2 = Tracks[1].GetCHODTime();
  
  // Plot time difference of the two tracks

  FillHisto("hDeltaTimeFromCHOD", CHODTime1 - CHODTime2);
  
  // Energy cuts, CUT 7: Cut on E/p in LKr
  
  Double_t EoP1  = Tracks[0].GetLKrEoP();
  Double_t EoP2  = Tracks[1].GetLKrEoP();
  Double_t MuEoP = 0.;
  Double_t PiEoP = 0.;

  if (Assoc == 1) {
    MuEoP = EoP1;
    PiEoP = EoP2;
  }
  else if (Assoc == 2) {
    MuEoP = EoP2;
    PiEoP = EoP1;
  }
  
  FillHisto("hEoP", EoP1);
  FillHisto("hEoP", EoP2);
  
  FillHisto("hEoPMuVsPi", PiEoP, MuEoP);

  if (MuEoP >= 0.2)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  if (PiEoP <= 0.2 || PiEoP >= 0.8)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  FillHisto("hCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 3
  FillHisto("hZvertexvsBeamlineDistAfterEnergyCuts",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 3
  FillHisto("hBeamlineDistvsTargetDist_TotMomAfterEnergyCuts", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 3

  // Veto cuts, CUT 8: LAV veto
  
  fLAVMatching->SetReferenceTime((CHODTime1 + CHODTime2) / 2);
  
  if (GetWithMC())
    fLAVMatching->SetTimeCuts(99999, 99999);     // LAV is not time-aligned in MC
  else
    fLAVMatching->SetTimeCuts(-10., 10.);
  
  if (fLAVMatching->LAVHasTimeMatching(LAVEvent))
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;
  
  // Veto cuts, CUT 9: SAV veto
  
  fSAVMatching->SetReferenceTime((CHODTime1 + CHODTime2) / 2);
  
  if (GetWithMC()) {
    fSAVMatching->SetIRCTimeCuts(99999, 99999);     // SAV is not time-aligned in MC
    fSAVMatching->SetSACTimeCuts(99999, 99999);
  } else {
    fSAVMatching->SetIRCTimeCuts(10.0, 10.0);
    fSAVMatching->SetSACTimeCuts(10.0, 10.0);
  }
  
  if (fSAVMatching->SAVHasTimeMatching(IRCEvent, SACEvent))
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;
  
  // Veto cuts, CUT 10: Residual LKr veto
  /*  
  Double_t LKrEnergyInTime = 0;
  Double_t LKrWeightedTime = 0;
  Int_t NLKrCand = LKrEvent->GetNCandidates();
  Int_t NLKrHits = LKrEvent->GetNHits();
  
  // Cut on candidates in time with CHOD, more than 20 cm distant from extrapolation points of both tracks, single candidate energy > 40 MeV
  
  for (Int_t i = 0; i < NLKrCand; i++) {
    TRecoLKrCandidate *LKrCand = (TRecoLKrCandidate*) LKrEvent->GetCandidate(i);
    Double_t X1 = LKrCand->GetClusterX() - Tracks[0].GetLKrClusterX();
    Double_t Y1 = LKrCand->GetClusterY() - Tracks[0].GetLKrClusterY();
    Double_t X2 = LKrCand->GetClusterX() - Tracks[1].GetLKrClusterX();
    Double_t Y2 = LKrCand->GetClusterY() - Tracks[1].GetLKrClusterY();
    Double_t LKrCandPos1 = sqrt(X1*X1 + Y1*Y1);
    Double_t LKrCandPos2 = sqrt(X2*X2 + Y2*Y2);
    Bool_t LKrPos = (fabs(LKrCandPos1) > 200. && fabs(LKrCandPos2) > 200.);
    Double_t LKrCandE = LKrCand->GetClusterEnergy();
    Double_t LKrCandDt1 = LKrCand->GetTime() - CHODTime1;
    Double_t LKrCandDt2 = LKrCand->GetTime() - CHODTime2;
    Bool_t LKrTime = (fabs(LKrCandDt1) < 5. && fabs(LKrCandDt2) < 5.);
    
    FillHisto("hSingleAddEnLKrCand", LKrCand->GetTime() - MUV3Time, LKrCandE/1000.);

    if (LKrPos && LKrCandE > 40. && LKrTime) {
      LKrEnergyInTime += LKrCandE;
      LKrWeightedTime += LKrCandE * LKrCand->GetTime();
    }
  }
  
  if (LKrEnergyInTime > 1000.) {
    LKrWeightedTime /= LKrEnergyInTime;
    FillHisto("hAddEnLKrCand", LKrWeightedTime - MUV3Time, LKrEnergyInTime/1000.);
    return;
  }
  
  LKrEnergyInTime = 0;
  LKrWeightedTime = 0;
  
  // Cut on total energy of hits in time with CHOD, more than 20 cm distant from extrapolation points of both tracks, single hit energy > 40 MeV
  
  TClonesArray& LKrHits = (*(LKrEvent->GetHits()));
  for (Int_t i = 0; i < NLKrHits; i++) {
    TRecoLKrHit* LKrHit = (TRecoLKrHit*) LKrHits[i];
    Double_t X1 = LKrHit->GetPosition().x() - Tracks[0].GetLKrClusterX();
    Double_t Y1 = LKrHit->GetPosition().y() - Tracks[0].GetLKrClusterY();
    Double_t X2 = LKrHit->GetPosition().x() - Tracks[1].GetLKrClusterX();
    Double_t Y2 = LKrHit->GetPosition().y() - Tracks[1].GetLKrClusterY();
    Double_t LKrHitPos1 = sqrt(X1*X1 + Y1*Y1);
    Double_t LKrHitPos2 = sqrt(X2*X2 + Y2*Y2);
    Bool_t LKrPos = (fabs(LKrHitPos1) > 200. && fabs(LKrHitPos2) > 200.);
    Double_t LKrHitE = LKrHit->GetEnergy();
    Double_t LKrHitTime = LKrHit->GetTime();
    Double_t LKrTime = (fabs(LKrHitTime - CHODTime1) < 5. && fabs(LKrHitTime - CHODTime2) < 5.);
    
    FillHisto("hSingleAddEnLKrHit", LKrHitTime - MUV3Time, LKrHitE/1000.);
    
    if (LKrPos && LKrHitE > 40. && LKrTime) {
      LKrEnergyInTime += LKrHitE;
      LKrWeightedTime += LKrHitE * LKrHitTime;
    }
  }
  
  if (LKrEnergyInTime > 1000.) {
    LKrWeightedTime /= LKrEnergyInTime;    
    FillHisto("hAddEnLKrHit", LKrWeightedTime - MUV3Time, LKrEnergyInTime/1000.);
    return;
  }
  */
  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;
  
  FillHisto("hCDAvsZVertex_TotMomToBeamlineAfterVetoes",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 4
  FillHisto("hZvertexvsBeamlineDistAfterVetoes",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 4
  FillHisto("hBeamlineDistvsTargetDist_TotMomAfterVetoes", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 4

  // Geometrical cuts, CUT 11: Cut on CDA of two tracks

  if (CDA >= 10.)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;
  
  // Geometrical cuts, CUT 12: Cut on two-track vertex wrt beamline

  if (BeamlineDist <= 100.)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;
  
  // Geometrical cuts, CUT 13: Cut on Z of two-track vertex

  if (Zvertex <= 102500. || Zvertex >= 180000.)
    return;

  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  // Geometrical cuts, CUT 14: Cut on CDA of each track wrt beam axis
  
  if (CDA1 <= 50. || CDA2 <= 50.)
    return;
  
  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;
  
  FillHisto("hCDAvsZVertex_TrackToBeamlineAfterCut",           Zvertex1 / 1000.,         CDA1 / 1000.);
  FillHisto("hCDAvsZVertex_TrackToBeamlineAfterCut",           Zvertex2 / 1000.,         CDA2 / 1000.);
  FillHisto("hCDAvsZVertex_TrackToTrackAfterCut",               Zvertex / 1000.,          CDA / 1000.);
  FillHisto("hCDAvsZVertex_TotMomToBeamlineAfterGeomCuts",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 5
  FillHisto("hZvertexvsBeamlineDistAfterGeomCuts",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 5
  FillHisto("hBeamlineDistvsTargetDist_TotMomAfterGeomCuts", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 5
  
  // Plot CDA vs Zvertex for events surviving previous selections
  
  FillHisto("hCDAvsZVertex_TrackToBeamlineFinal",      Zvertex1 / 1000.,         CDA1 / 1000.);
  FillHisto("hCDAvsZVertex_TrackToBeamlineFinal",      Zvertex2 / 1000.,         CDA2 / 1000.);
  FillHisto("hCDAvsZVertex_TrackToTrackFinal",          Zvertex / 1000.,          CDA / 1000.);
  FillHisto("hCDAvsZVertex_TotMomToBeamlineFinal",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 6
  FillHisto("hZvertexvsBeamlineDistFinal",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 6
  FillHisto("hBeamlineDistvsTargetDist_TotMomFinal", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 6

  // Computation of invariant mass
  
  TVector3 threeMomPi;
  TVector3 threeMomMu;
  Double_t energyMu;
  Double_t energyPi;
  Double_t invMass = 0.;

  if (Assoc == 1) {
    threeMomPi = SpectrometerCand2->GetThreeMomentumBeforeMagnet();
    threeMomMu = SpectrometerCand1->GetThreeMomentumBeforeMagnet();
  }
  else if (Assoc == 2) {
    threeMomPi = SpectrometerCand1->GetThreeMomentumBeforeMagnet();
    threeMomMu = SpectrometerCand2->GetThreeMomentumBeforeMagnet();
  }

  energyPi = TMath::Sqrt(threeMomPi.Px()*threeMomPi.Px() + threeMomPi.Py()*threeMomPi.Py() + threeMomPi.Pz()*threeMomPi.Pz() + fMpi*fMpi);
  energyMu = TMath::Sqrt(threeMomMu.Px()*threeMomMu.Px() + threeMomMu.Py()*threeMomMu.Py() + threeMomMu.Pz()*threeMomMu.Pz() + fMmu*fMmu);
  invMass = TMath::Sqrt((energyPi + energyMu)*(energyPi + energyMu) - (threeMomPi + threeMomMu).Mag2());

  FillHisto("hInvMassReco", invMass);

  // Computation of coupling-related quantities of the only good HNL in each event

  if (GetWithMC()) {
    Event *evt = GetMCEvent();
    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);
      if (p->GetParentID() == -1 && p->GetPDGcode() == 999 && p->GetEndProcessName() == "good") {
	point1.SetXYZ(p->GetProdPos().X(), p->GetProdPos().Y(), p->GetProdPos().Z());
	point2.SetXYZ(0., 0., fLInitialFV);
	momentum1.SetXYZ(p->GetInitial4Momentum().Px(), p->GetInitial4Momentum().Py(), p->GetInitial4Momentum().Pz());
	MN = ComputeHNLMass(p);
	gammaTot = GammaTot(MN);
	HNLTau = tauN(MN);
	LReach = ComputeL(point1, point2, momentum1);
	ProdFactor = ComputeProd(p, MN);
	DecayFactor = ComputeDecay(MN);
	NReachProb = ComputeNReachProb(p, HNLTau, LReach);
	NDecayProb = ComputeNDecayProb(p, HNLTau, fLFV);

	if (p->GetParticleName().Contains("e") && !p->GetParticleName().Contains("nu_tau"))
	  LeptonUSquared = fUeSquared;
	else if (p->GetParticleName().Contains("mu") && !p->GetParticleName().Contains("nu_tau"))
	  LeptonUSquared = fUmuSquared;
	else if (p->GetParticleName() == "DS->Ntau" || p->GetParticleName().Contains("nu_tau") || (p->GetParticleName().Contains("tau") && (p->GetParticleName().Contains("rho") || p->GetParticleName().Contains("pi") || p->GetParticleName().Contains("K"))))
	  LeptonUSquared = fUtauSquared;

        if (p->GetProdPos().Z() < fTAXDistance)
	  DProdProb = fDBeProdProb;
	else if(p->GetProdPos().Z() >= fTAXDistance)
          DProdProb = fDCuProdProb;
	
        Weight = DProdProb*fDDecayProb*NReachProb*NDecayProb*DecayFactor*ProdFactor*LeptonUSquared;
	fSumGood += Weight;
      }
    }
  }
}

void HeavyNeutrino::EndOfBurstUser() {
  
  FillHisto("hNbursts", 0.5);
}

void HeavyNeutrino::EndOfJobUser() {
  
  // Retrieve histos

  fhNk3pi    = (TH1D*) fHisto.GetTH1("hNk3pi");
  fhNbursts  = (TH1D*) fHisto.GetTH1("hNbursts");
  fhNEvents  = (TH1D*) fHisto.GetTH1("hNEvents");
  fhN2tracks = (TH1D*) fHisto.GetTH1("hN2tracks");
  fhNtracks  = (TH1D*) fHisto.GetTH1("hNtracks");

  fhZDProd       = (TH1D*) fHisto.GetTH1("hZDProd");
  fhZDDecay      = (TH1D*) fHisto.GetTH1("hZDDecay");
  fhDTheta       = (TH1D*) fHisto.GetTH1("hDTheta");
  fhDLambda      = (TH1D*) fHisto.GetTH1("hDLambda");
  fhDPath        = (TH1D*) fHisto.GetTH1("hDPath");
  fhDMom         = (TH1D*) fHisto.GetTH1("hDMom");
  fhZHNLDecay    = (TH1D*) fHisto.GetTH1("hZHNLDecay");
  fhHNLGamma     = (TH1D*) fHisto.GetTH1("hHNLGamma");
  fhHNLDecayProb = (TH1D*) fHisto.GetTH1("hHNLDecayProb");
  fhHNLReachProb = (TH1D*) fHisto.GetTH1("hHNLReachProb");
  fhHNLTheta     = (TH1D*) fHisto.GetTH1("hHNLTheta");
  fhHNLMom       = (TH1D*) fHisto.GetTH1("hHNLMom");
  fhMomPi        = (TH1D*) fHisto.GetTH1("hMomPi");
  fhMomMu        = (TH1D*) fHisto.GetTH1("hMomMu");
  fhWeight       = (TH1D*) fHisto.GetTH1("hWeight");

  fhXYSpec0Reco = (TH2D*) fHisto.GetTH2("hXYSpec0Reco");
  fhXYSpec1Reco = (TH2D*) fHisto.GetTH2("hXYSpec1Reco");
  fhXYSpec2Reco = (TH2D*) fHisto.GetTH2("hXYSpec2Reco");
  fhXYSpec3Reco = (TH2D*) fHisto.GetTH2("hXYSpec3Reco");
  fhXYCHODReco  = (TH2D*) fHisto.GetTH2("hXYCHODReco");
  fhXYCHODTrue  = (TH2D*) fHisto.GetTH2("hXYCHODTrue");
  fhXYMUV3True  = (TH2D*) fHisto.GetTH2("hXYMUV3True");
  fhP1vsP2      = (TH2D*) fHisto.GetTH2("hP1vsP2");

  fhPhysicsEventsVsCuts = (TH1D*) fHisto.GetTH1("hPhysicsEventsVsCuts");
  
  fhCDAvsZVertex_TotMomToBeamlineInitial              = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TotMomToBeamlineInitial");
  fhCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack");
  fhCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts      = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts");
  fhCDAvsZVertex_TotMomToBeamlineAfterVetoes          = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TotMomToBeamlineAfterVetoes");
  fhCDAvsZVertex_TotMomToBeamlineAfterGeomCuts        = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TotMomToBeamlineAfterGeomCuts");
  fhCDAvsZVertex_TotMomToBeamlineFinal                = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TotMomToBeamlineFinal");
  
  fhZvertexvsBeamlineDistInitial              = (TH2D*) fHisto.GetTH2("hZvertexvsBeamlineDistInitial");
  fhZvertexvsBeamlineDistAfterDownstreamTrack = (TH2D*) fHisto.GetTH2("hZvertexvsBeamlineDistAfterDownstreamTrack");
  fhZvertexvsBeamlineDistAfterEnergyCuts      = (TH2D*) fHisto.GetTH2("hZvertexvsBeamlineDistAfterEnergyCuts");
  fhZvertexvsBeamlineDistAfterVetoes          = (TH2D*) fHisto.GetTH2("hZvertexvsBeamlineDistAfterVetoes");
  fhZvertexvsBeamlineDistAfterGeomCuts        = (TH2D*) fHisto.GetTH2("hZvertexvsBeamlineDistAfterGeomCuts");
  fhZvertexvsBeamlineDistFinal                = (TH2D*) fHisto.GetTH2("hZvertexvsBeamlineDistFinal");
  
  fhCDAvsZVertex_TrackToBeamlineInitial  = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TrackToBeamlineInitial");
  fhCDAvsZVertex_TrackToTrackInitial     = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TrackToTrackInitial");
  fhCDAvsZVertex_TrackToBeamlineAfterCut = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TrackToBeamlineAfterCut");
  fhCDAvsZVertex_TrackToTrackAfterCut    = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TrackToTrackAfterCut");
  fhCDAvsZVertex_TrackToBeamlineFinal    = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TrackToBeamlineFinal");
  fhCDAvsZVertex_TrackToTrackFinal       = (TH2D*) fHisto.GetTH2("hCDAvsZVertex_TrackToTrackFinal");
  
  fhBeamlineDistvsTargetDist_TotMomInitial              = (TH2D*) fHisto.GetTH2("hBeamlineDistvsTargetDist_TotMomInitial");
  fhBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack = (TH2D*) fHisto.GetTH2("hBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack");
  fhBeamlineDistvsTargetDist_TotMomAfterEnergyCuts      = (TH2D*) fHisto.GetTH2("hBeamlineDistvsTargetDist_TotMomAfterEnergyCuts");
  fhBeamlineDistvsTargetDist_TotMomAfterVetoes          = (TH2D*) fHisto.GetTH2("hBeamlineDistvsTargetDist_TotMomAfterVetoes");
  fhBeamlineDistvsTargetDist_TotMomAfterGeomCuts        = (TH2D*) fHisto.GetTH2("hBeamlineDistvsTargetDist_TotMomAfterGeomCuts");
  fhBeamlineDistvsTargetDist_TotMomFinal                = (TH2D*) fHisto.GetTH2("hBeamlineDistvsTargetDist_TotMomFinal");
  
  fhDeltaTimeFromCHOD     = (TH1D*) fHisto.GetTH1("hDeltaTimeFromCHOD");
  fhNMUV3CandAssocToTrack = (TH1D*) fHisto.GetTH1("hNMUV3CandAssocToTrack");
  fhNCHODCandAssocToTrack = (TH1D*) fHisto.GetTH1("hNCHODCandAssocToTrack");
  
  fhEoP       = (TH1D*) fHisto.GetTH1("hEoP");
  fhEoPMuVsPi = (TH2D*) fHisto.GetTH1("hEoPMuVsPi");
  
  fhSingleAddEnLKrHit  = (TH2D*) fHisto.GetTH2("hSingleAddEnLKrHit");
  fhSingleAddEnLKrCand = (TH2D*) fHisto.GetTH2("hSingleAddEnLKrCand");
  fhAddEnLKrHit        = (TH2D*) fHisto.GetTH2("hAddEnLKrHit");
  fhAddEnLKrCand       = (TH2D*) fHisto.GetTH2("hAddEnLKrCand");
  
  fhInvMassMC    = (TH1D*) fHisto.GetTH1("hInvMassMC");
  fhInvMassReco  = (TH1D*) fHisto.GetTH1("hInvMassReco");

  fhAcc   = (TH1D*) fHisto.GetTH1("hAcc");
  fhYield = (TH1D*) fHisto.GetTH1("hYield");

  // X axis title

  fhNk3pi   ->GetXaxis()->SetTitle("Number of k3pi");
  fhNbursts ->GetXaxis()->SetTitle("Number of bursts");
  fhNEvents ->GetXaxis()->SetTitle("Number of events");
  fhN2tracks->GetXaxis()->SetTitle("Number of two-track events");
  fhNtracks ->GetXaxis()->SetTitle("Number of tracks in each event");

  fhZDProd      ->GetXaxis()->SetTitle("Z [mm]");
  fhZDDecay     ->GetXaxis()->SetTitle("Z [mm]");
  fhDTheta      ->GetXaxis()->SetTitle("Theta [rad]");
  fhDLambda     ->GetXaxis()->SetTitle("Decay length [mm]");
  fhDPath       ->GetXaxis()->SetTitle("Z [mm]");
  fhDMom        ->GetXaxis()->SetTitle("P [GeV]");
  fhZHNLDecay   ->GetXaxis()->SetTitle("Z [m]");
  fhHNLGamma    ->GetXaxis()->SetTitle("Lorentz gamma");
  fhHNLDecayProb->GetXaxis()->SetTitle("Decay probability");
  fhHNLReachProb->GetXaxis()->SetTitle("Reach probability");
  fhHNLTheta    ->GetXaxis()->SetTitle("Theta [rad]");
  fhHNLMom      ->GetXaxis()->SetTitle("P [GeV]");
  fhMomPi       ->GetXaxis()->SetTitle("P [GeV]");
  fhMomMu       ->GetXaxis()->SetTitle("P [GeV]");
  fhWeight      ->GetXaxis()->SetTitle("Weight");

  fhXYSpec0Reco->GetXaxis()->SetTitle("X [m]");
  fhXYSpec1Reco->GetXaxis()->SetTitle("X [m]");
  fhXYSpec2Reco->GetXaxis()->SetTitle("X [m]");
  fhXYSpec3Reco->GetXaxis()->SetTitle("X [m]");  
  fhXYCHODReco ->GetXaxis()->SetTitle("X [m]");
  fhXYCHODTrue ->GetXaxis()->SetTitle("X [m]");
  fhXYMUV3True ->GetXaxis()->SetTitle("X [m]");
  fhP1vsP2     ->GetXaxis()->SetTitle("P1 [GeV]");

  fhPhysicsEventsVsCuts->GetXaxis()->SetTitle("Cut ID");
    
  fhCDAvsZVertex_TotMomToBeamlineInitial             ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts     ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZVertex_TotMomToBeamlineAfterGeomCuts       ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZVertex_TotMomToBeamlineAfterVetoes         ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZVertex_TotMomToBeamlineFinal               ->GetXaxis()->SetTitle("Z [m]");
  
  fhZvertexvsBeamlineDistInitial             ->GetXaxis()->SetTitle("Z [m]");
  fhZvertexvsBeamlineDistAfterDownstreamTrack->GetXaxis()->SetTitle("Z [m]");
  fhZvertexvsBeamlineDistAfterEnergyCuts     ->GetXaxis()->SetTitle("Z [m]");
  fhZvertexvsBeamlineDistAfterGeomCuts       ->GetXaxis()->SetTitle("Z [m]");
  fhZvertexvsBeamlineDistAfterVetoes         ->GetXaxis()->SetTitle("Z [m]");
  fhZvertexvsBeamlineDistFinal               ->GetXaxis()->SetTitle("Z [m]");
  
  fhCDAvsZVertex_TrackToBeamlineInitial ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZVertex_TrackToBeamlineAfterCut->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZVertex_TrackToBeamlineFinal   ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZVertex_TrackToTrackInitial    ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZVertex_TrackToTrackAfterCut   ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZVertex_TrackToTrackFinal      ->GetXaxis()->SetTitle("Z [m]");
  
  fhBeamlineDistvsTargetDist_TotMomInitial             ->GetXaxis()->SetTitle("Target distance [m]");
  fhBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack->GetXaxis()->SetTitle("Target distance [m]");
  fhBeamlineDistvsTargetDist_TotMomAfterEnergyCuts     ->GetXaxis()->SetTitle("Target distance [m]");
  fhBeamlineDistvsTargetDist_TotMomAfterGeomCuts       ->GetXaxis()->SetTitle("Target distance [m]");
  fhBeamlineDistvsTargetDist_TotMomAfterVetoes         ->GetXaxis()->SetTitle("Target distance [m]");
  fhBeamlineDistvsTargetDist_TotMomFinal               ->GetXaxis()->SetTitle("Target distance [m]");
  
  fhDeltaTimeFromCHOD    ->GetXaxis()->SetTitle("Time difference [ns]");
  fhNMUV3CandAssocToTrack->GetXaxis()->SetTitle("Number of candidates");
  fhNCHODCandAssocToTrack->GetXaxis()->SetTitle("Number of candidates");
  
  fhEoP      ->GetXaxis()->SetTitle("E/p");
  fhEoPMuVsPi->GetXaxis()->SetTitle("Pion E/p");
  
  fhSingleAddEnLKrHit ->GetXaxis()->SetTitle("LKr time-track time [ns]");
  fhSingleAddEnLKrCand->GetXaxis()->SetTitle("LKr time-track time [ns]");
  fhAddEnLKrHit       ->GetXaxis()->SetTitle("LKr time-track time [ns]");
  fhAddEnLKrCand      ->GetXaxis()->SetTitle("LKr time-track time [ns]");
  
  fhInvMassMC   ->GetXaxis()->SetTitle("Invariant mass [MeV]");
  fhInvMassReco ->GetXaxis()->SetTitle("Invariant mass [MeV]");

  fhAcc  ->GetXaxis()->SetTitle("Acceptance");
  fhYield->GetXaxis()->SetTitle("Yield");

  // Y axis title

  fhXYSpec0Reco->GetYaxis()->SetTitle("Y [m]");
  fhXYSpec1Reco->GetYaxis()->SetTitle("Y [m]");
  fhXYSpec2Reco->GetYaxis()->SetTitle("Y [m]");
  fhXYSpec3Reco->GetYaxis()->SetTitle("Y [m]");  
  fhXYCHODReco ->GetYaxis()->SetTitle("Y [m]");
  fhXYCHODTrue ->GetYaxis()->SetTitle("Y [m]");
  fhXYMUV3True ->GetYaxis()->SetTitle("Y [m]");
  fhP1vsP2     ->GetYaxis()->SetTitle("P2 [GeV]");

  fhCDAvsZVertex_TotMomToBeamlineInitial             ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts     ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZVertex_TotMomToBeamlineAfterGeomCuts       ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZVertex_TotMomToBeamlineAfterVetoes         ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZVertex_TotMomToBeamlineFinal               ->GetYaxis()->SetTitle("CDA [m]");
  
  fhZvertexvsBeamlineDistInitial             ->GetYaxis()->SetTitle("Distance [m]");
  fhZvertexvsBeamlineDistAfterDownstreamTrack->GetYaxis()->SetTitle("Distance [m]");
  fhZvertexvsBeamlineDistAfterEnergyCuts     ->GetYaxis()->SetTitle("Distance [m]");
  fhZvertexvsBeamlineDistAfterGeomCuts       ->GetYaxis()->SetTitle("Distance [m]");
  fhZvertexvsBeamlineDistAfterVetoes         ->GetYaxis()->SetTitle("Distance [m]");
  fhZvertexvsBeamlineDistFinal               ->GetYaxis()->SetTitle("Distance [m]");
  
  fhCDAvsZVertex_TrackToBeamlineInitial ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZVertex_TrackToBeamlineAfterCut->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZVertex_TrackToBeamlineFinal   ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZVertex_TrackToTrackInitial    ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZVertex_TrackToTrackAfterCut   ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZVertex_TrackToTrackFinal      ->GetYaxis()->SetTitle("CDA [m]");
  
  fhBeamlineDistvsTargetDist_TotMomInitial             ->GetYaxis()->SetTitle("Beam axis distance [m]");
  fhBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack->GetYaxis()->SetTitle("Beam axis distance [m]");
  fhBeamlineDistvsTargetDist_TotMomAfterEnergyCuts     ->GetYaxis()->SetTitle("Beam axis distance [m]");
  fhBeamlineDistvsTargetDist_TotMomAfterGeomCuts       ->GetYaxis()->SetTitle("Beam axis distance [m]");
  fhBeamlineDistvsTargetDist_TotMomAfterVetoes         ->GetYaxis()->SetTitle("Beam axis distance [m]");
  fhBeamlineDistvsTargetDist_TotMomFinal               ->GetYaxis()->SetTitle("Beam axis distance [m]");
  
  fhEoPMuVsPi->GetYaxis()->SetTitle("Muon E/p");
  
  fhSingleAddEnLKrHit ->GetYaxis()->SetTitle("LKr additional energy [MeV]");
  fhSingleAddEnLKrCand->GetYaxis()->SetTitle("LKr additional energy [MeV]");
  fhAddEnLKrHit       ->GetYaxis()->SetTitle("LKr additional energy [GeV]");
  fhAddEnLKrCand      ->GetYaxis()->SetTitle("LKr additional energy [GeV]");

  // Acceptance computation

  Double_t Acceptance = fSumGood/fSumAll; 
  Double_t MeanProb = fSumAll/fNevents;
  Double_t Yield = Acceptance*MeanProb;

  cout << endl << "Acceptance = " << fSumGood << "/" << fSumAll << " = " << Acceptance << endl;
  cout << "Mean probability = " << fSumAll << "/" << fNevents << " = " << MeanProb << endl;
  cout << "Yield = " << Acceptance << "*" << MeanProb << " = " << Yield << endl;

  FillHisto("hAcc", Acceptance);
  FillHisto("hYield", Yield);

  // Plot residual number of events after each cut

  const int NCuts = 25;
  const char *CutNames[NCuts]  = {"Total", "TriggerOK", "2 tracks", "Straw0 acc", "Straw1 acc", "Straw2 acc", "Straw3 acc", "Chi2", "Straw chambers", "Charge", "CHOD acc", "CHOD assoc", "LKr acc", "LKr assoc", "MUV3 acc", "MUV3 assoc", "Mu E/p", "Pi E/p", "LAV veto", "SAV veto", "LKr veto", "CDA tracks", "Beam distance", "Z vertex", "CDA beam"};

  for (Int_t i = 1; i <= NCuts; i++)
    fhPhysicsEventsVsCuts->GetXaxis()->SetBinLabel(i, CutNames[i-1]);
  
  fhPhysicsEventsVsCuts->GetXaxis()->LabelsOption("v");

  SaveAllPlots();

  return;
}

HeavyNeutrino::~HeavyNeutrino() {

  fhNk3pi    = nullptr;
  fhNbursts  = nullptr;
  fhNEvents  = nullptr;
  fhN2tracks = nullptr;
  fhNtracks  = nullptr;

  fhZDProd       = nullptr;
  fhZDDecay      = nullptr;
  fhDTheta       = nullptr;
  fhDLambda      = nullptr;
  fhDPath        = nullptr;
  fhDMom         = nullptr;
  fhZHNLDecay    = nullptr;
  fhHNLGamma     = nullptr;
  fhHNLDecayProb = nullptr;
  fhHNLReachProb = nullptr;
  fhHNLTheta     = nullptr;
  fhHNLMom       = nullptr;
  fhWeight       = nullptr;
  fhMomPi        = nullptr;
  fhMomMu        = nullptr;

  fhXYSpec0Reco = nullptr;
  fhXYSpec1Reco = nullptr;
  fhXYSpec2Reco = nullptr;
  fhXYSpec3Reco = nullptr;
  fhXYCHODReco  = nullptr;
  fhXYCHODTrue  = nullptr;
  fhXYMUV3True  = nullptr;
  fhP1vsP2      = nullptr;

  fhPhysicsEventsVsCuts = nullptr;

  fhCDAvsZVertex_TotMomToBeamlineInitial              = nullptr;
  fhCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack = nullptr;
  fhCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts      = nullptr;
  fhCDAvsZVertex_TotMomToBeamlineAfterGeomCuts        = nullptr;
  fhCDAvsZVertex_TotMomToBeamlineAfterVetoes          = nullptr;
  fhCDAvsZVertex_TotMomToBeamlineFinal                = nullptr;

  fhZvertexvsBeamlineDistInitial              = nullptr;
  fhZvertexvsBeamlineDistAfterDownstreamTrack = nullptr;
  fhZvertexvsBeamlineDistAfterEnergyCuts      = nullptr;
  fhZvertexvsBeamlineDistAfterGeomCuts        = nullptr;
  fhZvertexvsBeamlineDistAfterVetoes          = nullptr;
  fhZvertexvsBeamlineDistFinal                = nullptr;

  fhCDAvsZVertex_TrackToBeamlineInitial  = nullptr;
  fhCDAvsZVertex_TrackToTrackInitial     = nullptr;
  fhCDAvsZVertex_TrackToBeamlineAfterCut = nullptr;
  fhCDAvsZVertex_TrackToTrackAfterCut    = nullptr;
  fhCDAvsZVertex_TrackToBeamlineFinal    = nullptr;
  fhCDAvsZVertex_TrackToTrackFinal       = nullptr;

  fhBeamlineDistvsTargetDist_TotMomInitial              = nullptr;
  fhBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack = nullptr;
  fhBeamlineDistvsTargetDist_TotMomAfterEnergyCuts      = nullptr;
  fhBeamlineDistvsTargetDist_TotMomAfterGeomCuts        = nullptr;
  fhBeamlineDistvsTargetDist_TotMomAfterVetoes          = nullptr;
  fhBeamlineDistvsTargetDist_TotMomFinal                = nullptr;

  fhDeltaTimeFromCHOD     = nullptr;
  fhNMUV3CandAssocToTrack = nullptr;
  fhNCHODCandAssocToTrack = nullptr;

  fhEoP       = nullptr;
  fhEoPMuVsPi = nullptr;

  fhSingleAddEnLKrHit  = nullptr;
  fhSingleAddEnLKrCand = nullptr;
  fhAddEnLKrHit        = nullptr;
  fhAddEnLKrCand       = nullptr;

  fhInvMassMC   = nullptr;
  fhInvMassReco = nullptr;

  fhAcc   = nullptr;
  fhYield = nullptr;    
}
