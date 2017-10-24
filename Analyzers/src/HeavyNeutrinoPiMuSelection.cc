#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include "HeavyNeutrinoPiMuSelection.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "GeometricAcceptance.hh"
#include "TriggerConditions.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

#define TRIGGER_L0_PHYSICS_TYPE 1
#define MCTriggerMask 0xFF
#define LabelSize 0.05

/// \class HeavyNeutrinoPiMuSelection

HeavyNeutrinoPiMuSelection::HeavyNeutrinoPiMuSelection(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrinoPiMuSelection") {
  
  RequestAllMCTrees();
  RequestAllRecoTrees();
  RequestL0Data();
  
  fCDAcomp     = new TwoLinesCDA();
  fDistcomp    = new PointLineDistance();
  fLAVMatching = new LAVMatching();
  fSAVMatching = new SAVMatching();
  
  AddParam("SumGood", &fSumGood, 0.);
  AddParam("Nevents", &fNevents, 0.);
  AddParam("SumAll", &fSumAll, 0.);
  AddParam("DBeProdProb", &fDBeProdProb, 0.);
  AddParam("DCuProdProb", &fDCuProdProb, 0.);
  AddParam("DDecayProb", &fDDecayProb, 0.);
  AddParam("NPiMu", &fNPiMu, 0.);
}

void HeavyNeutrinoPiMuSelection::InitHist() {

  BookHisto("hNk3pi",    new TH1D("Nk3pi", "Total number of K3pi events", 1, 0., 1.));
  BookHisto("hNbursts",  new TH1D("Nbursts", "Total number of processed bursts", 1, 0., 1.));
  BookHisto("hNEvents",  new TH1D("NEvents", "Number of total processed events" , 1, 0., 1.));
  BookHisto("hNtracks",  new TH1D("Ntracks",               "Number of tracks", 4, -0.5, 3.5));
  BookHisto("hN2tracks", new TH1D("N2tracks",      "Number of two-tracks events", 1, 0., 1.));

  BookHisto("hZDProd",        new TH1D("ZDProd",  "Z of D meson production point", 20000., -250., 33000.));
  BookHisto("hZDDecay",       new TH1D("ZDDecay",  "Z of D meson decay point",     20000., -250., 33000.));
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

  BookHisto("hXYSpec0Reco", new TH2D("XYSpec0Reco", "Two-track reconstructed events at CH1",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYSpec1Reco", new TH2D("XYSpec1Reco", "Two-track reconstructed events at CH2",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYSpec2Reco", new TH2D("XYSpec2Reco", "Two-track reconstructed events at CH3",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYSpec3Reco", new TH2D("XYSpec3Reco", "Two-track reconstructed events at CH4",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYCHODReco",  new TH2D("XYCHODReco",  "Two-track reconstructed events at CHOD", 100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYCHODTrue",  new TH2D("XYCHODTrue",  "X,Y of HNL daughters at CHOD, from MC", 100, -2., 2., 100, -2., 2.));
  BookHisto("hXYMUV3True",  new TH2D("XYMUV3True",  "X,Y of HNL daughters at MUV3, from MC", 100, -2., 2., 100, -2., 2.));
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
  
  BookHisto("hDeltaTimeFromCHOD",     new TH1D("DeltaTimeFromCHOD",     "Time difference of two tracks (CHOD candidates)",         200, -20., 20.));  
  BookHisto("hNMUV3CandAssocToTrack", new TH1D("NMUV3CandAssocToTrack", "Number of MUV3 candidates associated to each track", 4, -0.5, 3.5));
  BookHisto("hNCHODCandAssocToTrack", new TH1D("NCHODCandAssocToTrack", "Number of CHOD candidates associated to each track", 10, 0., 10.));
  
  BookHisto("hEoP", new TH1D("EoP", "E/p in LKr", 100, 0., 1.2));
  BookHisto("hEoPMuVsPi", new TH2D("EoPMuVsPi", "Muon E/p vs pion E/p in LKr", 100, 0., 1.2, 100, 0., 0.1));
  
  BookHisto("hSingleAddEnLKrHit",  new TH2D("SingleAddEnLKrHit", "Single hit energy vs time, for additional LKr hits",              250, -100., 100, 250, 0., 100.));
  BookHisto("hSingleAddEnLKrCand", new TH2D("SingleAddEnLKrCand", "Single candidate energy vs time, for additional LKr candidates", 250, -100., 100, 250, 0., 100.));
  BookHisto("hAddEnLKrHit",        new TH2D("AddEnLKrHit", "Additional energy vs time, for LKr hits",                               250, -100., 100, 250, 0., 10.));
  BookHisto("hAddEnLKrCand",       new TH2D("AddEnLKrCand", "Additional energy vs time, for LKr candidates",                        250, -100., 100, 250, 0., 10.));
  
  BookHisto("hInvMassMC",   new TH1D("InvMassMC",   "Invariant mass MC",                  100, 999., 1001.));
  BookHisto("hInvMassReco", new TH1D("InvMassReco", "Invariant mass Reco", 50, 960., 1040.));

  BookHisto("hAcc",   new TH1D("Acc",   "Acceptance", 50, 1.E-6, 5.E-5));
  BookHisto("hYield", new TH1D("Yield", "Yield",      50, 1.E-20, 1.E-18));
}

void HeavyNeutrinoPiMuSelection::Process(Int_t) {

  TRecoLKrEvent*          LKrEvent  = (TRecoLKrEvent*)          GetEvent("LKr");
  TRecoLAVEvent*          LAVEvent  = (TRecoLAVEvent*)          GetEvent("LAV");
  TRecoIRCEvent*          IRCEvent  = (TRecoIRCEvent*)          GetEvent("IRC");
  TRecoSACEvent*          SACEvent  = (TRecoSACEvent*)          GetEvent("SAC");

  // Counter for cuts

  Int_t CutID = 0;
  
  FillHisto("hPhysicsEventsVsCuts", CutID);
  CutID++;

  // Plots of kineParts from MC

  Double_t zCHOD = 239009.0;
  Double_t zMUV3 = 246800.0;
  TLorentzVector mom1;
  TLorentzVector mom2;
  Double_t mass1 = 0.;
  Double_t mass2 = 0.;
  Double_t MN = 0.;
  Double_t HNLTau = 0.;
  Double_t gammaTot = 0.;
  Double_t NDecayProb = 0.;
  Double_t NReachProb = 0.;
  Double_t LFV = 77500.;
  Double_t LReach = 0.;
  Double_t LInitialFV = 102500.;
  Double_t USquared = 1.E-6;
  Double_t UeSquared = USquared/20.8;
  Double_t UmuSquared = 16.*UeSquared;
  Double_t UtauSquared = 3.8*UeSquared;
  Double_t LeptonUSquared = 0.;
  Double_t BR = 0.;
  Double_t Weight = 0.;
  Int_t counter = 0;
  TVector3 point1;
  TVector3 point2;
  TVector3 momentum1;
  Double_t p1,p2;
  Double_t DProdProb = 0.;
  Double_t fTAXDistance = 24685.;

  fDBeProdProb = 0.00069;
  fDCuProdProb = fDBeProdProb*TMath::Power((29./4.),1./3.); // ACu/ABe
  fDDecayProb = 1.;
  fNPiMu = 0.;

  if (GetWithMC()) {
    Event *evt = GetMCEvent();
    counter = 0;
    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);
      if (p->GetParentID() == 0) {
	counter++;
	Double_t xCHOD = p->xAt(zCHOD);                      
	Double_t yCHOD = p->yAt(zCHOD); 
	Double_t xMUV3 = p->xAt(zMUV3);
	Double_t yMUV3 = p->yAt(zMUV3);
	FillHisto("hXYCHODTrue", xCHOD/1000., yCHOD/1000.);
	FillHisto("hXYMUV3True", xMUV3/1000., yMUV3/1000.);
	if (p->GetParticleName() == "pi+" || p->GetParticleName() == "pi-") {
	  mom1 = p->GetInitial4Momentum();
	  mass1 = TMath::Sqrt(mom1.E()*mom1.E() - mom1.P()*mom1.P());
	  p1 = TMath::Sqrt(mom1.Px()*mom1.Px()+mom1.Py()*mom1.Py()+mom1.Pz()*mom1.Pz());
	}
	else if (p->GetParticleName() == "mu+" || p->GetParticleName() == "mu-") {
	  mom2 = p->GetInitial4Momentum();
	  mass2 = TMath::Sqrt(mom2.E()*mom2.E() - mom2.P()*mom2.P());
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

    // Computation of coupling-related quantities of all HNLs in the MC event

    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);      
      if (p->GetParentID() == -1 && p->GetPDGcode() == 999) {
	fNevents++;
	point1.SetXYZ(p->GetProdPos().X(), p->GetProdPos().Y(), p->GetProdPos().Z());
	point2.SetXYZ(0., 0., LInitialFV);
	momentum1.SetXYZ(p->GetInitial4Momentum().Px(), p->GetInitial4Momentum().Py(), p->GetInitial4Momentum().Pz());
	MN = ComputeHNLMass(p);
	gammaTot = GammaTot(mass1, mass2, MN, USquared);
	HNLTau = Tau(gammaTot);
	LReach = ComputeL(point1, point2, momentum1);
	BR = ComputeBR(p, MN);
	NReachProb = ComputeNReachProb(p, HNLTau, LReach);
	NDecayProb = ComputeNDecayProb(p, HNLTau, LFV);
	fNPiMu = GammaPiMu(MN)*UmuSquared/(GammaPiE(MN)*UeSquared + GammaPiMu(MN)*UmuSquared);

	if (p->GetParticleName() == "HNLD+e+" || p->GetParticleName() == "HNLD-e-" || p->GetParticleName() == "HNLDS+e+" || p->GetParticleName() == "HNLDS-e-") 
	  LeptonUSquared = UeSquared;

	else if (p->GetParticleName() == "HNLD+mu+" || p->GetParticleName() == "HNLD-mu-" || p->GetParticleName() == "HNLDS+mu+" || p->GetParticleName() == "HNLDS-mu-") 
	  LeptonUSquared = UmuSquared;

	if (p->GetProdPos().Z() < fTAXDistance)
	  DProdProb = fDBeProdProb;
	else if (p->GetProdPos().Z() >= fTAXDistance)
	  DProdProb = fDCuProdProb;

	Weight = DProdProb*fDDecayProb*NReachProb*NDecayProb*fNPiMu*BR*LeptonUSquared;

	if (Weight >= 0. && Weight <= 10)
	  fSumAll += Weight;

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
  
  // Count K3Pi
  
  Bool_t k3pi = *(Bool_t*) GetOutput("K3piSelection.EventSelected");
  Int_t RunNumber = GetWithMC() ? 0 : GetEventHeader()->GetRunID();

  if (k3pi && 0x10)
    FillHisto("hNk3pi", 0.5);

  // L0 data
  
  L0TPData *L0TPData = GetL0Data();
  UChar_t L0DataType = GetWithMC() ? 0x11 : L0TPData->GetDataType();
  UInt_t L0TriggerFlags = GetWithMC() ? 0xFF : L0TPData->GetTriggerFlags();
  Bool_t PhysicsTriggerOK = (L0DataType & TRIGGER_L0_PHYSICS_TYPE);
  Bool_t TriggerFlagsOK = L0TriggerFlags & MCTriggerMask;
  Bool_t TriggerOK = PhysicsTriggerOK && TriggerFlagsOK;
      
  // Select physics triggers and L0 trigger conditions
  // IMPORTANTE!!!!!!!!!!!!!!!! aggiungi le trigger condition del 2017!!!!!!!!!!!!!!!!
  
  Int_t ID1 = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-Q2-MO1");
  Int_t ID2 = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-Q2-M1");
  Int_t ID3 = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-QX-MO1");
  Int_t ID4 = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-QX-M1");
  Bool_t On1 = TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0TPData, ID1);
  Bool_t On2 = TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0TPData, ID2);
  Bool_t On3 = TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0TPData, ID3);
  Bool_t On4 = TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0TPData, ID4);
  Bool_t TriggerStreamOK = (On1 || On2 || On3 || On4);
  
  if (!TriggerOK)
    return;
  if (!TriggerStreamOK && !GetWithMC())
    return;

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

  // Features of the tracks
    
  Int_t Charge1 = Tracks[0].GetCharge();
  Int_t Charge2 = Tracks[1].GetCharge();
  Double_t ChiSquare1 = Tracks[0].GetChi2();
  Double_t ChiSquare2 = Tracks[1].GetChi2();
  TRecoSpectrometerCandidate* SpectrometerCand1 = Tracks[0].GetSpectrometerCandidate();
  TRecoSpectrometerCandidate* SpectrometerCand2 = Tracks[1].GetSpectrometerCandidate();
  TVector3 Mom1 = SpectrometerCand1->GetThreeMomentumBeforeMagnet();
  TVector3 Mom2 = SpectrometerCand2->GetThreeMomentumBeforeMagnet();
  TVector3 TotMom = Mom1 + Mom2;
  
  // (X,Y) of reconstructed tracks for all Spectrometer chambers
  
  Double_t zStraw[4] = {183508.0, 194066.0, 204459.0, 218885.0};
  Double_t xStrawChamberCentre[4] = {101.2, 114.4, 92.4, 52.8};
  Double_t rMinStraw = 60.0; 
  Double_t rMaxStraw = 1010.0;
  Double_t zCHODPlane = 239009.0;
  Double_t rMinCHOD   = 120.0;
  Double_t rMaxCHOD   = 1110.0;
  
  for (UInt_t i = 0; i < Tracks.size(); i++) {
    TRecoSpectrometerCandidate* Cand = Tracks[i].GetSpectrometerCandidate();
    for (Int_t j = 0; j < 4; j++) {
      Double_t x = Cand->xAt(zStraw[j]);
      Double_t y = Cand->yAt(zStraw[j]);
      Double_t r = sqrt(x*x + y*y); 
      Double_t rShifted = sqrt(pow(x-xStrawChamberCentre[j],2) + y*y); 
      if (rShifted > rMinStraw && r < rMaxStraw) {
	TString name = Form("hXYSpec%dReco",j);
	FillHisto(name, x/1000., y/1000.);
      }
    }
    Double_t x = Cand->xAtAfterMagnet(zCHODPlane);
    Double_t y = Cand->yAtAfterMagnet(zCHODPlane);
    Double_t r = sqrt(x*x+y*y);
    Double_t r1 = rMinCHOD;
    Double_t r2 = rMaxCHOD;
    if (r > r1 && r < r2)
	FillHisto("hXYCHODReco", x/1000., y/1000.);
  }

  // Compute CDA of track1,2 wrt kaon axis, between each other and of two-track total momentum wrt beam axis
  
  fCDAcomp->SetLine1Point1(0.0, 0.0, 102000.0);     // kaon axis
  fCDAcomp->SetDir1(1.2E-3, 0., 1.);
  
  fCDAcomp->SetLine2Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);     // track1
  fCDAcomp->SetDir2(Mom1);
  fCDAcomp->ComputeVertexCDA();
  
  Double_t CDA1 = fCDAcomp->GetCDA();
  Double_t Zvertex1 = fCDAcomp->GetVertex().z();     // kaon axis-track1
  
  fCDAcomp->SetLine2Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand2->yAtBeforeMagnet(10.0), 10.0);     // track2
  fCDAcomp->SetDir2(Mom2);
  fCDAcomp->ComputeVertexCDA();
  
  Double_t CDA2 = fCDAcomp->GetCDA();
  Double_t Zvertex2 = fCDAcomp->GetVertex().z();     // kaon axis-track2
  
  fCDAcomp->SetLine1Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);     // track1
  fCDAcomp->SetDir1(Mom1);
  fCDAcomp->SetLine2Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand2->yAtBeforeMagnet(10.0), 10.0);     // track2
  fCDAcomp->SetDir2(Mom2);
  fCDAcomp->ComputeVertexCDA();
  
  Double_t CDA = fCDAcomp->GetCDA();
  TVector3 Vertex = fCDAcomp->GetVertex();
  Double_t Zvertex = fCDAcomp->GetVertex().z();     // track1-track2
  
  fCDAcomp->SetLine1Point1(0.0, 0.0, 102000.0);     // beam axis
  fCDAcomp->SetDir1(0., 0., 1.);
  
  fCDAcomp->SetLine2Point1(Vertex);     // total momentum
  fCDAcomp->SetDir2(TotMom);
  fCDAcomp->ComputeVertexCDA();
  
  Double_t CDAMom = fCDAcomp->GetCDA();
  Double_t ZvertexMom = fCDAcomp->GetVertex().z();     // beam axis-total momentum
  
  // Compute distance of two-track momentum wrt target (extrapolation at target)
  
  fDistcomp->SetLineDir(TotMom);     // total momentum
  fDistcomp->SetLinePoint1(Vertex);
  fDistcomp->SetPoint(0., 0., 0.);     // target
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
      Double_t x1 = SpectrometerCand1->xAt(zStraw[i]);
      Double_t y1 = SpectrometerCand1->yAt(zStraw[i]);
      Double_t r1 = sqrt(x1*x1 + y1*y1);
      Double_t rShifted1 = sqrt(pow(x1-xStrawChamberCentre[i],2) + y1*y1);
      Double_t x2 = SpectrometerCand2->xAt(zStraw[i]);
      Double_t y2 = SpectrometerCand2->yAt(zStraw[i]);
      Double_t r2 = sqrt(x2*x2 + y2*y2);
      Double_t rShifted2 = sqrt(pow(x2-xStrawChamberCentre[i],2) + y2*y2);
      Bool_t inAcc = false;

      if ((rShifted1 > rMinStraw && r1 < rMaxStraw) && (rShifted2 > rMinStraw && r2 < rMaxStraw))
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

  // Plot the number of candidates associated to each track, for MUV3, CHOD and NewCHOD
  
  FillHisto("hNMUV3CandAssocToTrack", Tracks[0].GetNMUV3AssociationRecords());
  FillHisto("hNMUV3CandAssocToTrack", Tracks[1].GetNMUV3AssociationRecords());
  FillHisto("hNCHODCandAssocToTrack", Tracks[0].GetNCHODAssociationRecords());
  FillHisto("hNCHODCandAssocToTrack", Tracks[1].GetNCHODAssociationRecords());
  
  // Downstream track selection, CUT 4: Extrapolation and association to CHOD

  Bool_t CHODAssoc = (Tracks[0].CHODAssociationExists() && Tracks[1].CHODAssociationExists());
  Double_t x1 = SpectrometerCand1->xAtAfterMagnet(zCHODPlane);
  Double_t y1 = SpectrometerCand1->yAtAfterMagnet(zCHODPlane);
  Double_t r1 = sqrt(x1*x1+y1*y1);
  Double_t x2 = SpectrometerCand2->xAtAfterMagnet(zCHODPlane);
  Double_t y2 = SpectrometerCand2->yAtAfterMagnet(zCHODPlane);
  Double_t r2 = sqrt(x2*x2+y2*y2);
  Bool_t inAcc = false;

  if ((r1 > rMinCHOD && r1 < rMaxCHOD) && (r2 > rMinCHOD && r2 < rMaxCHOD))
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
  Int_t Assoc = 0;
  
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
  
  Double_t EoP1 = Tracks[0].GetLKrEoP();
  Double_t EoP2 = Tracks[1].GetLKrEoP();
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
    fSAVMatching->SetIRCTimeCuts(99999, 99999);
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
  
  TVector3 threeMomPi;
  TVector3 threeMomMu;
  Double_t energyMu;
  Double_t energyPi;
  Double_t massPi = 139.57;
  Double_t massMu = 105.66;
  Double_t invMass = 0.;

  if (Assoc == 1) {
    threeMomPi = SpectrometerCand2->GetThreeMomentumBeforeMagnet();
    threeMomMu = SpectrometerCand1->GetThreeMomentumBeforeMagnet();
  }
  else if (Assoc == 2) {
    threeMomPi = SpectrometerCand1->GetThreeMomentumBeforeMagnet();
    threeMomMu = SpectrometerCand2->GetThreeMomentumBeforeMagnet();
  }

  energyPi = TMath::Sqrt(threeMomPi.Px()*threeMomPi.Px() + threeMomPi.Py()*threeMomPi.Py() + threeMomPi.Pz()*threeMomPi.Pz() + massPi*massPi);
  energyMu = TMath::Sqrt(threeMomMu.Px()*threeMomMu.Px() + threeMomMu.Py()*threeMomMu.Py() + threeMomMu.Pz()*threeMomMu.Pz() + massMu*massMu);
  invMass = TMath::Sqrt((energyPi + energyMu)*(energyPi + energyMu) - (threeMomPi + threeMomMu).Mag2());

  FillHisto("hInvMassReco", invMass);

  // Computation of coupling-related quantities of good HNL in the MC event         

  if (GetWithMC()) {
    Event *evt = GetMCEvent();
    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);
      if (p->GetParentID() == -1 && p->GetPDGcode() == 999 && p->GetEndProcessName() == "good") {
	point1.SetXYZ(p->GetProdPos().X(), p->GetProdPos().Y(), p->GetProdPos().Z());
	point2.SetXYZ(0., 0., LInitialFV);
	momentum1.SetXYZ(p->GetInitial4Momentum().Px(), p->GetInitial4Momentum().Py(), p->GetInitial4Momentum().Pz());
	MN = ComputeHNLMass(p);
	gammaTot = GammaTot(mass1, mass2, MN, USquared);
	HNLTau = Tau(gammaTot);
	LReach = ComputeL(point1, point2, momentum1);
	BR = ComputeBR(p, MN);
	NReachProb = ComputeNReachProb(p, HNLTau, LReach);
	NDecayProb = ComputeNDecayProb(p, HNLTau, LFV);
	fNPiMu = GammaPiMu(MN)*UmuSquared/(GammaPiE(MN)*UeSquared + GammaPiMu(MN)*UmuSquared);

	if (p->GetParticleName() == "HNLD+e+" || p->GetParticleName() == "HNLD-e-" || p->GetParticleName() == "HNLDS+e+" || p->GetParticleName() == "HNLDS-e-") 
	  LeptonUSquared = UeSquared;

	else if (p->GetParticleName() == "HNLD+mu+" || p->GetParticleName() == "HNLD-mu-" || p->GetParticleName() == "HNLDS+mu+" || p->GetParticleName() == "HNLDS-mu-") 
	  LeptonUSquared = UmuSquared;

        if (p->GetProdPos().Z() < fTAXDistance)
	  DProdProb = fDBeProdProb;
	else if(p->GetProdPos().Z() >= fTAXDistance)
          DProdProb = fDCuProdProb;
	
        Weight = DProdProb*fDDecayProb*NReachProb*NDecayProb*fNPiMu*BR*LeptonUSquared;
	
	if (Weight >= 0. && Weight <= 10)
	  fSumGood += Weight;
      }
    }
  }
}

void HeavyNeutrinoPiMuSelection::EndOfBurstUser() {
  
  FillHisto("hNbursts", 0.5);
}

void HeavyNeutrinoPiMuSelection::EndOfJobUser() {
  
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

  const int NCuts = 25;
  const char *CutNames[NCuts]  = {"Total", "TriggerOK", "2 tracks", "Straw0 acc", "Straw1 acc", "Straw2 acc", "Straw3 acc", "Chi2", "Straw chambers", "Charge", "CHOD acc", "CHOD assoc", "LKr acc", "LKr assoc", "MUV3 acc", "MUV3 assoc", "Mu E/p", "Pi E/p", "LAV veto", "SAV veto", "LKr veto", "CDA tracks", "Beam distance", "Z vertex", "CDA beam"};

  for (Int_t i = 1; i <= NCuts; i++)
    fhPhysicsEventsVsCuts->GetXaxis()->SetBinLabel(i, CutNames[i-1]);
  
  fhPhysicsEventsVsCuts->GetXaxis()->LabelsOption("v");

  SaveAllPlots();

  return;
}

Double_t HeavyNeutrinoPiMuSelection::ComputeHNLMass(KinePart* p) {

  Double_t MN = TMath::Sqrt(p->GetInitial4Momentum().E()*p->GetInitial4Momentum().E() - p->GetInitial4Momentum().P()*p->GetInitial4Momentum().P());

  return MN;
}

Double_t HeavyNeutrinoPiMuSelection::ComputeL(TVector3 p1, TVector3 p2, TVector3 mom1) {
 
  TVector3 r(p1.x() + mom1.Px()/mom1.Pz()*(p2.z()-p1.z()), p1.y() + mom1.Py()/mom1.Pz()*(p2.z()-p1.z()), p2.z());
  Double_t x = r.x()-p1.x();
  Double_t y = r.y()-p1.y();
  Double_t z = r.z()-p1.z();
  Double_t L = TMath::Sqrt(x*x + y*y + z*z);

  return L;
}

Double_t HeavyNeutrinoPiMuSelection::ComputeNDecayProb(KinePart* p, Double_t tau, Double_t l) {

  Double_t cLight = 299.9; // mm/ns
  Double_t NProb = 0.;
  NProb = 1. - TMath::Exp(-l/(p->GetInitial4Momentum().Beta()*p->GetInitial4Momentum().Gamma()*cLight*tau));

  return NProb;
}


Double_t HeavyNeutrinoPiMuSelection::ComputeNReachProb(KinePart* p, Double_t tau, Double_t l) {

  Double_t cLight = 299.9; // mm/ns
  Double_t NProb = 0.;
  NProb = TMath::Exp(-l/(p->GetInitial4Momentum().Beta()*p->GetInitial4Momentum().Gamma()*cLight*tau));

  return NProb;
}

Double_t HeavyNeutrinoPiMuSelection::GammaPiE(Double_t MN) {

  Double_t Mpi = 139.6;
  Double_t Me = 0.511;
  Double_t GF = 1.1E-11; // MeV^-2
  Double_t fPi = 130.41; // MeV
  Double_t gamma_pi_e = 0.;

  gamma_pi_e = GF*GF*fPi*fPi/(8*TMath::Pi()*MN*MN)*(Me*Me + MN*MN)*(MN*MN - Mpi*Mpi - Me*Me)*TMath::Power(MN*MN*MN*MN + Mpi*Mpi*Mpi*Mpi + Me*Me*Me*Me - 2.*MN*Me*MN*Me - 2.*MN*Mpi*MN*Mpi - 2.*Mpi*Me*Mpi*Me,0.5)/(2*MN);
  
  return gamma_pi_e;
}

Double_t HeavyNeutrinoPiMuSelection::GammaPiMu(Double_t MN) {

  Double_t Mpi = 139.6;
  Double_t Mmu = 105.;
  Double_t GF = 1.1E-11; 
  Double_t fPi = 130.41; 
  Double_t gamma_pi_mu = 0.;

  gamma_pi_mu = GF*GF*fPi*fPi/(8*TMath::Pi()*MN*MN)*(Mmu*Mmu + MN*MN)*(MN*MN - Mpi*Mpi - Mmu*Mmu)*TMath::Power(MN*MN*MN*MN + Mpi*Mpi*Mpi*Mpi + Mmu*Mmu*Mmu*Mmu - 2.*MN*Mmu*MN*Mmu - 2.*MN*Mpi*MN*Mpi - 2.*Mpi*Mmu*Mpi*Mmu,0.5)/(2*MN);
  
  return gamma_pi_mu;
}

Double_t HeavyNeutrinoPiMuSelection::GammaENuNu(Double_t MN) {

  Double_t GF = 1.1E-11;
  Double_t Me = 0.511;
  Double_t r  = 4*Me*Me/(MN*MN);
  Double_t a = TMath::Power(1-r,0.5)*(1./3.-7*r/6.-r*r/24.-r*r*r/16.);
  Double_t b = r*r*(1-r*r/16.)*TMath::ATanH(TMath::Power(1-r,0.5));
  Double_t f = a+b;
  Double_t gamma_e_nu_nu = 0.;

  gamma_e_nu_nu = GF*GF*TMath::Power(MN,5)*f/(64.*TMath::Power(TMath::Pi(),3));
  
  return gamma_e_nu_nu;
}

Double_t HeavyNeutrinoPiMuSelection::GammaTot(Double_t mass1, Double_t mass2, Double_t MN, Double_t USquared) {

  Double_t UeSquared = USquared/20.8;
  Double_t UmuSquared = 16.*UeSquared;
  Double_t Mpi = 139.6;
  Double_t Me = 0.511;
  Double_t Mmu = 105.;
  Double_t gammaTot = 0.;

  if ((mass1 < 1. || mass2 < 1.) && MN < Mpi + Me)       // N->eenu                               
    gammaTot = GammaENuNu(MN)*UeSquared;
  else if ((mass1 < 1. || mass2 < 1.) && MN < Mpi + Mmu)  // N->epi                               
    gammaTot = GammaENuNu(MN)*UeSquared + GammaPiE(MN)*UeSquared;
  else if (((mass1 > 1. && mass1 < 120.) || (mass2 > 1. && mass2 < 120.)) && MN >= Mpi + Mmu)  // N->pimu                                                                      
    gammaTot = GammaENuNu(MN)*UeSquared + GammaPiE(MN)*UeSquared + GammaPiMu(MN)*UmuSquared;

  // For now

  gammaTot = GammaPiE(MN)*UeSquared + GammaPiMu(MN)*UmuSquared;

  return gammaTot;
}

Double_t HeavyNeutrinoPiMuSelection::Tau(Double_t gammaTot) {

  Double_t hc = 200.E-12; // MeV*mm
  Double_t cLight = 299.9;
  Double_t tau = 0.;

  tau = hc/(gammaTot*cLight);

  return tau;
}

Double_t HeavyNeutrinoPiMuSelection::NIntoPiMuBR(Double_t USquared, Double_t MN, Double_t tau) {

  Double_t UeSquared = USquared/20.8;
  Double_t UmuSquared = 16.*UeSquared;
  Double_t hc = 200.E-12;
  Double_t cLight = 299.9;
  Double_t nIntoPiMuBR = 0.;

  nIntoPiMuBR = GammaPiMu(MN)*UmuSquared*tau*cLight/hc;

  return nIntoPiMuBR;
}

Double_t HeavyNeutrinoPiMuSelection::PhaseSpace(Double_t mesonMass, Double_t MN, Double_t leptonMass) {

  Double_t phaseSpace = 0.;

  phaseSpace = TMath::Power(mesonMass*mesonMass-MN*MN-leptonMass*leptonMass,2) - 4.*MN*MN*leptonMass*leptonMass;

  return phaseSpace;
}

Double_t HeavyNeutrinoPiMuSelection::PhaseSpaceFactor(Double_t mesonMass, Double_t MN, Double_t leptonMass, Double_t phaseSpace) {

  Double_t factor = 0.;

  if(phaseSpace > 0.)
    factor = (mesonMass*mesonMass*(MN*MN + leptonMass*leptonMass)-TMath::Power(MN*MN-leptonMass*leptonMass,2))*TMath::Power(phaseSpace,0.5)/(leptonMass*leptonMass*TMath::Power(mesonMass*mesonMass-leptonMass*leptonMass,2));

  return factor;
}

Double_t HeavyNeutrinoPiMuSelection::ComputeBR(KinePart* p, Double_t MN) {

  Double_t Me = 0.511;
  Double_t Mmu = 105.;
  Double_t DMass = 1870.;
  Double_t DSMass = 1969.;
  Double_t De2BR = 9.2E-9;
  Double_t Dm2BR = 3.74E-4;
  Double_t DSe2BR = 1.4E-7;
  Double_t DSm2BR = 5.56E-3;
  Double_t mesonMass = 0.;

  if (p->GetParticleName() == "HNLD+e+" || p->GetParticleName() == "HNLD-e-")
    mesonMass = DMass;
 
  else if (p->GetParticleName() == "HNLD+mu+" || p->GetParticleName() == "HNLD-mu-")
    mesonMass = DMass;

  else if(p->GetParticleName() == "HNLDS+e+" || p->GetParticleName() == "HNLDS-e-")
    mesonMass = DSMass;
  
  else if (p->GetParticleName() == "HNLDS+mu+" || p->GetParticleName() == "HNLDS-mu-") 
    mesonMass = DSMass;
   
  Double_t phaseSpaceE = PhaseSpace(mesonMass, MN, Me);
  Double_t phaseSpaceMu = PhaseSpace(mesonMass, MN, Mmu);
  Double_t FactorE = PhaseSpaceFactor(mesonMass, MN, Me, phaseSpaceE);
  Double_t FactorMu = PhaseSpaceFactor(mesonMass, MN, Mmu, phaseSpaceMu);
  Double_t brt = 0.;

  if (mesonMass == DMass)
    brt = De2BR*FactorE + Dm2BR*FactorMu;

  else if (mesonMass == DSMass)
    brt = DSe2BR*FactorE + DSm2BR*FactorMu;

  return brt;
}

Double_t HeavyNeutrinoPiMuSelection::ComputeTotalBR(KinePart* p, Double_t USquared, Double_t MN) {

  Double_t UeSquared = USquared/20.8;
  Double_t UmuSquared = 16.*UeSquared;
  Double_t Me = 0.511;
  Double_t Mmu = 105.;
  Double_t DMass = 1870.;
  Double_t DSMass = 1969.;
  Double_t De2BR = 9.2E-9;
  Double_t Dm2BR = 3.82E-4;
  Double_t DSe2BR = 1.4E-7;
  Double_t DSm2BR = 5.9E-3;
  Double_t mesonMass = 0.;
  Double_t totalBR = 0.;

  if (p->GetParticleName() == "HNLD+e+" || p->GetParticleName() == "HNLD-e-" || p->GetParticleName() == "HNLD+mu+" || p->GetParticleName() == "HNLD-mu-")
    mesonMass = DMass;

  else if(p->GetParticleName() == "HNLDS+e+" || p->GetParticleName() == "HNLDS-e-" || p->GetParticleName() == "HNLDS+mu+" || p->GetParticleName() == "HNLDS-mu-")
    mesonMass = DSMass;

  Double_t phaseSpaceE = PhaseSpace(mesonMass, MN, Me);
  Double_t phaseSpaceMu = PhaseSpace(mesonMass, MN, Mmu);
  Double_t FactorE = PhaseSpaceFactor(mesonMass, MN, Me, phaseSpaceE);
  Double_t FactorMu = PhaseSpaceFactor(mesonMass, MN, Mmu, phaseSpaceMu);

  if (mesonMass == DMass)
    totalBR = De2BR*FactorE*UeSquared + Dm2BR*FactorMu*UmuSquared;

  else if (mesonMass == DSMass)
    totalBR = DSe2BR*FactorE*UeSquared + DSm2BR*FactorMu*UmuSquared;

  return totalBR;
}

HeavyNeutrinoPiMuSelection::~HeavyNeutrinoPiMuSelection() {
  
  delete fCDAcomp;
  delete fDistcomp;
  delete fLAVMatching;
  delete fSAVMatching;
}
