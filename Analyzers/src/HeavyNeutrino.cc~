// ---------------------------------------------------------------                                  
//                                                                                            
// History:                                                                                  
//                                                                                               
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018                              
//                                                                                        
// ---------------------------------------------------------------                         

/// \class HeavyNeutrino              
/// \Brief                                                                         
/// Heavy neutral lepton selection
/// \EndBrief                                                                            
/// \Detailed                                                                               
/// Set of cuts for HNL studies, working both on MC samples and data.
/// Several histograms are plotted after each series of cuts.
/// A boolean with the outcome of the selection (whether the HNL passes or not the whole selection)
/// is stored as output and can be retrieved from another analyzer in the following way:
/// \code
/// Bool_t IsGood = *(Bool_t*)GetOutput("HeavyNeutrino".Output);
/// \endcode
/// This analyzer makes use of two ToolsLib, called HNLFunctions and HNLWeight.
/// The value of the squared HNL coupling and the values of the ratios between specific-flavour 
/// couplings can be either set as external parameters from command line or taken as default values.
/// For example, if the user sets USquared = 1.E-10, and UeSquaredRatio = 5., UmuSquaredRatio = 1., 
/// UtauSquaredRatio = 3.5, the specific-flavour coupling values will be: 
/// UeSquared = 5.25E-11, UmuSquared = 1.05E-11, UtauSquared = 3.68E-11.
/// Default values are: USquared = 1.E-6, UeSquaredRatio = 1., UmuSquaredRatio = 16., 
/// UtauSquaredRatio = 3.8.
/// The values of the beginning of the fiducial volume and its length (must be identical to the ones set
/// in the MC macro for production) can be either set as external parameters from command line or 
/// taken as default values.
/// For example, if the user assigns 100000. to the beginning of the FV and 80000. to its length, the FV
/// will begin at 100 m from the target centre and end at 180 m.
/// Default values are: InitialFV = 102500., LFV = 77500.
///
/// \author Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch)
/// \EndDetailed               

#include <stdlib.h>
#include <iostream>
#include <string>
#include <cmath>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "MCInfo.hh"
#include "Persistency.hh"
#include "BeamParameters.hh"
#include "DetectorAcceptance.hh"
#include "SpectrometerTrackVertex.hh"
#include "DownstreamTrack.hh"
#include "GeometricAcceptance.hh"
#include "TriggerConditions.hh"
#include "HNLFunctions.hh"
#include "HNLWeight.hh"
#include "HeavyNeutrino.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

#define TRIGGER_L0_PHYSICS_TYPE 1
#define MCTriggerMask 0xFF
#define LabelSize 0.05

/// \class HeavyNeutrino

HeavyNeutrino::HeavyNeutrino(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrino") {

  if (!GetIsTree()) return;

  RequestAllMCTrees();
  RequestAllRecoTrees();
  RequestL0Data();
  RequestL1Data();

  AddParam("USquared", &fUSquared, 1.E-6);
  AddParam("UeSquaredRatio", &fUeSquaredRatio, 1.);
  AddParam("UmuSquaredRatio", &fUmuSquaredRatio, 16.);
  AddParam("UtauSquaredRatio", &fUtauSquaredRatio, 3.8);
  AddParam("InitialFV", &fInitialFV, 102500.);
  AddParam("LFV", &fLFV, 77500.);

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

  fhNk3pi            = nullptr;
  fhNbursts          = nullptr;
  fhNEvents          = nullptr;
  fhN2tracks         = nullptr;
  fhNtracks          = nullptr;
  fhMomPi            = nullptr;
  fhMomMu            = nullptr;
  fhXYSpec0Reco      = nullptr;
  fhXYSpec1Reco      = nullptr;
  fhXYSpec2Reco      = nullptr;
  fhXYSpec3Reco      = nullptr;
  fhXYCHODReco       = nullptr;
  fhXYCHODTrue       = nullptr;
  fhXYMUV3True       = nullptr;
  fhCuts             = nullptr;
  fhCDAvsZ_In        = nullptr;
  fhCDAvsZ_Track     = nullptr;
  fhCDAvsZ_Energy    = nullptr;
  fhCDAvsZ_Geom      = nullptr;
  fhCDAvsZ_Vetoes    = nullptr;
  fhCDAvsZ_Fin       = nullptr;
  fhZvsBeam_In       = nullptr;
  fhZvsBeam_Track    = nullptr;
  fhZvsBeam_Energy   = nullptr;
  fhZvsBeam_Geom     = nullptr;
  fhZvsBeam_Vetoes   = nullptr;
  fhZvsBeam_Fin      = nullptr;
  fhBeamvsTar_In     = nullptr;
  fhBeamvsTar_Track  = nullptr;
  fhBeamvsTar_Energy = nullptr;
  fhBeamvsTar_Geom   = nullptr;
  fhBeamvsTar_Vetoes = nullptr;
  fhBeamvsTar_Fin    = nullptr;
  fhNMUV3Cand        = nullptr;
  fhEoP              = nullptr;
  fhEoPMuVsPi        = nullptr;
  fhInvMassReco      = nullptr;
}

void HeavyNeutrino::InitOutput() {

  RegisterOutput("Output", &fPassSelection);
}

void HeavyNeutrino::InitHist() {

  BookHisto("hNk3pi",    new TH1D("Nk3pi",    "Total number of K3pi events",       1, 0., 1.));
  BookHisto("hNbursts",  new TH1D("Nbursts",  "Total number of processed bursts",  1, 0., 1.));
  BookHisto("hNEvents",  new TH1D("NEvents",  "Number of total processed events" , 1, 0., 1.));
  BookHisto("hNtracks",  new TH1D("Ntracks",  "Number of tracks",                  4, -0.5, 3.5));
  BookHisto("hN2tracks", new TH1D("N2tracks", "Number of two-tracks events",       1, 0., 1.));
  BookHisto("hMomPi",    new TH1D("MomPi",    "Pion momentum",                     100, -0.5, 200.));
  BookHisto("hMomMu",    new TH1D("MomMu",    "Muon momentum",                     100, -0.5, 200.));
  BookHisto("hCuts",     new TH1D("PhysicsEventsVsCuts", "Physics events passing the selection cuts", 35, 0., 35.));

  BookHisto("hXYSpec0Reco", new TH2D("XYSpec0Reco", "Two-track reconstructed events at CH1",  100, -1.5, 1.5, 100, -1.5, 1.5)); 
  BookHisto("hXYSpec1Reco", new TH2D("XYSpec1Reco", "Two-track reconstructed events at CH2",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYSpec2Reco", new TH2D("XYSpec2Reco", "Two-track reconstructed events at CH3",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYSpec3Reco", new TH2D("XYSpec3Reco", "Two-track reconstructed events at CH4",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYCHODReco",  new TH2D("XYCHODReco",  "Two-track reconstructed events at CHOD", 100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("hXYCHODTrue",  new TH2D("XYCHODTrue",  "X,Y of HNL daughters at CHOD, from MC",  100, -2., 2., 100, -2., 2.));
  BookHisto("hXYMUV3True",  new TH2D("XYMUV3True",  "X,Y of HNL daughters at MUV3, from MC",  100, -2., 2., 100, -2., 2.));

  BookHisto("hCDAvsZ_In",     new TH2D("CDAvsZ_In",     "Two-track total momentum wrt beam axis, before all cuts",          200, 100., 190., 100, 0., 0.5));
  BookHisto("hCDAvsZ_Track",  new TH2D("CDAvsZ_Track",  "Two-track total momentum wrt beam axis, after track-quality cuts", 200, 100., 190., 100, 0., 0.5));
  BookHisto("hCDAvsZ_Energy", new TH2D("CDAvsZ_Energy", "Two-track total momentum wrt beam axis, after energy cuts",        200, 100., 190., 100, 0., 0.5));
  BookHisto("hCDAvsZ_Vetoes", new TH2D("CDAvsZ_Vetoes", "Two-track total momentum wrt beam axis, after veto cuts",          200, 100., 190., 100, 0., 0.5));
  BookHisto("hCDAvsZ_Geom",   new TH2D("CDAvsZ_Geom",   "Two-track total momentum wrt beam axis, after geometrical cuts",   200, 100., 190., 100, 0., 0.5));
  BookHisto("hCDAvsZ_Fin",    new TH2D("CDAvsZ_Fin",    "Two-track total momentum wrt beam axis, after all cuts",           200, 100., 190., 100, 0., 0.5));
  
  BookHisto("hZvsBeam_In",     new TH2D("ZvsBeam_In",     "Two track vertex wrt beam axis, before all cuts",          200, 100., 190., 100, 0., 1.));
  BookHisto("hZvsBeam_Track",  new TH2D("ZvsBeam_Track",  "Two track vertex wrt beam axis, after track-quality cuts", 200, 100., 190., 100, 0., 1.));
  BookHisto("hZvsBeam_Energy", new TH2D("ZvsBeam_Energy", "Two track vertex wrt beam axis, after energy cuts",        200, 100., 190., 100, 0., 1.));
  BookHisto("hZvsBeam_Vetoes", new TH2D("ZvsBeam_Vetoes", "Two track vertex wrt beam axis, after veto cuts",          200, 100., 190., 100, 0., 1.));
  BookHisto("hZvsBeam_Geom",   new TH2D("ZvsBeam_Geom",   "Two track vertex wrt beam axis, after geometrical cuts",   200, 100., 190., 100, 0., 1.));
  BookHisto("hZvsBeam_Fin",    new TH2D("ZvsBeam_Fin",    "Two track vertex wrt beam axis, after all cuts",           200, 100., 190., 100, 0., 1.));
  
  BookHisto("hBeamvsTar_In",     new TH2D("BeamvsTar_In",     "Two-track total momentum wrt beam axis, before all cuts",          100, 0., 1., 200, 0., 27.));
  BookHisto("hBeamvsTar_Track",  new TH2D("BeamvsTar_Track",  "Two-track total momentum wrt beam axis, after track-quality cuts", 100, 0., 1., 200, 0., 27.));
  BookHisto("hBeamvsTar_Energy", new TH2D("BeamvsTar_Energy", "Two-track total momentum wrt beam axis, after energy cuts",        100, 0., 1., 200, 0., 27.));
  BookHisto("hBeamvsTar_Vetoes", new TH2D("BeamvsTar_Vetoes", "Two-track total momentum wrt beam axis, after veto cuts",          100, 0., 1., 200, 0., 27.));
  BookHisto("hBeamvsTar_Geom",   new TH2D("BeamvsTar_Geom",   "Two-track total momentum wrt beam axis, after geometrical cuts",   100, 0., 1., 200, 0., 27.));
  BookHisto("hBeamvsTar_Fin",    new TH2D("BeamvsTar_Fin",    "Two-track total momentum wrt beam axis, after all cuts",           100, 0., 1., 200, 0., 27.));
  
  BookHisto("hNMUV3Cand",   new TH1D("NMUV3Cand", "Number of MUV3 candidates associated to each track", 4, -0.5, 3.5));
  BookHisto("hEoP",         new TH1D("EoP", "E/p in LKr", 100, 0., 1.2));
  BookHisto("hEoPMuVsPi",   new TH2D("EoPMuVsPi", "Muon E/p vs pion E/p in LKr", 100, 0., 1.2, 100, 0., 0.1));  
  BookHisto("hInvMassReco", new TH1D("InvMassReco", "Invariant mass Reco", 50, 960., 1040.));
}

void HeavyNeutrino::Process(Int_t) {

  fPassSelection = false;

  //TRecoLKrEvent*          LKrEvent  = (TRecoLKrEvent*)          GetEvent("LKr");
  TRecoLAVEvent*          LAVEvent  = (TRecoLAVEvent*)          GetEvent("LAV");
  TRecoIRCEvent*          IRCEvent  = (TRecoIRCEvent*)          GetEvent("IRC");
  TRecoSACEvent*          SACEvent  = (TRecoSACEvent*)          GetEvent("SAC");

  // Counter for cuts

  Int_t CutID = 0;
  
  FillHisto("hCuts", CutID);
  CutID++;

  Int_t counter = 0;
  TLorentzVector mom1;
  TLorentzVector mom2;
  Double_t p1,p2;
  Double_t Weight = 0.;

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
	}
      }
    }
  }
  
  // Retrieve weight associated to each HNL
  
  if (GetWithMC()) {
    Event *evt = GetMCEvent();
    
    std::vector<std::map<std::string, Double_t>> Weights = ComputeWeight(evt, fUSquared, fUeSquaredRatio, fUmuSquaredRatio, fUtauSquaredRatio, fLInitialFV, fLFV);

    for (UInt_t i = 0; i < Weights.size(); i++) {
      if ((Bool_t)(Weights[i]["IsGood"]) == true) {
	Weight = Weights[i]["Weight"];    
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
  
  FillHisto("hCuts", CutID);
  CutID++;

  // Select two-track events
  
  std::vector<DownstreamTrack> Tracks = *(std::vector<DownstreamTrack>*) GetOutput("DownstreamTrackBuilder.Output");

  FillHisto("hNtracks", Tracks.size());
  
  if (Tracks.size() != 2)
    return;

  FillHisto("hN2tracks", 0.5);
  FillHisto("hCuts", CutID);
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

  // Compute CDA of track1 wrt track2 and of vertex wrt beam axis
  
  fCDAcomp->SetLine1Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);
  fCDAcomp->SetDir1(Mom1);
  fCDAcomp->SetLine2Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand2->yAtBeforeMagnet(10.0), 10.0);
  fCDAcomp->SetDir2(Mom2);
  fCDAcomp->ComputeVertexCDA();
  
  Double_t CDA     = fCDAcomp->GetCDA();
  TVector3 Vertex  = fCDAcomp->GetVertex();
  Double_t Zvertex = fCDAcomp->GetVertex().z();
  
  fCDAcomp->SetLine1Point1(0.0, 0.0, 102000.0);
  fCDAcomp->SetDir1(0., 0., 1.);
  fCDAcomp->SetLine2Point1(Vertex);
  fCDAcomp->SetDir2(TotMom);
  fCDAcomp->ComputeVertexCDA();
  
  Double_t CDAMom = fCDAcomp->GetCDA();
  
  // Compute distance of two-track momentum wrt target/TAXs
  
  fDistcomp->SetLineDir(TotMom);    
  fDistcomp->SetLinePoint1(Vertex);
  fDistcomp->SetPoint(0., 0., 0.);  
  fDistcomp->ComputeDistance();
  
  Double_t TargetDist = fDistcomp->GetDistance();
  Double_t Extrap     = TargetDist;
  /*
  fDistcomp->SetLineDir(TotMom);    
  fDistcomp->SetLinePoint1(Vertex);
  fDistcomp->SetPoint(0., 0., fTAXDistance + fTAXLength/2.);  
  fDistcomp->ComputeDistance();
  
  Double_t TAXDist = fDistcomp->GetDistance();
  Double_t Extrap = 0.;

  if (TargetDist <= TAXDist)
    Extrap = TargetDist;
  else
    Extrap = TAXDist;
  */

  // Compute distance of two-track vertex wrt beam axis
  
  fDistcomp->SetLinePoint1(0., 0., 102000.);
  fDistcomp->SetLineDir(0., 0., 1.);
  fDistcomp->SetPoint(Vertex);
  
  fDistcomp->ComputeDistance();
  
  Double_t BeamlineDist = fDistcomp->GetDistance();
  
  // Reference plot - 1

  FillHisto("hCDAvsZ_In",   Zvertex/1000.,       CDAMom/1000., Weight);
  FillHisto("hZvsBeam_In",  Zvertex/1000., BeamlineDist/1000., Weight);
  FillHisto("hBeamvsTar_In", Extrap/1000., BeamlineDist/1000., Weight);

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

      FillHisto("hCuts", CutID);
      CutID++;
  }

  // Track selection, CUT 2: Chi2 and momentum cuts

  if (ChiSquare1 >= 20. || ChiSquare2 >= 20.)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  if (SpectrometerCand1->GetNChambers() < 3 || SpectrometerCand2->GetNChambers() < 3)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Track selection, CUT 3: Opposite-charged tracks
  
  if (Charge1 + Charge2 != 0)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Plot the number of candidates associated to each track, for MUV3
  
  FillHisto("hNMUV3Cand", Tracks[0].GetNMUV3AssociationRecords());
  FillHisto("hNMUV3Cand", Tracks[1].GetNMUV3AssociationRecords());
  
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

  FillHisto("hCuts", CutID);
  CutID++;

  if (!CHODAssoc)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Downstream track selection, CUT 5: Extrapolation and association to LKr

  Bool_t LKrAssoc = (Tracks[0].LKrAssociationExists() && Tracks[1].LKrAssociationExists());

  if (!GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kLKr) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kLKr))
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  if (!LKrAssoc)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Downstream track selection, CUT 6: Extrapolation and association to MUV3
    
  if (!GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kMUV3) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kMUV3))
    return;

  FillHisto("hCuts", CutID);
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

  FillHisto("hCuts", CutID);
  CutID++;

  // Reference plot - 2

  FillHisto("hCDAvsZ_Track",   Zvertex/1000.,       CDAMom/1000., Weight);
  FillHisto("hZvsBeam_Track",  Zvertex/1000., BeamlineDist/1000., Weight);
  FillHisto("hBeamvsTar_Track", Extrap/1000., BeamlineDist/1000., Weight);
  
  // Compute time of MUV3 and CHOD candidates for better resolution wrt Spectrometer tracks
  
  Double_t MUV3Time;
  
  if (Assoc == 1)
    MUV3Time = Tracks[0].GetMUV3Time(0);
  else if (Assoc == 2)
    MUV3Time = Tracks[1].GetMUV3Time(0);
  
  Double_t CHODTime1 = Tracks[0].GetCHODTime();
  Double_t CHODTime2 = Tracks[1].GetCHODTime();
  
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

  FillHisto("hCuts", CutID);
  CutID++;

  if (PiEoP <= 0.2 || PiEoP >= 0.8)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Reference plot - 3 

  FillHisto("hCDAvsZ_Energy",   Zvertex/1000.,       CDAMom/1000., Weight);
  FillHisto("hZvsBeam_Energy",  Zvertex/1000., BeamlineDist/1000., Weight);
  FillHisto("hBeamvsTar_Energy", Extrap/1000., BeamlineDist/1000., Weight);

  // Veto cuts, CUT 8: LAV veto
  
  fLAVMatching->SetReferenceTime((CHODTime1 + CHODTime2) / 2);
  
  if (GetWithMC())
    fLAVMatching->SetTimeCuts(99999, 99999);     // LAV is not time-aligned in MC
  else
    fLAVMatching->SetTimeCuts(-10., 10.);
  
  if (fLAVMatching->LAVHasTimeMatching(LAVEvent))
    return;

  FillHisto("hCuts", CutID);
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

  FillHisto("hCuts", CutID);
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

    if (LKrPos && LKrCandE > 40. && LKrTime) {
      LKrEnergyInTime += LKrCandE;
      LKrWeightedTime += LKrCandE * LKrCand->GetTime();
    }
  }
  
  if (LKrEnergyInTime > 1000.) {
    LKrWeightedTime /= LKrEnergyInTime;
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
    
    if (LKrPos && LKrHitE > 40. && LKrTime) {
      LKrEnergyInTime += LKrHitE;
      LKrWeightedTime += LKrHitE * LKrHitTime;
    }
  }
  
  if (LKrEnergyInTime > 1000.) {
    LKrWeightedTime /= LKrEnergyInTime;    
    return;
  }
  */
  FillHisto("hCuts", CutID);
  CutID++;

  // Reference plot - 4 

  FillHisto("hCDAvsZ_Vetoes",   Zvertex/1000.,       CDAMom/1000., Weight);
  FillHisto("hZvsBeam_Vetoes",  Zvertex/1000., BeamlineDist/1000., Weight);
  FillHisto("hBeamvsTar_Vetoes", Extrap/1000., BeamlineDist/1000., Weight);

  // Geometrical cuts, CUT 11: Cut on CDA of two tracks

  if (CDA >= 10.)
    return;

  FillHisto("hCuts", CutID);
  CutID++;
  
  // Geometrical cuts, CUT 12: Cut on two-track vertex wrt beamline

  if (BeamlineDist <= 100.)
    return;

  FillHisto("hCuts", CutID);
  CutID++;
  
  // Geometrical cuts, CUT 13: Cut on Z of two-track vertex

  if (Zvertex <= fInitialFV || Zvertex >= (fInitialFV + fLFV))
    return;

  FillHisto("hCuts", CutID);
  CutID++;
  
  // Reference plot - 5 

  FillHisto("hCDAvsZ_Geom",   Zvertex/1000.,       CDAMom/1000., Weight);
  FillHisto("hZvsBeam_Geom",  Zvertex/1000., BeamlineDist/1000., Weight);
  FillHisto("hBeamvsTar_Geom", Extrap/1000., BeamlineDist/1000., Weight);

  // Reference plot - 6

  FillHisto("hCDAvsZ_Fin",   Zvertex/1000.,       CDAMom/1000., Weight);
  FillHisto("hZvsBeam_Fin",  Zvertex/1000., BeamlineDist/1000., Weight);
  FillHisto("hBeamvsTar_Fin", Extrap/1000., BeamlineDist/1000., Weight);

  fPassSelection = true;

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
}

void HeavyNeutrino::EndOfBurstUser() {
  
  FillHisto("hNbursts", 0.5);
}

void HeavyNeutrino::EndOfJobUser() {
  
  // Retrieve histos

  fhNk3pi            = (TH1D*) fHisto.GetTH1("hNk3pi");
  fhNbursts          = (TH1D*) fHisto.GetTH1("hNbursts");
  fhNEvents          = (TH1D*) fHisto.GetTH1("hNEvents");
  fhN2tracks         = (TH1D*) fHisto.GetTH1("hN2tracks");
  fhNtracks          = (TH1D*) fHisto.GetTH1("hNtracks");
  fhMomPi            = (TH1D*) fHisto.GetTH1("hMomPi");
  fhMomMu            = (TH1D*) fHisto.GetTH1("hMomMu");
  fhXYSpec0Reco      = (TH2D*) fHisto.GetTH2("hXYSpec0Reco");
  fhXYSpec1Reco      = (TH2D*) fHisto.GetTH2("hXYSpec1Reco");
  fhXYSpec2Reco      = (TH2D*) fHisto.GetTH2("hXYSpec2Reco");
  fhXYSpec3Reco      = (TH2D*) fHisto.GetTH2("hXYSpec3Reco");
  fhXYCHODReco       = (TH2D*) fHisto.GetTH2("hXYCHODReco");
  fhXYCHODTrue       = (TH2D*) fHisto.GetTH2("hXYCHODTrue");
  fhXYMUV3True       = (TH2D*) fHisto.GetTH2("hXYMUV3True");
  fhCuts             = (TH1D*) fHisto.GetTH1("hCuts");
  fhCDAvsZ_In        = (TH2D*) fHisto.GetTH2("hCDAvsZ_In");
  fhCDAvsZ_Track     = (TH2D*) fHisto.GetTH2("hCDAvsZ_Track");
  fhCDAvsZ_Energy    = (TH2D*) fHisto.GetTH2("hCDAvsZ_Energy");
  fhCDAvsZ_Vetoes    = (TH2D*) fHisto.GetTH2("hCDAvsZ_Vetoes");
  fhCDAvsZ_Geom      = (TH2D*) fHisto.GetTH2("hCDAvsZ_Geom");
  fhCDAvsZ_Fin       = (TH2D*) fHisto.GetTH2("hCDAvsZ_Fin");
  fhZvsBeam_In       = (TH2D*) fHisto.GetTH2("hZvsBeam_In");
  fhZvsBeam_Track    = (TH2D*) fHisto.GetTH2("hZvsBeam_Track");
  fhZvsBeam_Energy   = (TH2D*) fHisto.GetTH2("hZvsBeam_Energy");
  fhZvsBeam_Vetoes   = (TH2D*) fHisto.GetTH2("hZvsBeam_Vetoes");
  fhZvsBeam_Geom     = (TH2D*) fHisto.GetTH2("hZvsBeam_Geom");
  fhZvsBeam_Fin      = (TH2D*) fHisto.GetTH2("hZvsBeam_Fin");
  fhBeamvsTar_In     = (TH2D*) fHisto.GetTH2("hBeamvsTar_In");
  fhBeamvsTar_Track  = (TH2D*) fHisto.GetTH2("hBeamvsTar_Track");
  fhBeamvsTar_Energy = (TH2D*) fHisto.GetTH2("hBeamvsTar_Energy");
  fhBeamvsTar_Vetoes = (TH2D*) fHisto.GetTH2("hBeamvsTar_Vetoes");
  fhBeamvsTar_Geom   = (TH2D*) fHisto.GetTH2("hBeamvsTar_Geom");
  fhBeamvsTar_Fin    = (TH2D*) fHisto.GetTH2("hBeamvsTar_Fin");
  fhNMUV3Cand        = (TH1D*) fHisto.GetTH1("hNMUV3Cand");
  fhEoP              = (TH1D*) fHisto.GetTH1("hEoP");
  fhEoPMuVsPi        = (TH2D*) fHisto.GetTH1("hEoPMuVsPi");
  fhInvMassReco      = (TH1D*) fHisto.GetTH1("hInvMassReco");

  // X axis title

  fhNk3pi           ->GetXaxis()->SetTitle("Number of k3pi");
  fhNbursts         ->GetXaxis()->SetTitle("Number of bursts");
  fhNEvents         ->GetXaxis()->SetTitle("Number of events");
  fhN2tracks        ->GetXaxis()->SetTitle("Number of two-track events");
  fhNtracks         ->GetXaxis()->SetTitle("Number of tracks in each event");
  fhMomPi           ->GetXaxis()->SetTitle("P [GeV]");
  fhMomMu           ->GetXaxis()->SetTitle("P [GeV]");
  fhXYSpec0Reco     ->GetXaxis()->SetTitle("X [m]");
  fhXYSpec1Reco     ->GetXaxis()->SetTitle("X [m]");
  fhXYSpec2Reco     ->GetXaxis()->SetTitle("X [m]");
  fhXYSpec3Reco     ->GetXaxis()->SetTitle("X [m]");  
  fhXYCHODReco      ->GetXaxis()->SetTitle("X [m]");
  fhXYCHODTrue      ->GetXaxis()->SetTitle("X [m]");
  fhXYMUV3True      ->GetXaxis()->SetTitle("X [m]");
  fhCuts            ->GetXaxis()->SetTitle("Cut ID");    
  fhCDAvsZ_In       ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZ_Track    ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZ_Energy   ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZ_Geom     ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZ_Vetoes   ->GetXaxis()->SetTitle("Z [m]");
  fhCDAvsZ_Fin      ->GetXaxis()->SetTitle("Z [m]");
  fhZvsBeam_In      ->GetXaxis()->SetTitle("Z [m]");
  fhZvsBeam_Track   ->GetXaxis()->SetTitle("Z [m]");
  fhZvsBeam_Energy  ->GetXaxis()->SetTitle("Z [m]");
  fhZvsBeam_Geom    ->GetXaxis()->SetTitle("Z [m]");
  fhZvsBeam_Vetoes  ->GetXaxis()->SetTitle("Z [m]");
  fhZvsBeam_Fin     ->GetXaxis()->SetTitle("Z [m]");
  fhBeamvsTar_In    ->GetXaxis()->SetTitle("Target/TAX distance [m]");
  fhBeamvsTar_Track ->GetXaxis()->SetTitle("Target/TAX distance [m]");
  fhBeamvsTar_Energy->GetXaxis()->SetTitle("Target/TAX distance [m]");
  fhBeamvsTar_Geom  ->GetXaxis()->SetTitle("Target/TAX distance [m]");
  fhBeamvsTar_Vetoes->GetXaxis()->SetTitle("Target/TAX distance [m]");
  fhBeamvsTar_Fin   ->GetXaxis()->SetTitle("Target/TAX distance [m]");
  fhNMUV3Cand       ->GetXaxis()->SetTitle("Number of candidates");
  fhEoP             ->GetXaxis()->SetTitle("E/p");
  fhEoPMuVsPi       ->GetXaxis()->SetTitle("Pion E/p");
  fhInvMassReco     ->GetXaxis()->SetTitle("Invariant mass [MeV]");

  // Y axis title

  fhXYSpec0Reco     ->GetYaxis()->SetTitle("Y [m]");
  fhXYSpec1Reco     ->GetYaxis()->SetTitle("Y [m]");
  fhXYSpec2Reco     ->GetYaxis()->SetTitle("Y [m]");
  fhXYSpec3Reco     ->GetYaxis()->SetTitle("Y [m]");  
  fhXYCHODReco      ->GetYaxis()->SetTitle("Y [m]");
  fhXYCHODTrue      ->GetYaxis()->SetTitle("Y [m]");
  fhXYMUV3True      ->GetYaxis()->SetTitle("Y [m]");
  fhCDAvsZ_In       ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZ_Track    ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZ_Energy   ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZ_Geom     ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZ_Vetoes   ->GetYaxis()->SetTitle("CDA [m]");
  fhCDAvsZ_Fin      ->GetYaxis()->SetTitle("CDA [m]");
  fhZvsBeam_In      ->GetYaxis()->SetTitle("Distance [m]");
  fhZvsBeam_Track   ->GetYaxis()->SetTitle("Distance [m]");
  fhZvsBeam_Energy  ->GetYaxis()->SetTitle("Distance [m]");
  fhZvsBeam_Geom    ->GetYaxis()->SetTitle("Distance [m]");
  fhZvsBeam_Vetoes  ->GetYaxis()->SetTitle("Distance [m]");
  fhZvsBeam_Fin     ->GetYaxis()->SetTitle("Distance [m]");
  fhBeamvsTar_In    ->GetYaxis()->SetTitle("Beam axis distance [m]");
  fhBeamvsTar_Track ->GetYaxis()->SetTitle("Beam axis distance [m]");
  fhBeamvsTar_Energy->GetYaxis()->SetTitle("Beam axis distance [m]");
  fhBeamvsTar_Geom  ->GetYaxis()->SetTitle("Beam axis distance [m]");
  fhBeamvsTar_Vetoes->GetYaxis()->SetTitle("Beam axis distance [m]");
  fhBeamvsTar_Fin   ->GetYaxis()->SetTitle("Beam axis distance [m]");
  fhEoPMuVsPi       ->GetYaxis()->SetTitle("Muon E/p");

  // Plot residual number of events after each cut

  const int NCuts = 25;
  const char *CutNames[NCuts]  = {"Total", "TriggerOK", "2 tracks", "Straw0 acc", "Straw1 acc", "Straw2 acc", "Straw3 acc", "Chi2", "Straw chambers", "Charge", "CHOD acc", "CHOD assoc", "LKr acc", "LKr assoc", "MUV3 acc", "MUV3 assoc", "Mu E/p", "Pi E/p", "LAV veto", "SAV veto", "LKr veto", "CDA tracks", "Beam distance", "Z vertex"};

  for (Int_t i = 1; i <= NCuts; i++)
    fhCuts->GetXaxis()->SetBinLabel(i, CutNames[i-1]);
  
  fhCuts->GetXaxis()->LabelsOption("v");

  SaveAllPlots();

  return;
}

HeavyNeutrino::~HeavyNeutrino() {

  fhNk3pi            = nullptr;
  fhNbursts          = nullptr;
  fhNEvents          = nullptr;
  fhN2tracks         = nullptr;
  fhNtracks          = nullptr;
  fhMomPi            = nullptr;
  fhMomMu            = nullptr;
  fhXYSpec0Reco      = nullptr;
  fhXYSpec1Reco      = nullptr;
  fhXYSpec2Reco      = nullptr;
  fhXYSpec3Reco      = nullptr;
  fhXYCHODReco       = nullptr;
  fhXYCHODTrue       = nullptr;
  fhXYMUV3True       = nullptr;
  fhCuts             = nullptr;
  fhCDAvsZ_In        = nullptr;
  fhCDAvsZ_Track     = nullptr;
  fhCDAvsZ_Energy    = nullptr;
  fhCDAvsZ_Geom      = nullptr;
  fhCDAvsZ_Vetoes    = nullptr;
  fhCDAvsZ_Fin       = nullptr;
  fhZvsBeam_In       = nullptr;
  fhZvsBeam_Track    = nullptr;
  fhZvsBeam_Energy   = nullptr;
  fhZvsBeam_Geom     = nullptr;
  fhZvsBeam_Vetoes   = nullptr;
  fhZvsBeam_Fin      = nullptr;
  fhBeamvsTar_In     = nullptr;
  fhBeamvsTar_Track  = nullptr;
  fhBeamvsTar_Energy = nullptr;
  fhBeamvsTar_Geom   = nullptr;
  fhBeamvsTar_Vetoes = nullptr;
  fhBeamvsTar_Fin    = nullptr;
  fhNMUV3Cand        = nullptr;
  fhEoP              = nullptr;
  fhEoPMuVsPi        = nullptr;
  fhInvMassReco      = nullptr;
}
