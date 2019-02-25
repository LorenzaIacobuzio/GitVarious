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
/// A boolean with the outcome of the selection (whether the HNL passes or not the whole selection) is stored as output and can be retrieved from another analyzer in the following way:
/// \code
/// Bool_t IsGood = *(Bool_t*)GetOutput("HeavyNeutrino".Output);
/// \endcode
/// This analyzer makes use of two ToolsLib, called HNLFunctions and HNLWeight.
/// Several parameters can be set: the value of the squared HNL coupling and the values of the ratios between specific-flavour couplings; the values of the beginning of the fiducial volume and its length; the value of the HNL mass for which the reconstructed invariant mass is computed; the HNL decay mode: 0 for pi-mu final states, 1 for pi-e, 2 for rho-mu and 3 for rho-e; a boolean is set to true if the signal regions are to be kept blinded; a boolean is set to true if MC samples are analysed for studying sidebands (it just changes the weight of the distributions to 1); a boolean that ensures plots are filled only for a certain mass and coupling (set by the user).

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
#include "CalorimeterCluster.hh"
#include "SpectrometerCHODAssociationOutput.hh"
#include "SpectrometerNewCHODAssociationOutput.hh"
#include "SpectrometerLKrAssociationOutput.hh"
#include "SpectrometerMUV3AssociationOutput.hh"
#include "HNLFunctions.hh"
#include "HNLWeight.hh"
#include "HeavyNeutrino.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

#define TRIGGER_L0_PHYSICS_TYPE 1
#define MCTriggerMask 0xFF

/// \class HeavyNeutrino

HeavyNeutrino::HeavyNeutrino(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrino") {

  fReadingData = GetIsTree();
  if (!fReadingData) return;

  RequestAllMCTrees();
  RequestAllRecoTrees();
  RequestL0Data();
  RequestL1Data();
  RequestL0SpecialTrigger();
  RequestBeamSpecialTrigger();

  AddParam("USquared", &fUSquared, 1.E-6); // change accordingly
  AddParam("UInitialeSquaredRatio", &fInitialUeSquaredRatio, 1.); // change accordingly
  AddParam("UInitialmuSquaredRatio", &fInitialUmuSquaredRatio, 16.); // change accordingly
  AddParam("UInitialtauSquaredRatio", &fInitialUtauSquaredRatio, 3.8); // change accordingly
  AddParam("InitialFV", &fInitialFV, 102425.); // keep
  AddParam("LFV", &fLFV, 77575.); // keep
  AddParam("Mode", &fMode, 0);
  AddParam("MassForReco", &fMassForReco, 1.);
  AddParam("BlindRegion", &fBlindRegion, false);
  AddParam("MCsample", &fMCsample, false);
  AddParam("EnableChecks", &fEnableChecks, false);
  AddParam("MassForChecks", &fMassForChecks, 1.6);
  AddParam("CouplingForChecks", &fCouplingForChecks, 1.E-6);

  fUeSquared = fUSquared/(fInitialUeSquaredRatio + fInitialUmuSquaredRatio + fInitialUtauSquaredRatio)*fInitialUeSquaredRatio;
  fUmuSquared = fUSquared/(fInitialUeSquaredRatio + fInitialUmuSquaredRatio + fInitialUtauSquaredRatio)*fInitialUmuSquaredRatio;
  fUtauSquared = fUSquared/(fInitialUeSquaredRatio + fInitialUmuSquaredRatio + fInitialUtauSquaredRatio)*fInitialUtauSquaredRatio;

  fCDAcomp = new TwoLinesCDA();
  fDistcomp = new PointLineDistance();

  fzCHOD = GeometricAcceptance::GetInstance()->GetZCHODVPlane();
  fzMUV3 = GeometricAcceptance::GetInstance()->GetZMUV3();
  fzStraw[0] = GeometricAcceptance::GetInstance()->GetZStraw(0);
  fzStraw[1] = GeometricAcceptance::GetInstance()->GetZStraw(1);
  fzStraw[2] = GeometricAcceptance::GetInstance()->GetZStraw(2);
  fzStraw[3] = GeometricAcceptance::GetInstance()->GetZStraw(3);
  fxStrawChamberCentre[0] = GeometricAcceptance::GetInstance()->GetXStrawChamberCentre(0);
  fxStrawChamberCentre[1] = GeometricAcceptance::GetInstance()->GetXStrawChamberCentre(1);
  fxStrawChamberCentre[2] = GeometricAcceptance::GetInstance()->GetXStrawChamberCentre(2);
  fxStrawChamberCentre[3] = GeometricAcceptance::GetInstance()->GetXStrawChamberCentre(3);
  frMinStraw = GeometricAcceptance::GetInstance()->GetStrawRmin();
  frMaxStraw = GeometricAcceptance::GetInstance()->GetStrawRmax();
  fzCHODPlane = GeometricAcceptance::GetInstance()->GetZCHODVPlane();
  frMinCHOD = GeometricAcceptance::GetInstance()->GetCHODRmin();
  frMaxCHOD = GeometricAcceptance::GetInstance()->GetCHODRmax();

  // Parameters for L0 trigger conditions

  fStream = {"RICH-Q2-MO1", "RICH-Q2-M1", "RICH-Q2-MO1-LKr10", "RICH-Q2-M1-LKr20", "RICH-Q2-MO2-nLKr20", "RICH-Q2-MO2", "RICH-Q2-M2", "RICH-QX-LKr20", "RICH-LKr20", "RICH-Q2-nMUV-LKr20", "RICH-Q2-MO1-LKr20",  "RICH-Q2-MO2-nLKr30", "RICH-QX-MO2"};

  for (UInt_t i = 0; i < fStream.size(); i++) {
    fID.push_back(TriggerConditions::GetInstance()->GetL0TriggerID(fStream[i]));
  }

  // Needed for POT computation

  fNPOTT10 = 0.;
  fNPOTFit = 0.;
  fNKTot = 0.;
  fNK3Pi = 0.;
  fNK = 0.;
  fNPOT = 0.;
  fBurstCounter = 0;
}

void HeavyNeutrino::InitOutput() {

  if (fReadingData) {

    RegisterOutput("Output", &fPassSelection);
    
    OpenNewTree("Passed", "Events");
    
    AddBranch("Passed", "Weight", &Weight);
    AddBranch("Passed", "CHODTime1", &CHODTime1);
    AddBranch("Passed", "CHODTime2", &CHODTime2);
    AddBranch("Passed", "CDA", &CDA);
    AddBranch("Passed", "Zvertex", &Zvertex);
    AddBranch("Passed", "CDALine", &CDALine);
    AddBranch("Passed", "ZCDALine", &ZCDALine);
    AddBranch("Passed", "BeamlineDist", &BeamlineDist);
    AddBranch("Passed", "xSR", &xSR);
    AddBranch("Passed", "ySR", &ySR);
    AddBranch("Passed", "MuEoP", &MuEoP);
    AddBranch("Passed", "PiEoP", &PiEoP);
    AddBranch("Passed", "R", &R);
    AddBranch("Passed", "energyPi", &energyPi);
    AddBranch("Passed", "energyMu", &energyMu);
    AddBranch("Passed", "Mom1", &Mom1);
    AddBranch("Passed", "Mom2", &Mom2);
    AddBranch("Passed", "TotMom", &TotMom);
    AddBranch("Passed", "Vertex", &Vertex);
    AddBranch("Passed", "threeMomPi", &threeMomPi);
    AddBranch("Passed", "threeMomMu", &threeMomMu);
    AddBranch("Passed", "Target", &Target);
    AddBranch("Passed", "Assoc", &Assoc);
  }
  else {
    ImportAllInputHistogram("HeavyNeutrino", false, "HeavyNeutrino");
  }

  return;
}

void HeavyNeutrino::StartOfRunUser() {

  return;
}

void HeavyNeutrino::InitHist() {

  if (fReadingData) {

    BookHisto("hNk3pi",    new TH1D("Nk3pi", "Total number of K3pi events", 1, 0., 1.));
    BookHisto("hNbursts",  new TH1D("Nbursts", "Total number of processed bursts", 1, 0., 1.));
    BookHisto("hNEvents",  new TH1D("NEvents", "Number of total processed events" , 1, 0., 1.));
    BookHisto("hNtracks",  new TH1D("Ntracks", "Number of tracks", 10, -0.5, 9.5));
    BookHisto("hMomPi",    new TH1D("MomPi", "Pion momentum", 100, -0.5, 200.));
    BookHisto("hMomMu",    new TH1D("MomMu", "Muon momentum", 100, -0.5, 200.));
    BookHisto("hCuts",     new TH1D("Cuts", "Physics events after cuts", 37, 0., 37.));

    // X,Y distributions

    BookHisto("hXYSpec0Reco", new TH2D("XYSpec0Reco", "Two-track reconstructed events at CH1", 100, -1.5, 1.5, 100, -1.5, 1.5)); 
    BookHisto("hXYSpec1Reco", new TH2D("XYSpec1Reco", "Two-track reconstructed events at CH2", 100, -1.5, 1.5, 100, -1.5, 1.5));
    BookHisto("hXYSpec2Reco", new TH2D("XYSpec2Reco", "Two-track reconstructed events at CH3", 100, -1.5, 1.5, 100, -1.5, 1.5));
    BookHisto("hXYSpec3Reco", new TH2D("XYSpec3Reco", "Two-track reconstructed events at CH4", 100, -1.5, 1.5, 100, -1.5, 1.5));
    BookHisto("hXYCHODReco",  new TH2D("XYCHODReco",  "Two-track reconstructed events at CHOD", 100, -1.5, 1.5, 100, -1.5, 1.5));
    BookHisto("hXYCHODTrue",  new TH2D("XYCHODTrue",  "X,Y of HNL daughters at CHOD, from MC", 100, -2., 2., 100, -2., 2.));
    BookHisto("hXYMUV3True",  new TH2D("XYMUV3True",  "X,Y of HNL daughters at MUV3, from MC", 100, -2., 2., 100, -2., 2.));
    BookHisto("hXYSpec0Mu",   new TH2D("XYSpec0Mu",  "X,Y of muon daughter at CH1, after track-quality cuts", 100, -2., 2., 100, -2., 2.));
    BookHisto("hXYSpec0Pi",   new TH2D("XYSpec0Pi",  "X,Y of pion daughter at CH1, after track-quality cuts", 100, -2., 2., 100, -2., 2.));

    // Beam vs Z
  
    BookHisto("hZvsBeam_In",     new TH2D("ZvsBeam_In", "Two-track vertex, before all cuts",             50, 100., 200., 50, 0., 1.));
    BookHisto("hZvsBeam_Track",  new TH2D("ZvsBeam_Track", "Two-track vertex, after track-quality cuts", 50, 100., 200., 50, 0., 1.));
    BookHisto("hZvsBeam_Energy", new TH2D("ZvsBeam_Energy", "Two-track vertex, after energy cuts",       50, 100., 200., 50, 0., 1.));
    BookHisto("hZvsBeam_Vetoes", new TH2D("ZvsBeam_Vetoes", "Two-track vertex, after veto cuts",         50, 100., 200., 50, 0., 1.));
    BookHisto("hZvsBeam_Geom",   new TH2D("ZvsBeam_Geom", "Two-track vertex, after geometrical cuts",    50, 100., 200., 50, 0., 1.));
    BookHisto("hZvsBeam_Fin",    new TH2D("ZvsBeam_Fin", "Two-track vertex, after all cuts",             50, 100., 200., 50, 0., 1.));

    // Unique signal region (3-segment line)

    BookHisto("hSR_In",     new TH2D("SR_In",     "Signal region, before all cuts",          500, -50., 50., 50, 0., 0.1));
    BookHisto("hSR_Track",  new TH2D("SR_Track",  "Signal region, after track-quality cuts", 500, -50., 50., 50, 0., 0.1));
    BookHisto("hSR_Energy", new TH2D("SR_Energy", "Signal region, after energy cuts",        500, -50., 50., 50, 0., 0.1));
    BookHisto("hSR_Vetoes", new TH2D("SR_Vetoes", "Signal region, after veto cuts",          500, -50., 50., 50, 0., 0.1));
    BookHisto("hSR_Geom",   new TH2D("SR_Geom",   "Signal region, after geometrical cuts",   500, -50., 50., 50, 0., 0.1));
    BookHisto("hSR_Fin",    new TH2D("SR_Fin",    "Signal region, after all cuts",           500, -50., 50., 50, 0., 0.1));

    BookHisto("hSRTar_In",     new TH2D("SRTar_In",     "Signal region (target), before all cuts",          500, -50., 50., 50, 0., 0.1));
    BookHisto("hSRTar_Track",  new TH2D("SRTar_Track",  "Signal region (target), after track-quality cuts", 500, -50., 50., 50, 0., 0.1));
    BookHisto("hSRTar_Energy", new TH2D("SRTar_Energy", "Signal region (target), after energy cuts",        500, -50., 50., 50, 0., 0.1));
    BookHisto("hSRTar_Vetoes", new TH2D("SRTar_Vetoes", "Signal region (target), after veto cuts",          500, -50., 50., 50, 0., 0.1));
    BookHisto("hSRTar_Geom",   new TH2D("SRTar_Geom",   "Signal region (target), after geometrical cuts",   500, -50., 50., 50, 0., 0.1));
    BookHisto("hSRTar_Fin",    new TH2D("SRTar_Fin",    "Signal region (target), after all cuts",           500, -50., 50., 50, 0., 0.1));

    BookHisto("hSRTAX_In",     new TH2D("SRTAX_In",     "Signal region (TAX), before all cuts",          500, -50., 50., 50, 0., 0.1));
    BookHisto("hSRTAX_Track",  new TH2D("SRTAX_Track",  "Signal region (TAX), after track-quality cuts", 500, -50., 50., 50, 0., 0.1));
    BookHisto("hSRTAX_Energy", new TH2D("SRTAX_Energy", "Signal region (TAX), after energy cuts",        500, -50., 50., 50, 0., 0.1));
    BookHisto("hSRTAX_Vetoes", new TH2D("SRTAX_Vetoes", "Signal region (TAX), after veto cuts",          500, -50., 50., 50, 0., 0.1));
    BookHisto("hSRTAX_Geom",   new TH2D("SRTAX_Geom",   "Signal region (TAX), after geometrical cuts",   500, -50., 50., 50, 0., 0.1));
    BookHisto("hSRTAX_Fin",    new TH2D("SRTAX_Fin",    "Signal region (TAX), after all cuts",           500, -50., 50., 50, 0., 0.1));

    // Unique signal region (one line)

    BookHisto("hCDAvsZCDATarget_In",     new TH2D("CDAvsZCDATarget_In", "N trajectory wrt target-TAX line, before all cuts",             500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDATarget_Track",  new TH2D("CDAvsZCDATarget_Track", "N trajectory wrt target-TAX line, after track-quality cuts", 500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDATarget_Energy", new TH2D("CDAvsZCDATarget_Energy", "N trajectory wrt target-TAX line, after energy cuts",       500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDATarget_Vetoes", new TH2D("CDAvsZCDATarget_Vetoes", "N trajectory wrt target-TAX line, after veto cuts",         500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDATarget_Geom",   new TH2D("CDAvsZCDATarget_Geom", "N trajectory wrt target-TAX line, after geometrical cuts",    500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDATarget_Fin",    new TH2D("CDAvsZCDATarget_Fin", "N trajectory wrt target-TAX line, after all cuts",             500, -50., 50., 50, 0., 0.1));

    BookHisto("hCDAvsZCDATAX_In",     new TH2D("CDAvsZCDATAX_In", "N trajectory wrt target-TAX line, before all cuts",             500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDATAX_Track",  new TH2D("CDAvsZCDATAX_Track", "N trajectory wrt target-TAX line, after track-quality cuts", 500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDATAX_Energy", new TH2D("CDAvsZCDATAX_Energy", "N trajectory wrt target-TAX line, after energy cuts",       500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDATAX_Vetoes", new TH2D("CDAvsZCDATAX_Vetoes", "N trajectory wrt target-TAX line, after veto cuts",         500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDATAX_Geom",   new TH2D("CDAvsZCDATAX_Geom", "N trajectory wrt target-TAX line, after geometrical cuts",    500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDATAX_Fin",    new TH2D("CDAvsZCDATAX_Fin", "N trajectory wrt target-TAX line, after all cuts",             500, -50., 50., 50, 0., 0.1));

    BookHisto("hCDAvsZCDAAll_In",     new TH2D("CDAvsZCDAAll_In", "N trajectory wrt target-TAX line, before all cuts",             500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDAAll_Track",  new TH2D("CDAvsZCDAAll_Track", "N trajectory wrt target-TAX line, after track-quality cuts", 500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDAAll_Energy", new TH2D("CDAvsZCDAAll_Energy", "N trajectory wrt target-TAX line, after energy cuts",       500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDAAll_Vetoes", new TH2D("CDAvsZCDAAll_Vetoes", "N trajectory wrt target-TAX line, after veto cuts",         500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDAAll_Geom",   new TH2D("CDAvsZCDAAll_Geom", "N trajectory wrt target-TAX line, after geometrical cuts",    500, -50., 50., 50, 0., 0.1));
    BookHisto("hCDAvsZCDAAll_Fin",    new TH2D("CDAvsZCDAAll_Fin", "N trajectory wrt target-TAX line, after all cuts",             500, -50., 50., 50, 0., 0.1));

    // Theta vs Z CDA (studies on SR resolution)

    BookHisto("hThetavsZCDA_Tar", new TH2D("ThetavsZCDA_Tar", "N angular distribution for target-produced events", 200, -15., 15., 50, 0., 0.005));
    BookHisto("hThetavsZCDA_TAX", new TH2D("ThetavsZCDA_TAX", "N angular distribution for TAX-produced events", 200, 0., 50., 50, 0., 0.005));


    // Sidebands 

    BookHisto("hCDA",  new TH1D("CDA", "Vertex-beamline distance for parasitic background", 200, 0., 1000.));
    BookHisto("hTime", new TH1D("Time", "Track time difference for combinatorial background", 500, -15., 15.));
    BookHisto("hZ", new TH1D("Z", "Z of vertex for prompt background", 300, 100., 190.));
    BookHisto("hBeamdistvsMass", new TH2D("BeamdistvsMass", "Vertex-beamline distance vs reconstructed mass", 200, 0., 2., 200, 0., 1000.));

    // Others

    BookHisto("hEoP",         new TH1D("EoP", "E/p in LKr", 100, 0., 1.2));
    BookHisto("hEoPMuVsPi",   new TH2D("EoPMuVsPi", "Muon E/p vs pion E/p in LKr", 100, 0., 1.4, 100, 0., 0.3));  
    BookHisto("hInvMassReco", new TH1D("InvMassReco", "Invariant mass Reco", 50, 0.96, 1.04));

    BookHisto("hKTAG",    new TH1D("KTAG", "Trigger time - KTAG candidate time", 100, -30., 30.));
    BookHisto("hCHOD",    new TH1D("CHOD", "Trigger time - CHOD candidate time", 100, -30., 30.));
    BookHisto("hNewCHOD", new TH1D("NewCHOD", "Trigger time - NewCHOD candidate time", 100, -30., 30.));
    BookHisto("hStraw",   new TH1D("Straw", "Trigger time - reference time", 100, -30., 30.));
    BookHisto("hLKr",     new TH1D("LKr", "Trigger time - LKr candidate time", 100, -30., 30.));
    BookHisto("hMUV3",    new TH1D("MUV3", "Trigger time - MUV3 candidate time", 100, -30., 30.));
    BookHisto("hLAV",     new TH1D("LAV", "Trigger time - LAV candidate time", 100, -150., 150.));
    BookHisto("hSAV",     new TH1D("SAV", "Trigger time - SAV candidate time", 100, -150., 150.));
    BookHisto("hCHANTI",  new TH1D("CHANTI", "Trigger time - CHANTI candidate time", 100, -150., 150.));

    BookHisto("hCHANTImult",   new TH1D("CHANTImult", "CHANTI multiplicity in time", 10, 0., 10.));
    BookHisto("hExtraLKrmult", new TH1D("ExtraLKrmult", "Residual LKr multiplicity in time", 10, 0., 10.));

    // POT

    BookHisto("hPOTT10", new TH1D("POTT10", "", 1, 0., 1.));
    BookHisto("hPOTFit", new TH1D("POTFit", "", 1, 0., 1.));

    BookHisto("hSpare1", new TH1D("Spare1", "", 50, 0., 30.));
    BookHisto("hSpare2", new TH2D("Spare2", "", 100, 110., 130., 100, 110., 130.));

    BookHisto("T10", new TGraph());
    fHisto.GetTGraph("T10")->SetNameTitle("T10", "N POT vs N K decays");
    BookHisto("POT1", new TGraph());
    fHisto.GetTGraph("POT1")->SetNameTitle("POT1", "N POT vs burst ID - T10 method");
    BookHisto("POT2", new TGraph());
    fHisto.GetTGraph("POT2")->SetNameTitle("POT2", "N POT vs burst ID - K3Pi method");
    BookHisto("NK", new TGraph());
    fHisto.GetTGraph("NK")->SetNameTitle("NK", "N kaon decays vs burst ID");
  }

  return;
}

void HeavyNeutrino::Process(Int_t) {

  if (!fReadingData) return;
    
  fPassSelection = false;

  // Compute number of processed events

  FillHisto("hNEvents", 0.5);

  TRecoCedarEvent* CedarEvent = (TRecoCedarEvent*)GetEvent("Cedar");
  TRecoCHODEvent* CHODEvent = (TRecoCHODEvent*)GetEvent("CHOD");
  TRecoLAVEvent* LAVEvent = (TRecoLAVEvent*)GetEvent("LAV");
  TRecoIRCEvent* IRCEvent = (TRecoIRCEvent*)GetEvent("IRC");
  TRecoSACEvent* SACEvent = (TRecoSACEvent*)GetEvent("SAC");
  TRecoCHANTIEvent* CHANTIEvent = (TRecoCHANTIEvent*)GetEvent("CHANTI");

  // Counter for cuts

  Int_t CutID = 0;
  
  FillHisto("hCuts", CutID);
  CutID++;

  Int_t counter = 0;
  TLorentzVector mom1;
  TLorentzVector mom2;
  Double_t p1,p2;
  Weight = 1.;

  // Some plots of KinePart quantities

  if (GetWithMC() && !fMCsample) {
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
	  FillHisto("hMomPi", p1/1000.);
	  FillHisto("hMomMu", p2/1000.);
	}
      }
    }
  }

  // Plots with only one mass and coupling (for checks)

  if (fEnableChecks)
    fUSquared = fCouplingForChecks;
  
  // Retrieve weight associated to good HNL

  Double_t MN = 0.;
  
  if (GetWithMC()) {
    Event *evt = GetMCEvent();

    std::vector<std::map<std::string, Double_t>> Weights = ComputeWeight(evt, fUSquared, fInitialUeSquaredRatio, fInitialUmuSquaredRatio, fInitialUtauSquaredRatio, fLInitialFV, fLFV, fMode);

    for (UInt_t i = 0; i < Weights.size(); i++) {
      if ((Bool_t)(Weights[i]["IsGood"]) == true) {
	Weight = Weights[i]["Weight"];
	MN = round(Weights[i]["Mass"])/1000.;
      }  
    }
  }

  // Plots with only one mass and coupling (for checks)
  
  if (fEnableChecks && MN != fMassForChecks)
    return;

  Int_t RunNumber = GetWithMC() ? 0 : GetEventHeader()->GetRunID();
  Bool_t ControlTrigger = TriggerConditions::GetInstance()->IsControlTrigger(GetL0Data());
  Double_t L0TPTime = ControlTrigger ? GetL0Data()->GetPrimitive(kL0TriggerSlot, kL0CHOD).GetFineTime() : GetL0Data()->GetPrimitive(kL0TriggerSlot, kL0RICH).GetFineTime();
  L0TPTime *= TdcCalib;
  Double_t L0Window = 1.22*5.;
  Double_t KTAGWindow = 0.96*5.;
  Double_t NewCHODWindow = 1.81*5.;
  Double_t LKrWindow = 1.96*5.;
  Double_t MUV3Window = 1.55*5.;
  Double_t CHANTIWindow = 2.60*5.;
  Bool_t k3pi = *(Bool_t*) GetOutput("K3piSelection.EventSelected");
  
  // K3Pi

  if (k3pi && 0x10) {
    fNK3Pi++;
    FillHisto("hNk3pi", 0.5);
  }

  // L0 data
  
  L0TPData *L0TPData = GetL0Data();
  UChar_t L0DataType = GetWithMC() ? 0x11 : L0TPData->GetDataType();
  UInt_t L0TriggerFlags = GetWithMC() ? 0xFF : L0TPData->GetTriggerFlags();
  Bool_t PhysicsTriggerOK = (L0DataType & TRIGGER_L0_PHYSICS_TYPE);
  Bool_t TriggerFlagsOK = L0TriggerFlags & MCTriggerMask;
  Bool_t TriggerOK = PhysicsTriggerOK && TriggerFlagsOK;

  // CUT: L0 + L1
  
  if (!TriggerOK) return;
  
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
      std::string L1Algo = (std::string)TriggerConditions::GetInstance()->GetL1TriggerConditionName(RunNumber, fID[i]);
      std::size_t foundKTAG = L1Algo.find("KTAG");
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

  // CUT: Select two-track events
  
  std::vector<DownstreamTrack> Tracks = *(std::vector<DownstreamTrack>*) GetOutput("DownstreamTrackBuilder.Output");

  FillHisto("hNtracks", Tracks.size());
  
  if (Tracks.size() != 2)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Track features
    
  Int_t Charge1 = Tracks[0].GetCharge();
  Int_t Charge2 = Tracks[1].GetCharge();
  Double_t ChiSquare1 = Tracks[0].GetChi2();
  Double_t ChiSquare2 = Tracks[1].GetChi2();
  TRecoSpectrometerCandidate* SpectrometerCand1 = Tracks[0].GetSpectrometerCandidate();
  TRecoSpectrometerCandidate* SpectrometerCand2 = Tracks[1].GetSpectrometerCandidate();
  Mom1 = SpectrometerCand1->GetThreeMomentumBeforeMagnet();
  Mom2 = SpectrometerCand2->GetThreeMomentumBeforeMagnet();
  TotMom = Mom1 + Mom2;
  CHODTime1 = Tracks[0].GetCHODTime();
  CHODTime2 = Tracks[1].GetCHODTime();

  // TIMING PLOTS

  // Track timing                                                             

  if (!GetWithMC()) {
    FillHisto("hStraw", L0TPTime - CHODTime1);
    FillHisto("hStraw", L0TPTime - CHODTime2);
  }

  // KTAG timing

  if (!GetWithMC()) {
    for(int i = 0; i < CedarEvent->GetNCandidates(); i++) {
      TRecoCedarCandidate* cand = (TRecoCedarCandidate*)CedarEvent->GetCandidate(i);
      if (cand->GetNSectors() >= 5)
	FillHisto("hKTAG", cand->GetTime() - L0TPTime);
    }
  }

  // CHOD timing
  
  if (!GetWithMC()) {
    for(int i = 0; i < CHODEvent->GetNCandidates(); i++) {
      TRecoVCandidate* cand = (TRecoVCandidate*)CHODEvent->GetCandidate(i);
      FillHisto("hCHOD", cand->GetTime() - L0TPTime);
    }
  }
  
  // NewCHOD timing
  
  if (!GetWithMC()) {
    std::vector<SpectrometerNewCHODAssociationOutput> SpecNewCHOD = *(std::vector<SpectrometerNewCHODAssociationOutput>*)GetOutput("SpectrometerNewCHODAssociation.Output");
    for (Int_t i = 0; i < SpecNewCHOD[0].GetNAssociationRecords(); i++) {
      Double_t dT = SpecNewCHOD[0].GetAssociationRecord(i)->GetRecoHitTime() - L0TPTime;
      FillHisto("hNewCHOD", dT);
    }
    
    for (Int_t i = 0; i < SpecNewCHOD[1].GetNAssociationRecords(); i++) {
      Double_t dT = SpecNewCHOD[1].GetAssociationRecord(i)->GetRecoHitTime() - L0TPTime;
      FillHisto("hNewCHOD", dT);
    }
  }

  // LKr timing

  if (!GetWithMC()) {
    std::vector<SpectrometerLKrAssociationOutput> SpecLKr = *(std::vector<SpectrometerLKrAssociationOutput>*)GetOutput("SpectrometerLKrAssociation.Output");
    for (UInt_t i = 0; i < SpecLKr[0].GetNAssociationRecords(); i++) {
      Double_t dT = SpecLKr[0].GetAssociationRecord(i)->GetLKrCandidate()->GetClusterTime() - L0TPTime;
      FillHisto("hLKr", dT);
    }

    for (UInt_t i = 0; i < SpecLKr[1].GetNAssociationRecords(); i++) {
      Double_t dT = SpecLKr[1].GetAssociationRecord(i)->GetLKrCandidate()->GetClusterTime() - L0TPTime;
      FillHisto("hLKr", dT);
    }

    Int_t inTime = 0;

    if (SpecLKr[0].GetNAssociationRecords() > 1) {
      for (UInt_t i = 0; i < SpecLKr[0].GetNAssociationRecords(); i++) {
	Double_t dT = SpecLKr[0].GetAssociationRecord(i)->GetLKrCandidate()->GetClusterTime() - L0TPTime;
	if (TMath::Abs(dT) <= LKrWindow)
	  inTime++;
      }
      FillHisto("hExtraLKrmult", inTime-1);
    }

    inTime = 0;

    if (SpecLKr[1].GetNAssociationRecords() > 1) {
      for (UInt_t i = 0; i < SpecLKr[1].GetNAssociationRecords(); i++) {
        Double_t dT = SpecLKr[1].GetAssociationRecord(i)->GetLKrCandidate()->GetClusterTime() - L0TPTime;
	if (TMath::Abs(dT) <= LKrWindow)
          inTime++;
      }
      FillHisto("hExtraLKrmult", inTime-1);
    }
  }

  // MUV3 timing

  if (!GetWithMC()) {
    std::vector<SpectrometerMUV3AssociationOutput> SpecMUV3 = *(std::vector<SpectrometerMUV3AssociationOutput>*)GetOutput("SpectrometerMUV3Association.Output");
    for (Int_t i = 0; i < SpecMUV3[0].GetNAssociationRecords(); i++) {
      Double_t dT = SpecMUV3[0].GetAssociationRecord(i)->GetMuonTime() - L0TPTime;
      FillHisto("hMUV3", dT);
    }
    
    for (Int_t i = 0; i < SpecMUV3[1].GetNAssociationRecords(); i++) {
      Double_t dT = SpecMUV3[1].GetAssociationRecord(i)->GetMuonTime() - L0TPTime;
      FillHisto("hMUV3", dT);
    }
  }

  // LAV timing

  if (!GetWithMC()) {
    for(int i = 0; i < LAVEvent->GetNCandidates(); i++) {
      TRecoVCandidate* cand = (TRecoVCandidate*)LAVEvent->GetCandidate(i);
      FillHisto("hLAV", cand->GetTime() - L0TPTime);
    }
  }

  // SAV timing

  if (!GetWithMC()) {
    for(int i = 0; i < IRCEvent->GetNCandidates(); i++) {
      TRecoVCandidate* cand = (TRecoVCandidate*)IRCEvent->GetCandidate(i);
      FillHisto("hSAV", cand->GetTime() - L0TPTime);
    }
    for(int i = 0; i < SACEvent->GetNCandidates(); i++) {
      TRecoVCandidate* cand = (TRecoVCandidate*)SACEvent->GetCandidate(i);
      FillHisto("hSAV", cand->GetTime() - L0TPTime);
    }
  }

  // CHANTI timing

  if (!GetWithMC()) {
    for(int i = 0; i < CHANTIEvent->GetNCandidates(); i++) {
      TRecoVCandidate* cand = (TRecoVCandidate*)CHANTIEvent->GetCandidate(i);
      FillHisto("hCHANTI", cand->GetTime() - L0TPTime);
    }
  }

  // CHANTI multiplicity

  if (!GetWithMC()) {
    Int_t inTime = 0;
    for(int i = 0; i < CHANTIEvent->GetNCandidates(); i++) {
      TRecoVCandidate* cand = (TRecoVCandidate*)CHANTIEvent->GetCandidate(i);
      if (TMath::Abs(cand->GetTime() - L0TPTime) <= CHANTIWindow)
	inTime++;
    }
    FillHisto("hCHANTImult", inTime);
  }

  // END OF TIMING PLOTS

  // Track timing                                                             

  if (!GetWithMC()) {
    if (TMath::Abs(L0TPTime - CHODTime1) >= L0Window || TMath::Abs(L0TPTime - CHODTime2) >= L0Window)
      return;
  }

  FillHisto("hCuts", CutID);
  CutID++;

  // KTAG timing

  if (!GetWithMC()) {
    for(int i = 0; i < CedarEvent->GetNCandidates(); i++) {
      TRecoCedarCandidate* cand = (TRecoCedarCandidate*)CedarEvent->GetCandidate(i);
      Double_t dT = cand->GetTime() - L0TPTime;
      if (cand->GetNSectors() >= 5 && TMath::Abs(dT) <= KTAGWindow)
	return;
    }
  }

  FillHisto("hCuts", CutID);
  CutID++;

  // (X,Y) of reconstructed tracks for all Spectrometer chambers
    
  for (UInt_t i = 0; i < Tracks.size(); i++) {
    TRecoSpectrometerCandidate* Cand = Tracks[i].GetSpectrometerCandidate();
    for (Int_t j = 0; j < 4; j++) {
      Double_t x = Cand->xAt(fzStraw[j]);
      Double_t y = Cand->yAt(fzStraw[j]);
      TString name = Form("hXYSpec%dReco",j);
      FillHisto(name, x/1000., y/1000.);
    }
    Double_t x = Cand->xAtAfterMagnet(fzCHODPlane);
    Double_t y = Cand->yAtAfterMagnet(fzCHODPlane);
    FillHisto("hXYCHODReco", x/1000., y/1000.);
  }

  // CDA of track1 wrt track2                                                                            

  fCDAcomp = new TwoLinesCDA();
  fCDAcomp->SetLine1Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);
  fCDAcomp->SetDir1(Mom1);
  fCDAcomp->SetLine2Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand2->yAtBeforeMagnet(10.0), 10.0);
  fCDAcomp->SetDir2(Mom2);
  fCDAcomp->ComputeVertexCDA();
  
  CDA = fCDAcomp->GetCDA();
  Vertex = fCDAcomp->GetVertex();
  Zvertex = fCDAcomp->GetVertex().z();

  // Distance of HNL wrt target-TAX line

  fCDAcomp = new TwoLinesCDA();

  fCDAcomp->SetLine1Point1(Vertex);
  fCDAcomp->SetDir1(TotMom);
  fCDAcomp->SetLine2Point1(0., 0., -26.5);
  fCDAcomp->SetLine2Point2(0., -22., 23230.);
  fCDAcomp->ComputeVertexCDA();

  CDALine = fCDAcomp->GetCDA();
  ZCDALine = fCDAcomp->GetVertex().z();

  // Distance of HNL vertex wrt beamline

  fDistcomp = new PointLineDistance();
  
  fDistcomp->SetLinePoint1(0., 0., 102000.);
  fDistcomp->SetLineDir(0., 0., 1.2E-3);
  fDistcomp->SetPoint(Vertex);
  
  fDistcomp->ComputeDistance();
  
  BeamlineDist = fDistcomp->GetDistance();

  // SR studies 1: target line + TAX line -> 2 SR
  // CDA of HNL wrt target line

  fCDAcomp = new TwoLinesCDA();

  fCDAcomp->SetLine1Point1(Vertex);
  fCDAcomp->SetDir1(TotMom);
  fCDAcomp->SetLine2Point1(0., 0., -26.5);
  fCDAcomp->SetDir2(0., 0., 1.);
  fCDAcomp->ComputeVertexCDA();

  Double_t ZVertexMomTarget = fCDAcomp->GetVertex().z();

  // CDA of HNL wrt TAX line

  fCDAcomp = new TwoLinesCDA();

  fCDAcomp->SetLine1Point1(Vertex);
  fCDAcomp->SetDir1(TotMom);
  fCDAcomp->SetLine2Point1(0., -22., 23230.);
  fCDAcomp->SetDir2(0., 0., 1.);
  fCDAcomp->ComputeVertexCDA();

  Double_t ZVertexMomTAX = fCDAcomp->GetVertex().z();

  // SR studies 2: 3-segment line -> 1 SR
  // CDA of hnl wrt target segment

  fCDAcomp = new TwoLinesCDA();

  fCDAcomp->SetLine1Point1(Vertex);
  fCDAcomp->SetDir1(TotMom);
  fCDAcomp->SetLine2Point1(0., 0., 14960.); // Acromat bending
  fCDAcomp->SetDir2(0., 0., 1.);
  fCDAcomp->ComputeVertexCDA();

  Double_t ZBrokenTarget = fCDAcomp->GetVertex().z();
  Double_t CDABrokenTarget = fCDAcomp->GetCDA();

  // CDA of hnl wrt proton line segment

  fCDAcomp = new TwoLinesCDA();

  fCDAcomp->SetLine1Point1(Vertex);
  fCDAcomp->SetDir1(TotMom);
  fCDAcomp->SetLine2Point1(0., 0., 14960.);
  fCDAcomp->SetLine2Point2(0., -22., 22760.); // Two big magnets
  fCDAcomp->ComputeVertexCDA();

  Double_t ZBrokenProtonLine = fCDAcomp->GetVertex().z();
  Double_t CDABrokenProtonLine = fCDAcomp->GetCDA();

  // CDA of hnl wrt TAX segment

  fCDAcomp = new TwoLinesCDA();

  fCDAcomp->SetLine1Point1(Vertex);
  fCDAcomp->SetDir1(TotMom);
  fCDAcomp->SetLine2Point1(0., -22., 22760.);
  fCDAcomp->SetDir2(0., 0., 1.);
  fCDAcomp->ComputeVertexCDA();

  Double_t ZBrokenTAX = fCDAcomp->GetVertex().z();
  Double_t CDABrokenTAX = fCDAcomp->GetCDA();

  xSR = -999.;
  ySR = -999.;
  Double_t minCDA = 999999.;
  Double_t xRes = 3.*12000.;

  if (ZBrokenTarget <= (14960. + xRes) && CDABrokenTarget < minCDA) {
    minCDA = CDABrokenTarget;
    xSR = ZBrokenTarget;
  }
  if (ZBrokenProtonLine > (14960. - xRes) && ZBrokenProtonLine <= (22760. + xRes) && CDABrokenProtonLine < minCDA) {
    minCDA = CDABrokenProtonLine;
    xSR = ZBrokenProtonLine;
  }
  if (ZBrokenTAX > (22760. - xRes) && CDABrokenTAX < minCDA) {
    minCDA = CDABrokenTAX;
    xSR = ZBrokenTAX;
  }

  ySR = minCDA;

  // MC studies on extrapolation to target/TAX

  Target = false;
  Double_t Theta = 0.;

  if (GetWithMC()) {
    Event *evt = GetMCEvent();

    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);
      if (p->GetParentID() == -1 && p->GetPDGcode() == 999 && p->GetEndProcessName() == "good") {  
	Theta = p->GetInitial4Momentum().Theta();
	if (p->GetProdPos().Z() >= -400. && p->GetProdPos().Z() <= 400.)
	  Target = true;
	else
	  Target = false;
      }
    }
  }

  // Reference plot - 1

  Double_t xMin = -15000.;
  Double_t xMax = 40000.;
  Double_t yMin = 0.;
  Double_t yMax = 50.;

  FillHisto("hZvsBeam_In", Zvertex/1000., BeamlineDist/1000., Weight);

  if (GetWithMC()) {
    if (Target == true) {
      FillHisto("hCDAvsZCDATarget_In", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_In", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      FillHisto("hCDAvsZCDATAX_In", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_In", xSR/1000., ySR/1000., Weight);
      }
    }
  }

  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_In", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hCDAvsZCDAAll_In", ZCDALine/1000., CDALine/1000., Weight);

  // Track selection, CUT: Two tracks in Spectrometer acceptance

  if (!GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kSpectrometer, 0) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kSpectrometer, 1) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kSpectrometer, 2) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kSpectrometer, 3) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kSpectrometer, 0) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kSpectrometer, 1) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kSpectrometer, 2) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kSpectrometer, 3))
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Track selection, CUT: Chi2

  if (ChiSquare1 >= 20. || ChiSquare2 >= 20.)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Track selection, CUT: Chambers

  if (SpectrometerCand1->GetNChambers() < 3 || SpectrometerCand2->GetNChambers() < 3)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Track selection, CUT: Opposite-charged tracks
  
  if (Charge1 + Charge2 != 0)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Downstream track selection, CUT: Extrapolation and association to CHOD

  Bool_t CHODAssoc = (Tracks[0].CHODAssociationExists() && Tracks[1].CHODAssociationExists());

  if (!GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kCHOD) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kCHOD))
    return;
  
  FillHisto("hCuts", CutID);
  CutID++;

  if (!CHODAssoc)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // CUT: CHOD timing
  
  if (!GetWithMC()) {
    if (!CHODEvent->GetNCandidates() || CHODTime1 < 0. || CHODTime2 < 0. || Tracks[0].isCHODShowerLikeEvent() || Tracks[1].isCHODShowerLikeEvent())
      return;
  }

  FillHisto("hCuts", CutID);
  CutID++;

  // Downstream track selection, CUT: Extrapolation and association to NewCHOD

  Bool_t NewCHODAssoc = (Tracks[0].NewCHODAssociationExists() && Tracks[1].NewCHODAssociationExists());

  if (!GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kNewCHOD) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kNewCHOD))
    return;
  
  FillHisto("hCuts", CutID);
  CutID++;

  if (!NewCHODAssoc)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // CUT: NewCHOD timing

  if (!GetWithMC()) {
    Int_t inTime = 0;
    std::vector<SpectrometerNewCHODAssociationOutput> SpecNewCHOD = *(std::vector<SpectrometerNewCHODAssociationOutput>*)GetOutput("SpectrometerNewCHODAssociation.Output");
    for (Int_t i = 0; i < SpecNewCHOD[0].GetNAssociationRecords(); i++) {
      Double_t dT = SpecNewCHOD[0].GetAssociationRecord(i)->GetRecoHitTime() - L0TPTime;
      if (TMath::Abs(dT) <= NewCHODWindow)
	inTime++;
    }

    if (!inTime)
      return;

    inTime = 0;

    for (Int_t i = 0; i < SpecNewCHOD[1].GetNAssociationRecords(); i++) {
      Double_t dT = SpecNewCHOD[1].GetAssociationRecord(i)->GetRecoHitTime() - L0TPTime;
      if (TMath::Abs(dT) <= NewCHODWindow)
        inTime++;
    }    

    if (!inTime) 
      return;
  }

  FillHisto("hCuts", CutID);
  CutID++;
  
  // Downstream track selection, CUT: Extrapolation and association to LKr
  
  //Bool_t LKrAssoc = (Tracks[0].LKrAssociationExists() && Tracks[1].LKrAssociationExists());
  
  if (!GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kLKr) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kLKr))
    return;
  
  FillHisto("hCuts", CutID);
  CutID++;
  /*
  if (!LKrAssoc)
    return;

  FillHisto("hCuts", CutID);
  CutID++;
  */

  // Downstream track selection, CUT: Extrapolation and association to MUV3
    
  if (!GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kMUV3) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kMUV3))
    return;

  FillHisto("hCuts", CutID);
  CutID++;
  
  Bool_t Assoc1 = Tracks[0].MUV3AssociationExists();
  Bool_t Assoc2 = Tracks[1].MUV3AssociationExists();
  Int_t NoAssoc = 0;

  Assoc   = 0;
  
  if (Assoc1 && !Assoc2) {
    Assoc = 1;
    NoAssoc = 2;
  }
  else if (!Assoc1 && Assoc2) {
    Assoc = 2;
    NoAssoc = 1;
  }
  else {
    Assoc = 0;
    NoAssoc = 0;
  }

  if (Assoc != 1 && Assoc != 2)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // CUT: MUV3 timing

  if (!GetWithMC()) {
    Int_t inTime = 0;
    std::vector<SpectrometerMUV3AssociationOutput> SpecMUV3 = *(std::vector<SpectrometerMUV3AssociationOutput>*)GetOutput("SpectrometerMUV3Association.Output");
    for (Int_t i = 0; i < SpecMUV3[Assoc-1].GetNAssociationRecords(); i++) {
      Double_t dT = SpecMUV3[Assoc-1].GetAssociationRecord(i)->GetMuonTime() - L0TPTime;
      FillHisto("hMUV3", dT);
      if (TMath::Abs(dT) <= MUV3Window)
	inTime++;
    }
    
    if (!inTime)
      return;
  }

  FillHisto("hCuts", CutID);
  CutID++;
  
  // CUT: LKr timing
  
  if (!GetWithMC()) {
    Int_t inTime = 0;
    std::vector<SpectrometerLKrAssociationOutput> SpecLKr = *(std::vector<SpectrometerLKrAssociationOutput>*)GetOutput("SpectrometerLKrAssociation.Output");
    for (UInt_t i = 0; i < SpecLKr[NoAssoc-1].GetNAssociationRecords(); i++) {
      Double_t dT = SpecLKr[NoAssoc-1].GetAssociationRecord(i)->GetLKrCandidate()->GetClusterTime() - L0TPTime;
      if (TMath::Abs(dT) <= LKrWindow)
        inTime++;
    }
    
    if (!inTime)
      return;
  }
  
  FillHisto("hCuts", CutID);
  CutID++;
  
  // Downstream track selection, CUT: LAV12 acceptance

  if (!GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kLAV12) || !GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kLAV12))
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Reference plot - 2

  FillHisto("hZvsBeam_Track", Zvertex/1000., BeamlineDist/1000., Weight);

  if (GetWithMC()) {
    if (Target == true) {
      FillHisto("hCDAvsZCDATarget_Track", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_Track", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      FillHisto("hCDAvsZCDATAX_Track", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_Track", xSR/1000., ySR/1000., Weight);
      }
    }
  }
  
  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_Track", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hCDAvsZCDAAll_Track", ZCDALine/1000., CDALine/1000., Weight);

  // X,Y of mu/pi at Straw1

  if (Assoc == 1) {  
    Double_t x = SpectrometerCand1->xAt(fzStraw[0]);
    Double_t y = SpectrometerCand1->yAt(fzStraw[0]);
    Double_t r = sqrt(x*x + y*y);
    Double_t rShifted = sqrt(pow(x-fxStrawChamberCentre[0],2) + y*y);
    if (rShifted > frMinStraw && r < frMaxStraw) {
      FillHisto("hXYSpec0Mu", x/1000., y/1000.);
    }
    x = SpectrometerCand2->xAt(fzStraw[0]);
    y = SpectrometerCand2->yAt(fzStraw[0]);
    r = sqrt(x*x + y*y);
    rShifted = sqrt(pow(x-fxStrawChamberCentre[0],2) + y*y);
    if (rShifted > frMinStraw && r < frMaxStraw) {
      FillHisto("hXYSpec0Pi", x/1000., y/1000.);
    }
  }
  else if (Assoc == 2) {
    Double_t x = SpectrometerCand2->xAt(fzStraw[0]);
    Double_t y = SpectrometerCand2->yAt(fzStraw[0]);
    Double_t r = sqrt(x*x + y*y);
    Double_t rShifted = sqrt(pow(x-fxStrawChamberCentre[0],2) + y*y);
    if (rShifted > frMinStraw && r < frMaxStraw) {
      FillHisto("hXYSpec0Mu", x/1000., y/1000.);
    }
    x = SpectrometerCand1->xAt(fzStraw[0]);
    y = SpectrometerCand1->yAt(fzStraw[0]);
    r = sqrt(x*x + y*y);
    rShifted = sqrt(pow(x-fxStrawChamberCentre[0],2) + y*y);
    if (rShifted > frMinStraw && r < frMaxStraw) {
      FillHisto("hXYSpec0Pi", x/1000., y/1000.);
    }
  }
  
  // Energy cuts, CUT: Cut on E/p in LKr
  
  Double_t EoP1 = Tracks[0].GetLKrEoP();
  Double_t EoP2 = Tracks[1].GetLKrEoP();
  MuEoP = 0.;
  PiEoP = 0.;

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

  if (PiEoP >= 0.8)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Reference plot - 3 

  FillHisto("hZvsBeam_Energy", Zvertex/1000., BeamlineDist/1000., Weight);

  if (GetWithMC()) {
    if (Target == true) {
      FillHisto("hCDAvsZCDATarget_Energy", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_Energy", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      FillHisto("hCDAvsZCDATAX_Energy", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_Energy", xSR/1000., ySR/1000., Weight);
      }
    }
  }

  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_Energy", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hCDAvsZCDAAll_Energy", ZCDALine/1000., CDALine/1000., Weight);

  // Veto cuts, CUT: LAV veto
  
  if (!GetWithMC()) {
    fLAVMatching = *(LAVMatching**)GetOutput("PhotonVetoHandler.LAVMatching");
    if (fLAVMatching->LAVHasTimeMatching(LAVEvent))
      return;
  }

  if (GetWithMC()) {
    fLAVMatching = new LAVMatching();
    fLAVMatching->SetReferenceTime((CHODTime1+CHODTime2)/2.);
    fLAVMatching->SetTimeCuts(-99999, 99999);     // LAV is not time-aligned in MC
    if (fLAVMatching->LAVHasTimeMatching(LAVEvent))
      return;
  }

  FillHisto("hCuts", CutID);
  CutID++;
  
  // Veto cuts, CUT: SAV veto
  
  if (!GetWithMC()) {
    fSAVMatching = *(SAVMatching**)GetOutput("PhotonVetoHandler.SAVMatching");
    if (fSAVMatching->SAVHasTimeMatching(IRCEvent, SACEvent))
      return;
  }

  if (GetWithMC()) {
    fSAVMatching = new SAVMatching();
    fSAVMatching->SetReferenceTime((CHODTime1+CHODTime2)/2.);
    fSAVMatching->SetIRCTimeCuts(-99999, 99999);     // SAV is not time-aligned in MC            
    fSAVMatching->SetSACTimeCuts(-99999, 99999);
    if (fSAVMatching->SAVHasTimeMatching(IRCEvent, SACEvent))
      return;
  }
  
  FillHisto("hCuts", CutID);
  CutID++;

  // Veto cuts, CUT: CHANTI veto

  if (!GetWithMC()) {
    for(int i = 0; i < CHANTIEvent->GetNCandidates(); i++) {
      TRecoVCandidate* cand = (TRecoVCandidate*)CHANTIEvent->GetCandidate(i);
      Double_t dT = cand->GetTime() - L0TPTime;
      if (TMath::Abs(dT) <= CHANTIWindow)
	return;
    }
  }

  FillHisto("hCuts", CutID);
  CutID++;

  // Veto cuts, CUT: Residual LKr

  if (!GetWithMC()) {
    std::vector<SpectrometerLKrAssociationOutput> SpecLKr = *(std::vector<SpectrometerLKrAssociationOutput>*)GetOutput("SpectrometerLKrAssociation.Output");
    if (SpecLKr[0].GetNAssociationRecords() > 1) {
      for (UInt_t i = 0; i < SpecLKr[0].GetNAssociationRecords(); i++) {
	Double_t dT = SpecLKr[0].GetAssociationRecord(i)->GetLKrCandidate()->GetClusterTime() - L0TPTime;
	if (TMath::Abs(dT) <= LKrWindow)
	  return;
      }
    }

    if (SpecLKr[1].GetNAssociationRecords() > 1) {
      for (UInt_t i = 0; i < SpecLKr[1].GetNAssociationRecords(); i++) {
        Double_t dT = SpecLKr[1].GetAssociationRecord(i)->GetLKrCandidate()->GetClusterTime() - L0TPTime;
        if (TMath::Abs(dT) <= LKrWindow)
          return;
      }
    }
  }

  FillHisto("hCuts", CutID);
  CutID++;

  // Reference plot - 4 

  FillHisto("hZvsBeam_Vetoes", Zvertex/1000., BeamlineDist/1000., Weight);

  if (GetWithMC()) {
    if (Target == true) {
      FillHisto("hCDAvsZCDATarget_Vetoes", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_Vetoes", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      FillHisto("hCDAvsZCDATAX_Vetoes", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_Vetoes", xSR/1000., ySR/1000., Weight);
      }
    }
  }
  
  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_Vetoes", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hCDAvsZCDAAll_Vetoes", ZCDALine/1000., CDALine/1000., Weight);

  // Geometrical cuts, CUT: Cut on distance between tracks at CH1

  Double_t X = SpectrometerCand1->xAt(fzStraw[0]) - SpectrometerCand2->xAt(fzStraw[0]);
  Double_t Y = SpectrometerCand1->yAt(fzStraw[0]) - SpectrometerCand2->yAt(fzStraw[0]);
  R = TMath::Sqrt(X*X + Y*Y);
  
  if (R < 20.)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Geometrical cuts, CUT: Cut on CDA of two tracks

  if (CDA >= 10.)
    return;

  FillHisto("hCuts", CutID);
  CutID++;
    
  // Geometrical cuts, CUT: Cut on Z of two-track vertex

  if (Zvertex <= fInitialFV || Zvertex >= (fInitialFV + fLFV))
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Studies on SR resolution

  if (Target == true)
    FillHisto("hThetavsZCDA_Tar", ZVertexMomTarget/1000., Theta, Weight);
  else
    FillHisto("hThetavsZCDA_TAX", ZVertexMomTAX/1000., Theta, Weight);

  // Geometrical cuts, CUT: Cut on two-track vertex wrt beamline

  if (BeamlineDist <= 100.)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Reference plot - 5 

  FillHisto("hZvsBeam_Geom", Zvertex/1000., BeamlineDist/1000., Weight);

  if (GetWithMC()) {
    if (Target == true) {
      FillHisto("hCDAvsZCDATarget_Geom", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_Geom", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      FillHisto("hCDAvsZCDATAX_Geom", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_Geom", xSR/1000., ySR/1000., Weight);
      }
    }
  }
  
  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_Geom", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hCDAvsZCDAAll_Geom", ZCDALine/1000., CDALine/1000., Weight);

  // Reference plot - 6

  FillHisto("hZvsBeam_Fin", Zvertex/1000., BeamlineDist/1000., Weight);

  if (GetWithMC()) {
    if (Target == true) {
      FillHisto("hCDAvsZCDATarget_Fin", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_Fin", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      FillHisto("hCDAvsZCDATAX_Fin", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_Fin", xSR/1000., ySR/1000., Weight);
      }
    }
  }
  
  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_Fin", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hCDAvsZCDAAll_Fin", ZCDALine/1000., CDALine/1000., Weight);

  // Sideband studies for parasitic, combinatorial and prompt background

  Double_t x1CDA = 100.;
  Double_t x2CDA = 125.;
  Double_t x1Time = -6.;
  Double_t x2Time = -3.;
  Double_t x3Time = 3.;
  Double_t x4Time = 6.;
  Double_t x1Z = 100000.;
  Double_t x2Z = 120000.;

  if (BeamlineDist > x1CDA && BeamlineDist < x2CDA)
    FillHisto("hCDA", BeamlineDist, Weight);

  if ((CHODTime1 - CHODTime2 > x1Time && CHODTime1 - CHODTime2 < x2Time) || (CHODTime1 - CHODTime2 > x3Time && CHODTime1 - CHODTime2 < x4Time))
    FillHisto("hTime", CHODTime1 - CHODTime2, Weight);
  
  if (Zvertex > x1Z && Zvertex < x2Z)
    FillHisto("hZ", Zvertex/1000., Weight);

  // Computation of invariant mass
  
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

  if (TMath::Abs(invMass - fMassForReco*1000.) <= 10.)
    FillHisto("hInvMassReco", invMass/1000.);

  FillHisto("hBeamdistvsMass", invMass/1000., BeamlineDist, Weight);

  // Output of selection

  fPassSelection = true;
  FillTrees();

  return;
}

void HeavyNeutrino::ProcessEOBEvent() {

  if (!fReadingData) return;

  if (!GetWithMC()) {
    
    // POT computation (works for parasitic and dump mode)    
    
    fNPOT = GetBeamSpecialTrigger()->GetIntensityT10()*1.E11;
    
    if (fNPOT >= 0.)
      fNPOTT10 += fNPOT;
  }

  return;
}

void HeavyNeutrino::EndOfBurstUser() {

  if (!fReadingData) return;
  
  FillHisto("hNbursts", 0.5);

  if (!GetWithMC()) {
    
    // POT computation (works for parasitic mode only)           
    
    Double_t BRK3Pi = 0.05583;
    Double_t AccK3Pi = 0.1536;
    L0TPData *L0TPData = GetL0Data();
    
    if (fNPOT >= 0.) {
      for (UInt_t i = 0; i < fStream.size(); i++) {
	if (TriggerConditions::GetInstance()->L0TriggerOn(GetEventHeader()->GetRunID(), L0TPData, fID[i])) {
	  fNK = TriggerConditions::GetInstance()->GetL0TriggerDownscaling(GetEventHeader()->GetRunID(), TriggerConditions::GetInstance()->GetL0TriggerID("RICH-QX"))*fNK3Pi/(BRK3Pi*AccK3Pi*TriggerConditions::GetInstance()->GetL0TriggerDownscaling(GetEventHeader()->GetRunID(), TriggerConditions::GetInstance()->GetL0TriggerID(fStream[i])));
	}
      }

      fNKaons.push_back(fNK);
      fNKTot += fNK;
      fHisto.GetTGraph("T10")->SetPoint(fBurstCounter, fNK, fNPOT);
      fHisto.GetTGraph("POT1")->SetPoint(fBurstCounter, fBurstCounter, fNPOT);
      fHisto.GetTGraph("NK")->SetPoint(fBurstCounter, fBurstCounter, fNK);
      fBurstCounter++; 
      fNK3Pi = 0.; 
    }
  }
  
  return;
}

void HeavyNeutrino::EndOfJobUser() {

  if (fReadingData) {

    // X axis

    fHisto.GetTH1("hNk3pi")->GetXaxis()->SetTitle("Number of k3pi");
    fHisto.GetTH1("hNbursts")->GetXaxis()->SetTitle("Number of bursts");
    fHisto.GetTH1("hNEvents")->GetXaxis()->SetTitle("Number of events");
    fHisto.GetTH1("hNtracks")->GetXaxis()->SetTitle("Number of tracks");
    fHisto.GetTH1("hMomPi")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("hMomMu")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH2("hXYSpec0Reco")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("hXYSpec1Reco")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("hXYSpec2Reco")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("hXYSpec3Reco")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("hXYCHODReco")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("hXYCHODTrue")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("hXYMUV3True")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("hXYSpec0Mu")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("hXYSpec0Pi")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH1("hCuts")->GetXaxis()->SetTitle("Cut ID");

    fHisto.GetTH2("hZvsBeam_In")->GetXaxis()->SetTitle("Vertex position along Z [m]");
    fHisto.GetTH2("hZvsBeam_Track")->GetXaxis()->SetTitle("Vertex position along Z [m]");
    fHisto.GetTH2("hZvsBeam_Energy")->GetXaxis()->SetTitle("Vertex position along Z [m]");
    fHisto.GetTH2("hZvsBeam_Vetoes")->GetXaxis()->SetTitle("Vertex position along Z [m]");
    fHisto.GetTH2("hZvsBeam_Geom")->GetXaxis()->SetTitle("Vertex position along Z [m]");
    fHisto.GetTH2("hZvsBeam_Fin")->GetXaxis()->SetTitle("Vertex position along Z [m]");

    fHisto.GetTH2("hSR_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSR_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSR_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSR_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSR_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSR_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");

    fHisto.GetTH2("hSRTar_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTar_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTar_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTar_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTar_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTar_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");

    fHisto.GetTH2("hSRTAX_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTAX_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTAX_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTAX_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTAX_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTAX_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt segment line [m]");

    fHisto.GetTH2("hCDAvsZCDATarget_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATarget_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATarget_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATarget_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATarget_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATarget_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");

    fHisto.GetTH2("hCDAvsZCDATAX_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATAX_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATAX_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATAX_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATAX_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATAX_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");

    fHisto.GetTH2("hCDAvsZCDAAll_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDAAll_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDAAll_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDAAll_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDAAll_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDAAll_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt target-TAX line [m]");

    fHisto.GetTH2("hThetavsZCDA_Tar")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
    fHisto.GetTH2("hThetavsZCDA_TAX")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");

    fHisto.GetTH1("hCDA")->GetXaxis()->SetTitle("Vertex-beamline distance [mm]");
    fHisto.GetTH1("hTime")->GetXaxis()->SetTitle("Track time difference [ns]");
    fHisto.GetTH1("hZ")->GetXaxis()->SetTitle("Z of vertex [m]");
    fHisto.GetTH2("hBeamdistvsMass")->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");

    fHisto.GetTH1("hEoP")->GetXaxis()->SetTitle("E/p");
    fHisto.GetTH1("hEoPMuVsPi")->GetXaxis()->SetTitle("Pion E/p");
    fHisto.GetTH1("hInvMassReco")->GetXaxis()->SetTitle("Invariant mass [GeV/c^{2}]");
    fHisto.GetTH1("hKTAG")->GetXaxis()->SetTitle("Time difference [ns]");
    fHisto.GetTH1("hCHOD")->GetXaxis()->SetTitle("Time difference [ns]");
    fHisto.GetTH1("hNewCHOD")->GetXaxis()->SetTitle("Time difference [ns]");
    fHisto.GetTH1("hLKr")->GetXaxis()->SetTitle("Time difference [ns]");
    fHisto.GetTH1("hMUV3")->GetXaxis()->SetTitle("Time difference [ns]");
    fHisto.GetTH1("hStraw")->GetXaxis()->SetTitle("Time difference [ns]");
    fHisto.GetTH1("hLAV")->GetXaxis()->SetTitle("Time difference [ns]");
    fHisto.GetTH1("hSAV")->GetXaxis()->SetTitle("Time difference [ns]");
    fHisto.GetTH1("hCHANTI")->GetXaxis()->SetTitle("Time difference [ns]");
    fHisto.GetTH1("hCHANTImult")->GetXaxis()->SetTitle("CHANTI multiplicity");
    fHisto.GetTH1("hExtraLKrmult")->GetXaxis()->SetTitle("Residual LKr multiplicity");
    fHisto.GetTGraph("T10")->GetXaxis()->SetTitle("N of K decays per burst");
    fHisto.GetTGraph("POT1")->GetXaxis()->SetTitle("Burst ID");
    fHisto.GetTGraph("POT2")->GetXaxis()->SetTitle("Burst ID");

    // Y axis

    fHisto.GetTH2("hXYSpec0Reco")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("hXYSpec1Reco")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("hXYSpec2Reco")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("hXYSpec3Reco")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("hXYCHODReco")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("hXYCHODTrue")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("hXYMUV3True")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("hXYSpec0Mu")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("hXYSpec0Pi")->GetYaxis()->SetTitle("Y [m]");

    fHisto.GetTH2("hZvsBeam_In")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
    fHisto.GetTH2("hZvsBeam_Track")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
    fHisto.GetTH2("hZvsBeam_Energy")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
    fHisto.GetTH2("hZvsBeam_Vetoes")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
    fHisto.GetTH2("hZvsBeam_Geom")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
    fHisto.GetTH2("hZvsBeam_Fin")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");

    fHisto.GetTH2("hSR_In")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSR_Track")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSR_Energy")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSR_Vetoes")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSR_Geom")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSR_Fin")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");

    fHisto.GetTH2("hSRTar_In")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTar_Track")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTar_Energy")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTar_Vetoes")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTar_Geom")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTar_Fin")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");

    fHisto.GetTH2("hSRTAX_In")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTAX_Track")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTAX_Energy")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTAX_Vetoes")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTAX_Geom")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");
    fHisto.GetTH2("hSRTAX_Fin")->GetYaxis()->SetTitle("CDA of mother wrt segment line [m]");

    fHisto.GetTH2("hCDAvsZCDATarget_In")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATarget_Track")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATarget_Energy")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATarget_Vetoes")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATarget_Geom")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATarget_Fin")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");

    fHisto.GetTH2("hCDAvsZCDATAX_In")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATAX_Track")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATAX_Energy")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATAX_Vetoes")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATAX_Geom")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDATAX_Fin")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");

    fHisto.GetTH2("hCDAvsZCDAAll_In")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDAAll_Track")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDAAll_Energy")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDAAll_Vetoes")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDAAll_Geom")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");
    fHisto.GetTH2("hCDAvsZCDAAll_Fin")->GetYaxis()->SetTitle("CDA of mother wrt target-TAX line [m]");

    fHisto.GetTH2("hThetavsZCDA_Tar")->GetYaxis()->SetTitle("N theta [mrad]");
    fHisto.GetTH2("hThetavsZCDA_TAX")->GetYaxis()->SetTitle("N theta [mrad]");

    fHisto.GetTH2("hBeamdistvsMass")->GetYaxis()->SetTitle("Vertex-beamline distance [mm]");

    fHisto.GetTH1("hEoPMuVsPi")->GetYaxis()->SetTitle("Muon E/p");
    fHisto.GetTGraph("T10")->GetYaxis()->SetTitle("N POT");
    fHisto.GetTGraph("POT1")->GetYaxis()->SetTitle("N POT");
    fHisto.GetTGraph("POT2")->GetYaxis()->SetTitle("N POT");

    // Colz axis

    fHisto.GetTH2("hZvsBeam_In")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");
    fHisto.GetTH2("hZvsBeam_Track")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");
    fHisto.GetTH2("hZvsBeam_Energy")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");
    fHisto.GetTH2("hZvsBeam_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");
    fHisto.GetTH2("hZvsBeam_Geom")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");
    fHisto.GetTH2("hZvsBeam_Fin")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");

    fHisto.GetTH2("hSR_In")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSR_Track")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSR_Energy")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSR_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSR_Geom")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSR_Fin")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");

    fHisto.GetTH2("hSRTar_In")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSRTar_Track")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSRTar_Energy")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSRTar_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSRTar_Geom")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSRTar_Fin")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");

    fHisto.GetTH2("hSRTAX_In")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSRTAX_Track")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSRTAX_Energy")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSRTAX_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSRTAX_Geom")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("hSRTAX_Fin")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");

    fHisto.GetTH2("hCDAvsZCDATarget_In")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDATarget_Track")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDATarget_Energy")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDATarget_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDATarget_Geom")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDATarget_Fin")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");

    fHisto.GetTH2("hCDAvsZCDATAX_In")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDATAX_Track")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDATAX_Energy")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDATAX_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDATAX_Geom")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDATAX_Fin")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");

    fHisto.GetTH2("hCDAvsZCDAAll_In")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDAAll_Track")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDAAll_Energy")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDAAll_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDAAll_Geom")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");
    fHisto.GetTH2("hCDAvsZCDAAll_Fin")->GetZaxis()->SetTitle("Normalized to POT and 50 cm^{2}");

    fHisto.GetTH2("hThetavsZCDA_Tar")->GetZaxis()->SetTitle("Normalized to POT and 0.2 cm^{2}");
    fHisto.GetTH2("hThetavsZCDA_TAX")->GetZaxis()->SetTitle("Normalized to POT and 0.2 cm^{2}");

    fHisto.GetTH2("hBeamdistvsMass")->GetZaxis()->SetTitle("Normalized to POT, 10 MeV/c^{2} and 0.5 cm");

    if (!GetWithMC()) {
      fHisto.GetTH1("hKTAG")->Fit("gaus", "", "", -2., 2.);
      fHisto.GetTH1("hCHOD")->Fit("gaus", "", "", -4., 4.);
      fHisto.GetTH1("hNewCHOD")->Fit("gaus", "", "", -4., 4.);
      fHisto.GetTH1("hLKr")->Fit("gaus", "", "", -2., 6.);
      fHisto.GetTH1("hMUV3")->Fit("gaus", "", "", -4., 4.);
      fHisto.GetTH1("hStraw")->Fit("gaus", "", "", -2., 2.);
      fHisto.GetTH1("hLAV")->Fit("gaus", "", "", -10., 10.);
      fHisto.GetTH1("hSAV")->Fit("gaus", "", "", -10., 10.);
      fHisto.GetTH1("hCHANTI")->Fit("gaus", "", "", -10., 10.);
    
      // POT computation (works for parasitic mode only)

      TF1 *func = new TF1("func", "[0]*x");
      fHisto.GetTGraph("T10")->Fit(func);
      fHisto.GetTGraph("T10")->Draw("AP*");
    
      for (Int_t i = 0; i < fHisto.GetTGraph("POT1")->GetN(); i++) {
	fHisto.GetTGraph("POT2")->SetPoint(i, i, fHisto.GetTGraph("T10")->GetFunction("func")->GetParameter(0)*fNKaons[i]);
	Double_t x;
	Double_t y;
	fHisto.GetTGraph("POT1")->GetPoint(i,x,y);
      }
    
      fNPOTFit = fHisto.GetTGraph("T10")->GetFunction("func")->GetParameter(0)*fNKTot;
      FillHisto("hPOTT10", 0.5, fNPOTT10);
      FillHisto("hPOTFit", 0.5, fNPOTFit);
    
      cout << endl << "Total number of K decays, POT (T10 method) and POT (fit method) in this job: " << fNKTot << " " << fNPOTT10 << " " << fNPOTFit << endl;
    }

    // Plot residual number of events after each cut

    const int NCuts = 37;
    const char *CutNames[NCuts] = {"Total", "TriggerOK", "2 tracks", "Track time", "KTAG time", "Straw acc", "Chi2", "Straw chambers", "Charge", "CHOD acc", "CHOD assoc", "CHOD time", "NewCHOD acc", "NewCHOD assoc", "NewCHOD time", "LKr acc", /*"LKr assoc",*/ "MUV3 acc", "MUV3 assoc", "MUV3 time", "LKr time", "LAV12 acc", "Mu E/p", "Pi E/p", "LAV veto", "SAV veto", "CHANTI veto", "LKr veto", "Track dist CH1", "CDA tracks", "Z vertex", "Beam dist"};
  
    for (Int_t i = 1; i <= NCuts; i++)
      fHisto.GetTH1("hCuts")->GetXaxis()->SetBinLabel(i, CutNames[i-1]);
  
    fHisto.GetTH1("hCuts")->GetXaxis()->LabelsOption("v");
  }
  else {
    TTree* tree = static_cast<TTree*>(GetCurrentFile()->Get("Passed"))->CloneTree();
    tree->Write();
  }
  
  SaveAllPlots();

  return;
}