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

  if (!fReadingData) return;

  RegisterOutput("Output", &fPassSelection);

  return;
}

void HeavyNeutrino::InitHist() {

  if (fReadingData) {

    BookHisto("hNk3pi",    new TH1D("Nk3pi", "Total number of K3pi events", 1, 0., 1.));
    BookHisto("hNbursts",  new TH1D("Nbursts", "Total number of processed bursts", 1, 0., 1.));
    BookHisto("hNEvents",  new TH1D("NEvents", "Number of total processed events" , 1, 0., 1.));
    BookHisto("hNtracks",  new TH1D("Ntracks", "Number of tracks", 10, -0.5, 9.5));
    BookHisto("hN2tracks", new TH1D("N2tracks", "Number of two-tracks events", 1, 0., 1.));
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

    // CDA vs Z

    BookHisto("hCDAvsZ_In",     new TH2D("CDAvsZ_In", "N trajectory wrt beam axis, before all cuts",             50, 100., 200., 30, 0., 0.15));
    BookHisto("hCDAvsZ_Track",  new TH2D("CDAvsZ_Track", "N trajectory wrt beam axis, after track-quality cuts", 50, 100., 200., 30, 0., 0.15));
    BookHisto("hCDAvsZ_Energy", new TH2D("CDAvsZ_Energy", "N trajectory wrt beam axis, after energy cuts",       50, 100., 200., 30, 0., 0.15));
    BookHisto("hCDAvsZ_Vetoes", new TH2D("CDAvsZ_Vetoes", "N trajectory wrt beam axis, after veto cuts",         50, 100., 200., 30, 0., 0.15));
    BookHisto("hCDAvsZ_Geom",   new TH2D("CDAvsZ_Geom", "N trajectory wrt beam axis, after geometrical cuts",    50, 100., 200., 30, 0., 0.15));
    BookHisto("hCDAvsZ_Fin",    new TH2D("CDAvsZ_Fin", "N trajectory wrt beam axis, after all cuts",             50, 100., 200., 30, 0., 0.15));

    // CDA vs CDA

    BookHisto("hCDAvsCDA_Track",  new TH2D("CDAvsCDA_Track", "Track wrt beam axis, after track-quality cuts", 50, 0., 0.5, 50, 0., 0.5));
    BookHisto("hCDAvsCDA_Energy", new TH2D("CDAvsCDA_Energy", "Track wrt beam axis, after energy cuts",       50, 0., 0.5, 50, 0., 0.5));
    BookHisto("hCDAvsCDA_Vetoes", new TH2D("CDAvsCDA_Vetoes", "Track wrt beam axis, after veto cuts",         50, 0., 0.5, 50, 0., 0.5));
    BookHisto("hCDAvsCDA_Geom",   new TH2D("CDAvsCDA_Geom", "Track wrt beam axis, after geometrical cuts",    50, 0., 0.5, 50, 0., 0.5));
    BookHisto("hCDAvsCDA_Fin",    new TH2D("CDAvsCDA_Fin", "Track wrt beam axis, after all cuts",             50, 0., 0.5, 50, 0., 0.5));

    // Beam vs Z
  
    BookHisto("hZvsBeam_In",     new TH2D("ZvsBeam_In", "Two-track vertex, before all cuts",             50, 100., 200., 50, 0., 1.));
    BookHisto("hZvsBeam_Track",  new TH2D("ZvsBeam_Track", "Two-track vertex, after track-quality cuts", 50, 100., 200., 50, 0., 1.));
    BookHisto("hZvsBeam_Energy", new TH2D("ZvsBeam_Energy", "Two-track vertex, after energy cuts",       50, 100., 200., 50, 0., 1.));
    BookHisto("hZvsBeam_Vetoes", new TH2D("ZvsBeam_Vetoes", "Two-track vertex, after veto cuts",         50, 100., 200., 50, 0., 1.));
    BookHisto("hZvsBeam_Geom",   new TH2D("ZvsBeam_Geom", "Two-track vertex, after geometrical cuts",    50, 100., 200., 50, 0., 1.));
    BookHisto("hZvsBeam_Fin",    new TH2D("ZvsBeam_Fin", "Two-track vertex, after all cuts",             50, 100., 200., 50, 0., 1.));

    // Beam vs target
  
    BookHisto("hBeamvsTar_In",     new TH2D("BeamvsTar_In", "N trajectory (target), before all cuts",             50, 0., 0.5, 50, 0., 1.)); // target extrapolation for target-produced events
    BookHisto("hBeamvsTar_Track",  new TH2D("BeamvsTar_Track", "N trajectory (target), after track-quality cuts", 50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTar_Energy", new TH2D("BeamvsTar_Energy", "N trajectory (target), after energy cuts",       50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTar_Vetoes", new TH2D("BeamvsTar_Vetoes", "N trajectory (target), after veto cuts",         50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTar_Geom",   new TH2D("BeamvsTar_Geom", "N trajectory (target), after geometrical cuts",    50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTar_Fin",    new TH2D("BeamvsTar_Fin", "N trajectory (target), after all cuts",             50, 0., 0.5, 50, 0., 1.));

    BookHisto("hBeamvsTarMismatched_In",     new TH2D("BeamvsTarMismatched_In", "N trajectory (target), before all cuts",             50, 0., 0.5, 50, 0., 1.)); // target extrapolation for TAX-produced events
    BookHisto("hBeamvsTarMismatched_Track",  new TH2D("BeamvsTarMismatched_Track", "N trajectory (target), after track-quality cuts", 50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTarMismatched_Energy", new TH2D("BeamvsTarMismatched_Energy", "N trajectory (target), after energy cuts",       50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTarMismatched_Vetoes", new TH2D("BeamvsTarMismatched_Vetoes", "N trajectory (target), after veto cuts",         50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTarMismatched_Geom",   new TH2D("BeamvsTarMismatched_Geom", "N trajectory (target), after geometrical cuts",    50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTarMismatched_Fin",    new TH2D("BeamvsTarMismatched_Fin", "N trajectory (target), after all cuts",             50, 0., 0.5, 50, 0., 1.));

    BookHisto("hBeamvsTarAll_In",     new TH2D("BeamvsTarAll_In", "N trajectory (target), before all cuts",             50, 0., 0.5, 50, 0., 1.)); // target extrapolation for target- and TAX-produced events
    BookHisto("hBeamvsTarAll_Track",  new TH2D("BeamvsTarAll_Track", "N trajectory (target), after track-quality cuts", 50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTarAll_Energy", new TH2D("BeamvsTarAll_Energy", "N trajectory (target), after energy cuts",       50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTarAll_Vetoes", new TH2D("BeamvsTarAll_Vetoes", "N trajectory (target), after veto cuts",         50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTarAll_Geom",   new TH2D("BeamvsTarAll_Geom", "N trajectory (target), after geometrical cuts",    50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTarAll_Fin",    new TH2D("BeamvsTarAll_Fin", "N trajectory (target), after all cuts",             50, 0., 0.5, 50, 0., 1.));

    // Beam vs TAX

    BookHisto("hBeamvsTAX_In",     new TH2D("BeamvsTAX_In", "N trajectory (TAX), before all cuts",             50, 0., 0.5, 50, 0., 1.)); // TAX extrapolation for TAX-produced events
    BookHisto("hBeamvsTAX_Track",  new TH2D("BeamvsTAX_Track", "N trajectory (TAX), after track-quality cuts", 50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAX_Energy", new TH2D("BeamvsTAX_Energy", "N trajectory (TAX), after energy cuts",       50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAX_Vetoes", new TH2D("BeamvsTAX_Vetoes", "N trajectory (TAX), after veto cuts",         50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAX_Geom",   new TH2D("BeamvsTAX_Geom", "N trajectory (TAX), after geometrical cuts",    50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAX_Fin",    new TH2D("BeamvsTAX_Fin", "N trajectory (TAX), after all cuts",             50, 0., 0.5, 50, 0., 1.));

    BookHisto("hBeamvsTAXMismatched_In",     new TH2D("BeamvsTAXMismatched_In", "N trajectory (TAX), before all cuts",             50, 0., 0.5, 50, 0., 1.)); // TAX extrapolation for target-produced events
    BookHisto("hBeamvsTAXMismatched_Track",  new TH2D("BeamvsTAXMismatched_Track", "N trajectory (TAX), after track-quality cuts", 50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAXMismatched_Energy", new TH2D("BeamvsTAXMismatched_Energy", "N trajectory (TAX), after energy cuts",       50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAXMismatched_Vetoes", new TH2D("BeamvsTAXMismatched_Vetoes", "N trajectory (TAX), after veto cuts",         50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAXMismatched_Geom",   new TH2D("BeamvsTAXMismatched_Geom", "N trajectory (TAX), after geometrical cuts",    50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAXMismatched_Fin",    new TH2D("BeamvsTAXMismatched_Fin", "N trajectory (TAX), after all cuts",             50, 0., 0.5, 50, 0., 1.));

    BookHisto("hBeamvsTAXAll_In",     new TH2D("BeamvsTAXAll_In", "N trajectory (TAX), before all cuts",             50, 0., 0.5, 50, 0., 1.)); // TAX extrapolation for target- and TAX-produced events
    BookHisto("hBeamvsTAXAll_Track",  new TH2D("BeamvsTAXAll_Track", "N trajectory (TAX), after track-quality cuts", 50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAXAll_Energy", new TH2D("BeamvsTAXAll_Energy", "N trajectory (TAX), after energy cuts",       50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAXAll_Vetoes", new TH2D("BeamvsTAXAll_Vetoes", "N trajectory (TAX), after veto cuts",         50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAXAll_Geom",   new TH2D("BeamvsTAXAll_Geom", "N trajectory (TAX), after geometrical cuts",    50, 0., 0.5, 50, 0., 1.));
    BookHisto("hBeamvsTAXAll_Fin",    new TH2D("BeamvsTAXAll_Fin", "N trajectory (TAX), after all cuts",             50, 0., 0.5, 50, 0., 1.));

    // Separated signal regions

    BookHisto("hSignalRegionTar_In",     new TH2D("SignalRegionTar_In", "Signal region (target), before all cuts",             500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTar_Track",  new TH2D("SignalRegionTar_Track", "Signal region (target), after track-quality cuts", 500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTar_Energy", new TH2D("SignalRegionTar_Energy", "Signal region (target), after energy cuts",       500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTar_Vetoes", new TH2D("SignalRegionTar_Vetoes", "Signal region (target), after veto cuts",         500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTar_Geom",   new TH2D("SignalRegionTar_Geom", "Signal region (target), after geometrical cuts",    500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTar_Fin",    new TH2D("SignalRegionTar_Fin", "Signal region (target), after all cuts",             500, -100., 100., 80, 0., 0.2));

    BookHisto("hSignalRegionTarMismatched_In",     new TH2D("SignalRegionTarMismatched_In", "Signal region (target), before all cuts",             500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTarMismatched_Track",  new TH2D("SignalRegionTarMismatched_Track", "Signal region (target), after track-quality cuts", 500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTarMismatched_Energy", new TH2D("SignalRegionTarMismatched_Energy", "Signal region (target), after energy cuts",       500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTarMismatched_Vetoes", new TH2D("SignalRegionTarMismatched_Vetoes", "Signal region (target), after veto cuts",         500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTarMismatched_Geom",   new TH2D("SignalRegionTarMismatched_Geom", "Signal region (target), after geometrical cuts",    500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTarMismatched_Fin",    new TH2D("SignalRegionTarMismatched_Fin", "Signal region (target), after all cuts",             500, -100., 100., 80, 0., 0.2));

    BookHisto("hSignalRegionTarAll_In",     new TH2D("SignalRegionTarAll_In", "Signal region (target), before all cuts",             500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTarAll_Track",  new TH2D("SignalRegionTarAll_Track", "Signal region (target), after track-quality cuts", 500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTarAll_Energy", new TH2D("SignalRegionTarAll_Energy", "Signal region (target), after energy cuts",       500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTarAll_Vetoes", new TH2D("SignalRegionTarAll_Vetoes", "Signal region (target), after veto cuts",         500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTarAll_Geom",   new TH2D("SignalRegionTarAll_Geom", "Signal region (target), after geometrical cuts",    500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTarAll_Fin",    new TH2D("SignalRegionTarAll_Fin", "Signal region (target), after all cuts",             500, -100., 100., 80, 0., 0.2));

    BookHisto("hSignalRegionTAX_In",     new TH2D("SignalRegionTAX_In", "Signal region (TAX), before all cuts",             500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAX_Track",  new TH2D("SignalRegionTAX_Track", "Signal region (TAX), after track-quality cuts", 500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAX_Energy", new TH2D("SignalRegionTAX_Energy", "Signal region (TAX), after energy cuts",       500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAX_Vetoes", new TH2D("SignalRegionTAX_Vetoes", "Signal region (TAX), after veto cuts",         500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAX_Geom",   new TH2D("SignalRegionTAX_Geom", "Signal region (TAX), after geometrical cuts",    500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAX_Fin",    new TH2D("SignalRegionTAX_Fin", "Signal region (TAX), after all cuts",             500, -100., 100., 80, 0., 0.2));

    BookHisto("hSignalRegionTAXMismatched_In",     new TH2D("SignalRegionTAXMismatched_In", "Signal region (TAX), before all cuts",             500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAXMismatched_Track",  new TH2D("SignalRegionTAXMismatched_Track", "Signal region (TAX), after track-quality cuts", 500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAXMismatched_Energy", new TH2D("SignalRegionTAXMismatched_Energy", "Signal region (TAX), after energy cuts",       500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAXMismatched_Vetoes", new TH2D("SignalRegionTAXMismatched_Vetoes", "Signal region (TAX), after veto cuts",         500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAXMismatched_Geom",   new TH2D("SignalRegionTAXMismatched_Geom", "Signal region (TAX), after geometrical cuts",    500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAXMismatched_Fin",    new TH2D("SignalRegionTAXMismatched_Fin", "Signal region (TAX), after all cuts",             500, -100., 100., 80, 0., 0.2));

    BookHisto("hSignalRegionTAXAll_In",     new TH2D("SignalRegionTAXAll_In", "Signal region (TAX), before all cuts",             500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAXAll_Track",  new TH2D("SignalRegionTAXAll_Track", "Signal region (TAX), after track-quality cuts", 500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAXAll_Energy", new TH2D("SignalRegionTAXAll_Energy", "Signal region (TAX), after energy cuts",       500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAXAll_Vetoes", new TH2D("SignalRegionTAXAll_Vetoes", "Signal region (TAX), after veto cuts",         500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAXAll_Geom",   new TH2D("SignalRegionTAXAll_Geom", "Signal region (TAX), after geometrical cuts",    500, -100., 100., 80, 0., 0.2));
    BookHisto("hSignalRegionTAXAll_Fin",    new TH2D("SignalRegionTAXAll_Fin", "Signal region (TAX), after all cuts",             500, -100., 100., 80, 0., 0.2));

    // Unique signal region

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

    // CDA vs Z CDA

    BookHisto("hCDAvsZCDATarget_In",     new TH2D("CDAvsZCDATarget_In", "N trajectory wrt target-TAX line, before all cuts",             100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDATarget_Track",  new TH2D("CDAvsZCDATarget_Track", "N trajectory wrt target-TAX line, after track-quality cuts", 100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDATarget_Energy", new TH2D("CDAvsZCDATarget_Energy", "N trajectory wrt target-TAX line, after energy cuts",       100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDATarget_Vetoes", new TH2D("CDAvsZCDATarget_Vetoes", "N trajectory wrt target-TAX line, after veto cuts",         100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDATarget_Geom",   new TH2D("CDAvsZCDATarget_Geom", "N trajectory wrt target-TAX line, after geometrical cuts",    100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDATarget_Fin",    new TH2D("CDAvsZCDATarget_Fin", "N trajectory wrt target-TAX line, after all cuts",             100, -50., 50., 30, 0., 0.15));

    BookHisto("hCDAvsZCDATAX_In",     new TH2D("CDAvsZCDATAX_In", "N trajectory wrt target-TAX line, before all cuts",             100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDATAX_Track",  new TH2D("CDAvsZCDATAX_Track", "N trajectory wrt target-TAX line, after track-quality cuts", 100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDATAX_Energy", new TH2D("CDAvsZCDATAX_Energy", "N trajectory wrt target-TAX line, after energy cuts",       100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDATAX_Vetoes", new TH2D("CDAvsZCDATAX_Vetoes", "N trajectory wrt target-TAX line, after veto cuts",         100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDATAX_Geom",   new TH2D("CDAvsZCDATAX_Geom", "N trajectory wrt target-TAX line, after geometrical cuts",    100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDATAX_Fin",    new TH2D("CDAvsZCDATAX_Fin", "N trajectory wrt target-TAX line, after all cuts",             100, -50., 50., 30, 0., 0.15));

    BookHisto("hCDAvsZCDAAll_In",     new TH2D("CDAvsZCDAAll_In", "N trajectory wrt target-TAX line, before all cuts",             100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDAAll_Track",  new TH2D("CDAvsZCDAAll_Track", "N trajectory wrt target-TAX line, after track-quality cuts", 100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDAAll_Energy", new TH2D("CDAvsZCDAAll_Energy", "N trajectory wrt target-TAX line, after energy cuts",       100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDAAll_Vetoes", new TH2D("CDAvsZCDAAll_Vetoes", "N trajectory wrt target-TAX line, after veto cuts",         100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDAAll_Geom",   new TH2D("CDAvsZCDAAll_Geom", "N trajectory wrt target-TAX line, after geometrical cuts",    100, -50., 50., 30, 0., 0.15));
    BookHisto("hCDAvsZCDAAll_Fin",    new TH2D("CDAvsZCDAAll_Fin", "N trajectory wrt target-TAX line, after all cuts",             100, -50., 50., 30, 0., 0.15));

    // Theta vs Z CDA (studies on SR resolution)

    BookHisto("hThetavsZCDA_Tar", new TH2D("ThetavsZCDA_Tar", "N angular distribution for target-produced events", 200, -15., 15., 50, 0., 0.005));
    BookHisto("hThetavsZCDA_TAX", new TH2D("ThetavsZCDA_TAX", "N angular distribution for TAX-produced events", 200, 0., 50., 50, 0., 0.005));


    // Sidebands 

    BookHisto("hCDA",  new TH1D("CDA", "Vertex-beamline distance for parasitic background", 200, 0., 1000.));
    BookHisto("hTime", new TH1D("Time", "Track time difference for combinatorial background", 500, -15., 15.));
    BookHisto("hZ", new TH1D("Z", "Z of vertex for prompt background", 300, 100., 190.));
    BookHisto("hBeamdistvsMass", new TH2D("BeamdistvsMass", "Vertex-beamline distance vs reconstructed mass", 200, 0., 2., 200, 0., 1000.));

    // Others

    BookHisto("hNMUV3Cand",   new TH1D("NMUV3Cand", "MUV3 candidates for each track", 4, -0.5, 3.5));
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
  Double_t Weight = 1.;

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

  FillHisto("hN2tracks", 0.5);
  FillHisto("hCuts", CutID);
  CutID++;

  // Track features
    
  Int_t Charge1 = Tracks[0].GetCharge();
  Int_t Charge2 = Tracks[1].GetCharge();
  Double_t ChiSquare1 = Tracks[0].GetChi2();
  Double_t ChiSquare2 = Tracks[1].GetChi2();
  TRecoSpectrometerCandidate* SpectrometerCand1 = Tracks[0].GetSpectrometerCandidate();
  TRecoSpectrometerCandidate* SpectrometerCand2 = Tracks[1].GetSpectrometerCandidate();
  TVector3 Mom1 = SpectrometerCand1->GetThreeMomentumBeforeMagnet();
  TVector3 Mom2 = SpectrometerCand2->GetThreeMomentumBeforeMagnet();
  TVector3 TotMom = Mom1 + Mom2;
  Double_t CHODTime1 = Tracks[0].GetCHODTime();
  Double_t CHODTime2 = Tracks[1].GetCHODTime();

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
      //Double_t r = sqrt(x*x + y*y); 
      //Double_t rShifted = sqrt(pow(x-fxStrawChamberCentre[j],2) + y*y); 
      //if (rShifted > frMinStraw && r < frMaxStraw) {
	TString name = Form("hXYSpec%dReco",j);
	FillHisto(name, x/1000., y/1000.);
	//}
    }
    Double_t x = Cand->xAtAfterMagnet(fzCHODPlane);
    Double_t y = Cand->yAtAfterMagnet(fzCHODPlane);
    //Double_t r = sqrt(x*x+y*y);
    //Double_t r1 = frMinCHOD;
    //Double_t r2 = frMaxCHOD;
    //if (r > r1 && r < r2)
	FillHisto("hXYCHODReco", x/1000., y/1000.);
  }

  // CDA of track1 wrt track2

  fCDAcomp = new TwoLinesCDA();
  
  fCDAcomp->SetLine1Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);
  fCDAcomp->SetDir1(Mom1);
  fCDAcomp->SetLine2Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand2->yAtBeforeMagnet(10.0), 10.0);
  fCDAcomp->SetDir2(Mom2);
  fCDAcomp->ComputeVertexCDA();

  Double_t CDA = fCDAcomp->GetCDA();
  TVector3 Vertex = fCDAcomp->GetVertex();
  Double_t Zvertex = fCDAcomp->GetVertex().z();

  // CDA of HNL wrt beamline

  fCDAcomp = new TwoLinesCDA();

  fCDAcomp->SetLine1Point1(0.0, 0.0, 102000.);
  fCDAcomp->SetDir1(0., 0., 1.2E-3);
  fCDAcomp->SetLine2Point1(Vertex);
  fCDAcomp->SetDir2(TotMom);
  fCDAcomp->ComputeVertexCDA();

  Double_t CDAMom = fCDAcomp->GetCDA();

  // CDA of track1,2 wrt beamline

  fCDAcomp = new TwoLinesCDA();

  fCDAcomp->SetLine1Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);
  fCDAcomp->SetDir1(Mom1);
  fCDAcomp->SetLine2Point1(0.0, 0.0, 102000.);
  fCDAcomp->SetDir2(0., 0., 1.2E-3);
  fCDAcomp->ComputeVertexCDA();

  Double_t CDA1 = fCDAcomp->GetCDA();

  fCDAcomp = new TwoLinesCDA();

  fCDAcomp->SetLine1Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);
  fCDAcomp->SetDir1(Mom2);
  fCDAcomp->SetLine2Point1(0.0, 0.0, 102000.);
  fCDAcomp->SetDir2(0., 0., 1.2E-3);
  fCDAcomp->ComputeVertexCDA();

  Double_t CDA2 = fCDAcomp->GetCDA();

  // Distance of HNL wrt target-TAX line

  fCDAcomp = new TwoLinesCDA();

  fCDAcomp->SetLine1Point1(Vertex);
  fCDAcomp->SetDir1(TotMom);
  fCDAcomp->SetLine2Point1(0., 0., -26.5);
  fCDAcomp->SetLine2Point2(0., -22., 23230.);
  fCDAcomp->ComputeVertexCDA();

  Double_t CDALine = fCDAcomp->GetCDA();
  Double_t ZCDALine = fCDAcomp->GetVertex().z();

  // Distance of HNL vertex wrt beamline

  fDistcomp = new PointLineDistance();
  
  fDistcomp->SetLinePoint1(0., 0., 102000.);
  fDistcomp->SetLineDir(0., 0., 1.2E-3);
  fDistcomp->SetPoint(Vertex);
  
  fDistcomp->ComputeDistance();
  
  Double_t BeamlineDist = fDistcomp->GetDistance();

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

  // Distance of HNL wrt target/TAXs

  fDistcomp = new PointLineDistance();
  
  fDistcomp->SetLineDir(TotMom);    
  fDistcomp->SetLinePoint1(Vertex);
  fDistcomp->SetPoint(0., 0., -26.5); // Z = mean of Z of N prod point  
  fDistcomp->ComputeDistance();
  
  Double_t TargetDist = fDistcomp->GetDistance();

  fDistcomp = new PointLineDistance();

  fDistcomp->SetLineDir(TotMom);    
  fDistcomp->SetLinePoint1(Vertex);
  fDistcomp->SetPoint(0., -22., 23230.); // Z = mean of Z of N prod point
  fDistcomp->ComputeDistance();
  
  Double_t TAXDist = fDistcomp->GetDistance();

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

  Double_t xSR = -999.;
  Double_t ySR = -999.;
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

  Bool_t Target = false;
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

  // Reference plot - 1 (* useless)

  Double_t TargetXMin = -15000.;
  Double_t TargetXMax = 40000.;
  Double_t TargetYMin = 0.;
  Double_t TargetYMax = 50.;
  Double_t TAXXMin = -15000.;
  Double_t TAXXMax = 40000.;
  Double_t TAXYMin = 0.;
  Double_t TAXYMax = 50.;
  Double_t xMin = -15000.;
  Double_t xMax = 40000.;
  Double_t yMin = 0.;
  Double_t yMax = 50.;

  FillHisto("hCDAvsZ_In", Zvertex/1000., CDAMom/1000., Weight); //*
  FillHisto("hZvsBeam_In", Zvertex/1000., BeamlineDist/1000., Weight);

  if (GetWithMC()) {
    if (Target == true) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTar_In", ZVertexMomTarget/1000., TargetDist/1000., Weight); //*
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTAXMismatched_In", ZVertexMomTAX/1000., TAXDist/1000., Weight); //*
      }
      FillHisto("hBeamvsTar_In", TargetDist/1000., BeamlineDist/1000., Weight); //*
      FillHisto("hBeamvsTAXMismatched_In", TAXDist/1000., BeamlineDist/1000., Weight); //*
      FillHisto("hCDAvsZCDATarget_In", ZCDALine/1000., CDALine/1000., Weight); //*
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_In", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTarMismatched_In", ZVertexMomTarget/1000., TargetDist/1000., Weight); //*
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTAX_In", ZVertexMomTAX/1000., TAXDist/1000., Weight); //*
      }
      FillHisto("hBeamvsTAX_In", TAXDist/1000., BeamlineDist/1000., Weight); //*
      FillHisto("hBeamvsTarMismatched_In", TargetDist/1000., BeamlineDist/1000., Weight); //*
      FillHisto("hCDAvsZCDATAX_In", ZCDALine/1000., CDALine/1000., Weight); //*
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_In", xSR/1000., ySR/1000., Weight);
      }
    }
  }

  if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
    FillHisto("hSignalRegionTarAll_In", ZVertexMomTarget/1000., TargetDist/1000., Weight); //*
  }
  if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {
    FillHisto("hSignalRegionTAXAll_In", ZVertexMomTAX/1000., TAXDist/1000., Weight); //*
  }

  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_In", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hBeamvsTarAll_In", TargetDist/1000., BeamlineDist/1000., Weight); //*
  FillHisto("hBeamvsTAXAll_In", TAXDist/1000., BeamlineDist/1000., Weight); //*
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

  // Plot the number of candidates associated to each track, for MUV3
  
  FillHisto("hNMUV3Cand", Tracks[0].GetNMUV3AssociationRecords());
  FillHisto("hNMUV3Cand", Tracks[1].GetNMUV3AssociationRecords());
  
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
  Int_t Assoc   = 0;
  Int_t NoAssoc = 0;
  
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

  FillHisto("hCDAvsZ_Track", Zvertex/1000., CDAMom/1000., Weight);
  FillHisto("hZvsBeam_Track", Zvertex/1000., BeamlineDist/1000., Weight);

  if (GetWithMC()) {
    if (Target == true) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTar_Track", ZVertexMomTarget/1000., TargetDist/1000., Weight);
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTAXMismatched_Track", ZVertexMomTAX/1000., TAXDist/1000., Weight);
      }
      FillHisto("hBeamvsTar_Track", TargetDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hBeamvsTAXMismatched_Track", TAXDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hCDAvsZCDATarget_Track", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_Track", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTarMismatched_Track", ZVertexMomTarget/1000., TargetDist/1000., Weight);
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTAX_Track", ZVertexMomTAX/1000., TAXDist/1000., Weight);
      }
      FillHisto("hBeamvsTAX_Track", TAXDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hBeamvsTarMismatched_Track", TargetDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hCDAvsZCDATAX_Track", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_Track", xSR/1000., ySR/1000., Weight);
      }
    }
  }
  
  if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
    FillHisto("hSignalRegionTarAll_Track", ZVertexMomTarget/1000., TargetDist/1000., Weight);
  }
  if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {
    FillHisto("hSignalRegionTAXAll_Track", ZVertexMomTAX/1000., TAXDist/1000., Weight);
  }

  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_Track", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hBeamvsTarAll_Track", TargetDist/1000., BeamlineDist/1000., Weight);
  FillHisto("hBeamvsTAXAll_Track", TAXDist/1000., BeamlineDist/1000., Weight);
  FillHisto("hCDAvsZCDAAll_Track", ZCDALine/1000., CDALine/1000., Weight);

  if (Assoc == 1) 
    FillHisto("hCDAvsCDA_Track", CDA2/1000., CDA1/1000.); //*
  else if (Assoc == 2)
    FillHisto("hCDAvsCDA_Track", CDA1/1000., CDA2/1000.); //*
  /*
  // RICH hand-made mass

  Int_t iRICHCandClose = 999;                                                                     
  Double_t RICHMinDist = 999.;                                                              
  Bool_t MuonRICHOK = false;                                                                  
  Bool_t PionRICHOK = false;                                                             
  Double_t RICHMass = -9999.9;                                                
  Double_t RICHDist = 10.;
  Int_t BookedCand = 999;

  for (Int_t i = 0; i < 2; i++) {

    iRICHCandClose = 999;
    RICHMinDist = 999.;
    RICHMass = -9999.9;

    Double_t FocalLength = 17020.0;
    Double_t RefIndex = 1.0000627022;
    Double_t trackCenter = TMath::Sqrt(FocalLength*Tracks[i].GetSlopeXAfterMagnet()*FocalLength*Tracks[i].GetSlopeXAfterMagnet() + FocalLength*Tracks[i].GetSlopeYAfterMagnet()*FocalLength*Tracks[i].GetSlopeYAfterMagnet());

    if (RICHEvent->GetNRingCandidates() >= 2) {
      for (Int_t iRICH = 0; iRICH < RICHEvent->GetNRingCandidates(); iRICH++) {
	TRecoRICHCandidate* RICHCand =(TRecoRICHCandidate*)RICHEvent->GetRingCandidate(iRICH);
	Double_t RICHRingCenter = TMath::Sqrt(RICHCand->GetRingCenter().X()*RICHCand->GetRingCenter().X() + RICHCand->GetRingCenter().Y()*RICHCand->GetRingCenter().Y());

	if (TMath::Abs(RICHRingCenter - trackCenter) <= RICHMinDist) {
	  RICHMinDist = TMath::Abs(RICHRingCenter - trackCenter);
	  iRICHCandClose = iRICH;
	}
      }

      if (i == 0)
	BookedCand = iRICHCandClose;

      if (RICHMinDist <= RICHDist) {
	if (i == 0 || (i == 1 && BookedCand != iRICHCandClose)) {

	  TRecoRICHCandidate* RICHCandClose = (TRecoRICHCandidate*)RICHEvent->GetRingCandidate(iRICHCandClose);
	  Double_t RingRadius = RICHCandClose->GetRingRadius();
	  Double_t RICHMass2 = -9999.9;
	  Double_t cosTheta = (RingRadius/FocalLength);
	  Double_t n2 = RefIndex*RefIndex;
	  Double_t ptrack = Tracks[i].GetMomentum();
	  
	  RICHMass2 = ptrack*ptrack*((n2/(1+cosTheta*cosTheta)) - 1);
	  
	  if(RICHMass2 > 0) {
	    RICHMass = TMath::Sqrt(RICHMass2);
	  }
	  else {
	    RICHMass = -1.0*TMath::Sqrt(-1.0*RICHMass2);
	  }
	
	  if ((Assoc == 1 && i == 0) || (Assoc == 2 && i == 1))
	    FillHisto("hRICHMass1", RICHMass);
	  else if ((Assoc == 1 && i == 1) || (Assoc == 2 && i == 0))
	    FillHisto("hRICHMass2", RICHMass);
	}
      
	if (Assoc == 1 && i == 0 && RICHMass >= 80. && RICHMass <= 120.)
	  MuonRICHOK = true;
	else if (Assoc == 2 && i == 1 && RICHMass >= 80. && RICHMass <= 120.)
	  MuonRICHOK = true;
	else if (Assoc == 1 && i == 1 && RICHMass >= 120. && RICHMass <= 160.)
	  PionRICHOK = true;
	else if (Assoc == 2 && i == 0 && RICHMass >= 120. && RICHMass <= 160.)
	  PionRICHOK = true;
	else {
	  MuonRICHOK = false;
	  PionRICHOK = false;
	}
      }
    }
  }
    
    if (!MuonRICHOK || !PionRICHOK)
      return;
  */    
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

  if (PiEoP >= 0.8)
    return;

  FillHisto("hCuts", CutID);
  CutID++;

  // Reference plot - 3 

  FillHisto("hCDAvsZ_Energy", Zvertex/1000., CDAMom/1000., Weight);
  FillHisto("hZvsBeam_Energy", Zvertex/1000., BeamlineDist/1000., Weight);

  if (GetWithMC()) {
    if (Target == true) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTar_Energy", ZVertexMomTarget/1000., TargetDist/1000., Weight);
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
	FillHisto("hSignalRegionTAXMismatched_Energy", ZVertexMomTAX/1000., TAXDist/1000., Weight);
      }
      FillHisto("hBeamvsTar_Energy", TargetDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hBeamvsTAXMismatched_Energy", TAXDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hCDAvsZCDATarget_Energy", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_Energy", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTarMismatched_Energy", ZVertexMomTarget/1000., TargetDist/1000., Weight);
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
	FillHisto("hSignalRegionTAX_Energy", ZVertexMomTAX/1000., TAXDist/1000., Weight);
      }
      FillHisto("hBeamvsTAX_Energy", TAXDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hBeamvsTarMismatched_Energy", TargetDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hCDAvsZCDATAX_Energy", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_Energy", xSR/1000., ySR/1000., Weight);
      }
    }
  }

  if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
    FillHisto("hSignalRegionTarAll_Energy", ZVertexMomTarget/1000., TargetDist/1000., Weight);
  }
  if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
    FillHisto("hSignalRegionTAXAll_Energy", ZVertexMomTAX/1000., TAXDist/1000., Weight);
  }

  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_Energy", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hBeamvsTarAll_Energy", TargetDist/1000., BeamlineDist/1000., Weight);
  FillHisto("hBeamvsTAXAll_Energy", TAXDist/1000., BeamlineDist/1000., Weight);
  FillHisto("hCDAvsZCDAAll_Energy", ZCDALine/1000., CDALine/1000., Weight);

  if (Assoc == 1) 
    FillHisto("hCDAvsCDA_Energy", CDA2/1000., CDA1/1000.);
  else if (Assoc == 2)
    FillHisto("hCDAvsCDA_Energy", CDA1/1000., CDA2/1000.);

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

  FillHisto("hCDAvsZ_Vetoes", Zvertex/1000., CDAMom/1000., Weight);
  FillHisto("hZvsBeam_Vetoes", Zvertex/1000., BeamlineDist/1000., Weight);

  if (GetWithMC()) {
    if (Target == true) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTar_Vetoes", ZVertexMomTarget/1000., TargetDist/1000., Weight);
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
	FillHisto("hSignalRegionTAXMismatched_Vetoes", ZVertexMomTAX/1000., TAXDist/1000., Weight);
      }
      FillHisto("hBeamvsTar_Vetoes", TargetDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hBeamvsTAXMismatched_Vetoes", TAXDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hCDAvsZCDATarget_Vetoes", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_Vetoes", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTarMismatched_Vetoes", ZVertexMomTarget/1000., TargetDist/1000., Weight);
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
	FillHisto("hSignalRegionTAX_Vetoes", ZVertexMomTAX/1000., TAXDist/1000., Weight);
      }
      FillHisto("hBeamvsTAX_Vetoes", TAXDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hBeamvsTarMismatched_Vetoes", TargetDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hCDAvsZCDATAX_Vetoes", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_Vetoes", xSR/1000., ySR/1000., Weight);
      }
    }
  }
  
  if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
    FillHisto("hSignalRegionTarAll_Vetoes", ZVertexMomTarget/1000., TargetDist/1000., Weight);
  }
  if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
    FillHisto("hSignalRegionTAXAll_Vetoes", ZVertexMomTAX/1000., TAXDist/1000., Weight);
  }

  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_Vetoes", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hBeamvsTarAll_Vetoes", TargetDist/1000., BeamlineDist/1000., Weight);
  FillHisto("hBeamvsTAXAll_Vetoes", TAXDist/1000., BeamlineDist/1000., Weight);
  FillHisto("hCDAvsZCDAAll_Vetoes", ZCDALine/1000., CDALine/1000., Weight);

  if (Assoc == 1) 
    FillHisto("hCDAvsCDA_Vetoes", CDA2/1000., CDA1/1000.);
  else if (Assoc == 2)
    FillHisto("hCDAvsCDA_Vetoes", CDA1/1000., CDA2/1000.);

  // Geometrical cuts, CUT: Cut on distance between tracks at CH1

  Double_t X = SpectrometerCand1->xAt(fzStraw[0]) - SpectrometerCand2->xAt(fzStraw[0]);
  Double_t Y = SpectrometerCand1->yAt(fzStraw[0]) - SpectrometerCand2->yAt(fzStraw[0]);
  Double_t R = TMath::Sqrt(X*X + Y*Y);
  
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
  /*
  if (BeamlineDist <= 100.)
    return;
  */
  FillHisto("hCuts", CutID);
  CutID++;

  // Reference plot - 5 

  FillHisto("hCDAvsZ_Geom", Zvertex/1000., CDAMom/1000., Weight);
  FillHisto("hZvsBeam_Geom", Zvertex/1000., BeamlineDist/1000., Weight);

  if (GetWithMC()) {
    if (Target == true) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTar_Geom", ZVertexMomTarget/1000., TargetDist/1000., Weight);
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
	FillHisto("hSignalRegionTAXMismatched_Geom", ZVertexMomTAX/1000., TAXDist/1000., Weight);
      }
      FillHisto("hBeamvsTar_Geom", TargetDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hBeamvsTAXMismatched_Geom", TAXDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hCDAvsZCDATarget_Geom", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_Geom", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTarMismatched_Geom", ZVertexMomTarget/1000., TargetDist/1000., Weight);
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
	FillHisto("hSignalRegionTAX_Geom", ZVertexMomTAX/1000., TAXDist/1000., Weight);
      }
      FillHisto("hBeamvsTAX_Geom", TAXDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hBeamvsTarMismatched_Geom", TargetDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hCDAvsZCDATAX_Geom", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_Geom", xSR/1000., ySR/1000., Weight);
      }
    }
  }
  
  if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
    FillHisto("hSignalRegionTarAll_Geom", ZVertexMomTarget/1000., TargetDist/1000., Weight);
  }
  if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
    FillHisto("hSignalRegionTAXAll_Geom", ZVertexMomTAX/1000., TAXDist/1000., Weight);
  }

  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_Geom", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hBeamvsTarAll_Geom", TargetDist/1000., BeamlineDist/1000., Weight);
  FillHisto("hBeamvsTAXAll_Geom", TAXDist/1000., BeamlineDist/1000., Weight);
  FillHisto("hCDAvsZCDAAll_Geom", ZCDALine/1000., CDALine/1000., Weight);

  if (Assoc == 1) 
    FillHisto("hCDAvsCDA_Geom", CDA2/1000., CDA1/1000.);
  else if (Assoc == 2)
    FillHisto("hCDAvsCDA_Geom", CDA1/1000., CDA2/1000.);

  // Reference plot - 6

  FillHisto("hCDAvsZ_Fin", Zvertex/1000., CDAMom/1000., Weight);
  FillHisto("hZvsBeam_Fin", Zvertex/1000., BeamlineDist/1000., Weight);

  if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
    FillHisto("hSignalRegionTarAll_Geom", ZVertexMomTarget/1000., TargetDist/1000., Weight);
  }
  if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
    FillHisto("hSignalRegionTAXAll_Geom", ZVertexMomTAX/1000., TAXDist/1000., Weight);
  }

  if (GetWithMC()) {
    if (Target == true) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTar_Fin", ZVertexMomTarget/1000., TargetDist/1000., Weight);
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
	FillHisto("hSignalRegionTAXMismatched_Fin", ZVertexMomTAX/1000., TAXDist/1000., Weight);
      }
      FillHisto("hBeamvsTar_Fin", TargetDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hBeamvsTAXMismatched_Fin", TAXDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hCDAvsZCDATarget_Fin", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTar_Fin", xSR/1000., ySR/1000., Weight);
      }
    }
    else if (Target == false) {
      if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
	FillHisto("hSignalRegionTarMismatched_Fin", ZVertexMomTarget/1000., TargetDist/1000., Weight);
      }
      if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
	FillHisto("hSignalRegionTAX_Fin", ZVertexMomTAX/1000., TAXDist/1000., Weight);
      }
      FillHisto("hBeamvsTAX_Fin", TAXDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hBeamvsTarMismatched_Fin", TargetDist/1000., BeamlineDist/1000., Weight);
      FillHisto("hCDAvsZCDATAX_Fin", ZCDALine/1000., CDALine/1000., Weight);
      if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
	FillHisto("hSRTAX_Fin", xSR/1000., ySR/1000., Weight);
      }
    }
  }
  
  if ((fBlindRegion && (ZVertexMomTarget <= TargetXMin || ZVertexMomTarget >= TargetXMax) && (TargetDist <= TargetYMin || TargetDist >= TargetYMax)) || !fBlindRegion) {
    FillHisto("hSignalRegionTarAll_Fin", ZVertexMomTarget/1000., TargetDist/1000., Weight);
  }
  if ((fBlindRegion && (ZVertexMomTAX <= TAXXMin || ZVertexMomTAX >= TAXXMax) && (TAXDist <= TAXYMin || TAXDist >= TAXYMax)) || !fBlindRegion) {  
    FillHisto("hSignalRegionTAXAll_Fin", ZVertexMomTAX/1000., TAXDist/1000., Weight);
  }

  if ((fBlindRegion && (xSR <= xMin || xSR >= xMax) && (ySR <= yMin || ySR >= yMax)) || !fBlindRegion) {
    FillHisto("hSR_Fin", xSR/1000., ySR/1000., Weight);
  }

  FillHisto("hBeamvsTarAll_Fin", TargetDist/1000., BeamlineDist/1000., Weight);
  FillHisto("hBeamvsTAXAll_Fin", TAXDist/1000., BeamlineDist/1000., Weight);
  FillHisto("hCDAvsZCDAAll_Fin", ZCDALine/1000., CDALine/1000., Weight);

  if (Assoc == 1) 
    FillHisto("hCDAvsCDA_Fin", CDA2/1000., CDA1/1000.);
  else if (Assoc == 2)
    FillHisto("hCDAvsCDA_Fin", CDA1/1000., CDA2/1000.);

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

  if (TMath::Abs(invMass - fMassForReco*1000.) <= 10.)
    FillHisto("hInvMassReco", invMass/1000.);

  FillHisto("hBeamdistvsMass", invMass/1000., BeamlineDist, Weight);

  // Output of selection

  fPassSelection = true;

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

  if (!fReadingData) return;

  // X axis

  fHisto.GetTH1("hNk3pi")->GetXaxis()->SetTitle("Number of k3pi");
  fHisto.GetTH1("hNbursts")->GetXaxis()->SetTitle("Number of bursts");
  fHisto.GetTH1("hNEvents")->GetXaxis()->SetTitle("Number of events");
  fHisto.GetTH1("hN2tracks")->GetXaxis()->SetTitle("Number of two-track events");
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

  fHisto.GetTH2("hCDAvsZ_In")->GetXaxis()->SetTitle("Vertex position along Z [m]");
  fHisto.GetTH2("hCDAvsZ_Track")->GetXaxis()->SetTitle("Vertex position along Z [m]");
  fHisto.GetTH2("hCDAvsZ_Energy")->GetXaxis()->SetTitle("Vertex position along Z [m]");
  fHisto.GetTH2("hCDAvsZ_Vetoes")->GetXaxis()->SetTitle("Vertex position along Z [m]");
  fHisto.GetTH2("hCDAvsZ_Geom")->GetXaxis()->SetTitle("Vertex position along Z [m]");
  fHisto.GetTH2("hCDAvsZ_Fin")->GetXaxis()->SetTitle("Vertex position along Z [m]");

  fHisto.GetTH2("hCDAvsCDA_Track")->GetXaxis()->SetTitle("Pion CDA [m]");
  fHisto.GetTH2("hCDAvsCDA_Energy")->GetXaxis()->SetTitle("Pion CDA [m]");
  fHisto.GetTH2("hCDAvsCDA_Vetoes")->GetXaxis()->SetTitle("Pion CDA [m]");
  fHisto.GetTH2("hCDAvsCDA_Geom")->GetXaxis()->SetTitle("Pion CDA [m]");
  fHisto.GetTH2("hCDAvsCDA_Fin")->GetXaxis()->SetTitle("Pion CDA [m]");

  fHisto.GetTH2("hZvsBeam_In")->GetXaxis()->SetTitle("Vertex position along Z [m]");
  fHisto.GetTH2("hZvsBeam_Track")->GetXaxis()->SetTitle("Vertex position along Z [m]");
  fHisto.GetTH2("hZvsBeam_Energy")->GetXaxis()->SetTitle("Vertex position along Z [m]");
  fHisto.GetTH2("hZvsBeam_Vetoes")->GetXaxis()->SetTitle("Vertex position along Z [m]");
  fHisto.GetTH2("hZvsBeam_Geom")->GetXaxis()->SetTitle("Vertex position along Z [m]");
  fHisto.GetTH2("hZvsBeam_Fin")->GetXaxis()->SetTitle("Vertex position along Z [m]");

  fHisto.GetTH2("hBeamvsTar_In")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTar_Track")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTar_Energy")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTar_Vetoes")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTar_Geom")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTar_Fin")->GetXaxis()->SetTitle("Impact parameter [m]");

  fHisto.GetTH2("hBeamvsTarMismatched_In")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTarMismatched_Track")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTarMismatched_Energy")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTarMismatched_Vetoes")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTarMismatched_Geom")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTarMismatched_Fin")->GetXaxis()->SetTitle("Impact parameter [m]");

  fHisto.GetTH2("hBeamvsTarAll_In")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTarAll_Track")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTarAll_Energy")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTarAll_Vetoes")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTarAll_Geom")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTarAll_Fin")->GetXaxis()->SetTitle("Impact parameter [m]");

  fHisto.GetTH2("hBeamvsTAX_In")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAX_Track")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAX_Energy")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAX_Vetoes")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAX_Geom")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAX_Fin")->GetXaxis()->SetTitle("Impact parameter [m]");

  fHisto.GetTH2("hBeamvsTAXMismatched_In")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAXMismatched_Track")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAXMismatched_Energy")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAXMismatched_Vetoes")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAXMismatched_Geom")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAXMismatched_Fin")->GetXaxis()->SetTitle("Impact parameter [m]");

  fHisto.GetTH2("hBeamvsTAXAll_In")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAXAll_Track")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAXAll_Energy")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAXAll_Vetoes")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAXAll_Geom")->GetXaxis()->SetTitle("Impact parameter [m]");
  fHisto.GetTH2("hBeamvsTAXAll_Fin")->GetXaxis()->SetTitle("Impact parameter [m]");

  fHisto.GetTH2("hSignalRegionTar_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTar_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTar_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTar_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTar_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTar_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");

  fHisto.GetTH2("hSignalRegionTarMismatched_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTarMismatched_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTarMismatched_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTarMismatched_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTarMismatched_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTarMismatched_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");

  fHisto.GetTH2("hSignalRegionTarAll_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTarAll_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTarAll_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTarAll_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTarAll_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");
  fHisto.GetTH2("hSignalRegionTarAll_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt target line [m]");

  fHisto.GetTH2("hSignalRegionTAX_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAX_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAX_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAX_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAX_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAX_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");

  fHisto.GetTH2("hSignalRegionTAXMismatched_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");

  fHisto.GetTH2("hSignalRegionTAXAll_In")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAXAll_Track")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAXAll_Energy")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAXAll_Vetoes")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAXAll_Geom")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");
  fHisto.GetTH2("hSignalRegionTAXAll_Fin")->GetXaxis()->SetTitle("Z of CDA of mother wrt TAX line [m]");

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

  fHisto.GetTH1("hNMUV3Cand")->GetXaxis()->SetTitle("Number of candidates");
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

  fHisto.GetTH2("hCDAvsZ_In")->GetYaxis()->SetTitle("CDA [m]");
  fHisto.GetTH2("hCDAvsZ_Track")->GetYaxis()->SetTitle("CDA [m]");
  fHisto.GetTH2("hCDAvsZ_Energy")->GetYaxis()->SetTitle("CDA [m]");
  fHisto.GetTH2("hCDAvsZ_Vetoes")->GetYaxis()->SetTitle("CDA [m]");
  fHisto.GetTH2("hCDAvsZ_Geom")->GetYaxis()->SetTitle("CDA [m]");
  fHisto.GetTH2("hCDAvsZ_Fin")->GetYaxis()->SetTitle("CDA [m]");

  fHisto.GetTH2("hCDAvsCDA_Track")->GetYaxis()->SetTitle("Muon CDA [m]");
  fHisto.GetTH2("hCDAvsCDA_Energy")->GetYaxis()->SetTitle("Muon CDA [m]");
  fHisto.GetTH2("hCDAvsCDA_Vetoes")->GetYaxis()->SetTitle("Muon CDA [m]");
  fHisto.GetTH2("hCDAvsCDA_Geom")->GetYaxis()->SetTitle("Muon CDA [m]");
  fHisto.GetTH2("hCDAvsCDA_Fin")->GetYaxis()->SetTitle("Muon CDA [m]");

  fHisto.GetTH2("hZvsBeam_In")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hZvsBeam_Track")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hZvsBeam_Energy")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hZvsBeam_Vetoes")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hZvsBeam_Geom")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hZvsBeam_Fin")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");

  fHisto.GetTH2("hBeamvsTar_In")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTar_Track")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTar_Energy")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTar_Vetoes")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTar_Geom")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTar_Fin")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");

  fHisto.GetTH2("hBeamvsTarMismatched_In")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTarMismatched_Track")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTarMismatched_Energy")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTarMismatched_Vetoes")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTarMismatched_Geom")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTarMismatched_Fin")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");

  fHisto.GetTH2("hBeamvsTarAll_In")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTarAll_Track")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTarAll_Energy")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTarAll_Vetoes")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTarAll_Geom")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTarAll_Fin")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");

  fHisto.GetTH2("hBeamvsTAX_In")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAX_Track")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAX_Energy")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAX_Vetoes")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAX_Geom")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAX_Fin")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");

  fHisto.GetTH2("hBeamvsTAXMismatched_In")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAXMismatched_Track")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAXMismatched_Energy")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAXMismatched_Vetoes")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAXMismatched_Geom")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAXMismatched_Fin")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");

  fHisto.GetTH2("hBeamvsTAXAll_In")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAXAll_Track")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAXAll_Energy")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAXAll_Vetoes")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAXAll_Geom")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
  fHisto.GetTH2("hBeamvsTAXAll_Fin")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");

  fHisto.GetTH2("hSignalRegionTar_In")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTar_Track")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTar_Energy")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTar_Vetoes")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTar_Geom")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTar_Fin")->GetYaxis()->SetTitle("Impact parameter of mother [m]");

  fHisto.GetTH2("hSignalRegionTarMismatched_In")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTarMismatched_Track")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTarMismatched_Energy")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTarMismatched_Vetoes")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTarMismatched_Geom")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTarMismatched_Fin")->GetYaxis()->SetTitle("Impact parameter of mother [m]");

  fHisto.GetTH2("hSignalRegionTarAll_In")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTarAll_Track")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTarAll_Energy")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTarAll_Vetoes")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTarAll_Geom")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTarAll_Fin")->GetYaxis()->SetTitle("Impact parameter of mother [m]");

  fHisto.GetTH2("hSignalRegionTAX_In")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAX_Track")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAX_Energy")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAX_Vetoes")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAX_Geom")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAX_Fin")->GetYaxis()->SetTitle("Impact parameter of mother [m]");

  fHisto.GetTH2("hSignalRegionTAXMismatched_In")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Track")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Energy")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Vetoes")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Geom")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Fin")->GetYaxis()->SetTitle("Impact parameter of mother [m]");

  fHisto.GetTH2("hSignalRegionTAXAll_In")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAXAll_Track")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAXAll_Energy")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAXAll_Vetoes")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAXAll_Geom")->GetYaxis()->SetTitle("Impact parameter of mother [m]");
  fHisto.GetTH2("hSignalRegionTAXAll_Fin")->GetYaxis()->SetTitle("Impact parameter of mother [m]");

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

  fHisto.GetTH2("hCDAvsZ_In")->GetZaxis()->SetTitle("Normalized to POT and 100 cm^{2}");
  fHisto.GetTH2("hCDAvsZ_Track")->GetZaxis()->SetTitle("Normalized to POT and 100 cm^{2}");
  fHisto.GetTH2("hCDAvsZ_Energy")->GetZaxis()->SetTitle("Normalized to POT and 100 cm^{2}");
  fHisto.GetTH2("hCDAvsZ_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 100 cm^{2}");
  fHisto.GetTH2("hCDAvsZ_Geom")->GetZaxis()->SetTitle("Normalized to POT and 100 cm^{2}");
  fHisto.GetTH2("hCDAvsZ_Fin")->GetZaxis()->SetTitle("Normalized to POT and 100 cm^{2}");

  fHisto.GetTH2("hCDAvsCDA_Track")->GetZaxis()->SetTitle("Normalized to POT and 1 cm^{2}");
  fHisto.GetTH2("hCDAvsCDA_Energy")->GetZaxis()->SetTitle("Normalized to POT and 1 cm^{2}");
  fHisto.GetTH2("hCDAvsCDA_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 1 cm^{2}");
  fHisto.GetTH2("hCDAvsCDA_Geom")->GetZaxis()->SetTitle("Normalized to POT and 1 cm^{2}");
  fHisto.GetTH2("hCDAvsCDA_Fin")->GetZaxis()->SetTitle("Normalized to POT and 1 cm^{2}");

  fHisto.GetTH2("hZvsBeam_In")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");
  fHisto.GetTH2("hZvsBeam_Track")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");
  fHisto.GetTH2("hZvsBeam_Energy")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");
  fHisto.GetTH2("hZvsBeam_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");
  fHisto.GetTH2("hZvsBeam_Geom")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");
  fHisto.GetTH2("hZvsBeam_Fin")->GetZaxis()->SetTitle("Normalized to POT and 400 cm^{2}");

  fHisto.GetTH2("hBeamvsTar_In")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTar_Track")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTar_Energy")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTar_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTar_Geom")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTar_Fin")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");

  fHisto.GetTH2("hBeamvsTarMismatched_In")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTarMismatched_Track")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTarMismatched_Energy")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTarMismatched_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTarMismatched_Geom")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTarMismatched_Fin")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");

  fHisto.GetTH2("hBeamvsTarAll_In")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTarAll_Track")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTarAll_Energy")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTarAll_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTarAll_Geom")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTarAll_Fin")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");

  fHisto.GetTH2("hBeamvsTAX_In")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAX_Track")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAX_Energy")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAX_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAX_Geom")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAX_Fin")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");

  fHisto.GetTH2("hBeamvsTAXMismatched_In")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAXMismatched_Track")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAXMismatched_Energy")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAXMismatched_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAXMismatched_Geom")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAXMismatched_Fin")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");

  fHisto.GetTH2("hBeamvsTAXAll_In")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAXAll_Track")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAXAll_Energy")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAXAll_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAXAll_Geom")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");
  fHisto.GetTH2("hBeamvsTAXAll_Fin")->GetZaxis()->SetTitle("Normalized to POT and 2 cm^{2}");

  fHisto.GetTH2("hSignalRegionTar_In")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTar_Track")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTar_Energy")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTar_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTar_Geom")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTar_Fin")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");

  fHisto.GetTH2("hSignalRegionTarMismatched_In")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTarMismatched_Track")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTarMismatched_Energy")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTarMismatched_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTarMismatched_Geom")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTarMismatched_Fin")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");

  fHisto.GetTH2("hSignalRegionTarAll_In")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTarAll_Track")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTarAll_Energy")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTarAll_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTarAll_Geom")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTarAll_Fin")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");

  fHisto.GetTH2("hSignalRegionTAX_In")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAX_Track")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAX_Energy")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAX_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAX_Geom")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAX_Fin")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");

  fHisto.GetTH2("hSignalRegionTAXMismatched_In")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Track")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Energy")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Geom")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAXMismatched_Fin")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");

  fHisto.GetTH2("hSignalRegionTAXAll_In")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAXAll_Track")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAXAll_Energy")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAXAll_Vetoes")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAXAll_Geom")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");
  fHisto.GetTH2("hSignalRegionTAXAll_Fin")->GetZaxis()->SetTitle("Normalized to POT and 10 cm^{2}");

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
  
  SaveAllPlots();

  return;
}
