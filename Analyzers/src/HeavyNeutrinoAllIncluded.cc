// ---------------------------------------------------------------                                      
//                                                                                                      
// History:                                                                                             
//                                                                                                      
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018                               
//                                                                                                      
// ---------------------------------------------------------------                                      
/// \class HeavyNeutrinoAllInclusive
/// \Brief                                                                                              
/// Scan on the mass and coupling of the heavy neutral leptons. It includes selection, weight and yield 
/// computation as a function of the HNL mass                                                           
/// \EndBrief                                                                                           
/// \Detailed                                                                                           
/// If the analyzer is run on MC samples, it produces a weight to be associated                         
/// to each heavy neutral lepton in the sample. After applying selection cuts,                          
/// acceptance and yield per POT are also computed. Moreover, an expected exclusion plot is produced as 
/// a function of the HNL mass and its coupling to SM leptons. Several quantities are also plotted as 
/// a function of either the HNL mass or its coupling.
/// All couplings are expressed as Log(U2).                                                            
/// Thus, the analyzer is able to run either on MC samples of just one HNL mass or on samples
/// containing HNLs of different masses.                                                                
/// If the analyzer is run on data samples, only the selection is applied.                              
/// This analyzer makes use of an external library implemented in PhysicsObjects, called HNLFunctions.  
/// The values of the ratios between specific-flavour couplings can be either set as external           
/// parameters from command line or taken as default values.                                            
/// For example, if the user sets UeRatio = 5., UmuRatio = 1., UtauRatio = 3.5,                         
/// the specific-flavour coupling values will be: Ue2 = 5.25E-11, Umu2 = 1.05E-11, Utau2 = 3.68E-11.   
/// The values of the initial and final coupling and the scan step can be either set as external      
/// parameters from command line or taken as default values.                                           
/// For example, if the user assigns -3. to the starting coupling, -2. to the final one                
/// and 0.5 to the step, the scan will be performed for Log(U2) = [-3., -2.5, -2.], this meaning       
/// U2 = [1.E-3., 1.E-2.5, 1.E-2.].                                                                    
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
#include "HeavyNeutrinoAllInclusive.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "GeometricAcceptance.hh"
#include "TriggerConditions.hh"
#include "TF1.h"
#include "TF2.h"
#include "HNLFunctions.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

#define TRIGGER_L0_PHYSICS_TYPE 1
#define MCTriggerMask 0xFF
#define LabelSize 0.05

/// \class HeavyNeutrinoAllInclusive

HeavyNeutrinoAllInclusive::HeavyNeutrinoAllInclusive(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrinoAllInclusive") {

  if (!GetIsTree()) return;

  RequestAllMCTrees();
  RequestAllRecoTrees();
  RequestL0Data();
  RequestL1Data();

  AddParam("USquared", &fUSquared, 1.E-6);
  AddParam("UeSquaredRatio", &fUeSquaredRatio, 1.);
  AddParam("UmuSquaredRatio", &fUmuSquaredRatio, 16.);
  AddParam("UtauSquaredRatio", &fUtauSquaredRatio, 3.8);
  AddParam("CouplingStart", &fCouplingStart, -10.);
  AddParam("CouplingStop", &fCouplingStop, 0.);
  AddParam("CouplingStep", &fCouplingStep, 0.1);

  fN = round((std::abs(fCouplingStop-fCouplingStart))/fCouplingStep);
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

  // Selection histos

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

  // Scan histos

  fhReachCoupling     = nullptr;
  fhDecayCoupling     = nullptr;
  fhWeightCoupling    = nullptr;
  fgAccCoupling       = nullptr;
  fgYieldCoupling     = nullptr;
  fgGammaTotCoupling  = nullptr;
  fgTauCoupling       = nullptr;

  fhReachMass     = nullptr;
  fhDecayMass     = nullptr;
  fhWeightMass    = nullptr;
  fgAccMass       = nullptr;
  fgYieldMass     = nullptr;
  fgGammaTotMass  = nullptr;
  fgTauMass       = nullptr;

  fgExclusion = nullptr;
}

void HeavyNeutrinoAllInclusive::InitHist() {

  BookHisto("Selection/hNk3pi",    new TH1D("Nk3pi",    "Total number of K3pi events",       1, 0., 1.));
  BookHisto("Selection/hNbursts",  new TH1D("Nbursts",  "Total number of processed bursts",  1, 0., 1.));
  BookHisto("Selection/hNEvents",  new TH1D("NEvents",  "Number of total processed events" , 1, 0., 1.));
  BookHisto("Selection/hNtracks",  new TH1D("Ntracks",  "Number of tracks",                  4, -0.5, 3.5));
  BookHisto("Selection/hN2tracks", new TH1D("N2tracks", "Number of two-tracks events",       1, 0., 1.));

  BookHisto("Selection/hZDProd",        new TH1D("ZDProd", "Z of D meson production point", 20000., -250., 33000.));
  BookHisto("Selection/hZDDecay",       new TH1D("ZDDecay", "Z of D meson decay point",     20000., -250., 33000.));
  BookHisto("Selection/hDTheta",        new TH1D("DTheta",     "D meson theta",              100,  0., 0.3));
  BookHisto("Selection/hDLambda",       new TH1D("DLambda",    "D meson decay length",       100, -1., 40.));
  BookHisto("Selection/hDPath",         new TH1D("DPath",      "D meson path in Z",          100, -1., 50.));
  BookHisto("Selection/hDMom",          new TH1D("DMom",       "D meson momentum",           100, -1., 170.));

  BookHisto("Selection/hZHNLDecay",     new TH1D("ZHNLDecay",    "Z of HNL decay point",           100., 90., 190.));
  BookHisto("Selection/hHNLGamma",      new TH1D("HNLGamma",     "Lorentz gamma of HNL",           50., 0., 170.));
  BookHisto("Selection/hHNLDecayProb",  new TH1D("HNLDecayProb", "HNL decay probability",          100., 0., 0.0065));
  BookHisto("Selection/hHNLReachProb",  new TH1D("HNLReachProb", "HNL probability of reaching FV", 100., 0.99, 1.001));
  BookHisto("Selection/hHNLTheta",      new TH1D("HNLTheta",     "HNL theta",                      100., 0., 0.5));
  BookHisto("Selection/hHNLMom",        new TH1D("HNLMom",       "HNL momentum",                   100., -0.5, 200.));
  BookHisto("Selection/hWeight",        new TH1D("Weight",       "Weight",                         1000, 1.E-15, 1.E-12));
  BookHisto("Selection/hMomPi",         new TH1D("MomPi",        "Pion momentum",                  100, -0.5, 200.));
  BookHisto("Selection/hMomMu",         new TH1D("MomMu",        "Muon momentum",                  100, -0.5, 200.));

  BookHisto("Selection/hXYSpec0Reco", new TH2D("XYSpec0Reco",        "Two-track reconstructed events at CH1",  100, -1.5, 1.5, 100, -1.5, 1.5)); 
  BookHisto("Selection/hXYSpec1Reco", new TH2D("XYSpec1Reco",        "Two-track reconstructed events at CH2",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("Selection/hXYSpec2Reco", new TH2D("XYSpec2Reco",        "Two-track reconstructed events at CH3",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("Selection/hXYSpec3Reco", new TH2D("XYSpec3Reco",        "Two-track reconstructed events at CH4",  100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("Selection/hXYCHODReco",  new TH2D("XYCHODReco",         "Two-track reconstructed events at CHOD", 100, -1.5, 1.5, 100, -1.5, 1.5));
  BookHisto("Selection/hXYCHODTrue",  new TH2D("XYCHODTrue",         "X,Y of HNL daughters at CHOD, from MC",  100, -2., 2., 100, -2., 2.));
  BookHisto("Selection/hXYMUV3True",  new TH2D("XYMUV3True",         "X,Y of HNL daughters at MUV3, from MC",  100, -2., 2., 100, -2., 2.));
  BookHisto("Selection/hP1vsP2",      new TH2D("P1vsP2",      "Trimomentum of the two HNL daughters, from MC", 100, 0., 100., 100, 0., 100.));

  BookHisto("Selection/hPhysicsEventsVsCuts", new TH1D("PhysicsEventsVsCuts", "Physics events passing the selection cuts", 35, 0., 35.));

  BookHisto("Selection/hCDAvsZVertex_TotMomToBeamlineInitial",              new TH2D("CDAvsZVertex_TotMomToBeamlineInitial",              "Two-track total momentum wrt beam axis, before all cuts",            100, 0., 240., 100, 0., 0.5));
  BookHisto("Selection/hCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack", new TH2D("CDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack", "Two-track total momentum wrt beam axis, after downstream selection", 100, 0., 240., 100, 0., 0.5));
  BookHisto("Selection/hCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts",      new TH2D("CDAvsZVertex_TotMomToBeamlineAfterEnergyCuts",      "Two-track total momentum wrt beam axis, after energy cuts",          100, 0., 240., 100, 0., 0.5));
  BookHisto("Selection/hCDAvsZVertex_TotMomToBeamlineAfterVetoes",          new TH2D("CDAvsZVertex_TotMomToBeamlineAfterVetoes",          "Two-track total momentum wrt beam axis, after vetoes",               100, 0., 240., 100, 0., 0.5));
  BookHisto("Selection/hCDAvsZVertex_TotMomToBeamlineAfterGeomCuts",        new TH2D("CDAvsZVertex_TotMomToBeamlineAfterGeomCuts",        "Two-track total momentum wrt beam axis, after geometrical cuts",     100, 0., 240., 100, 0., 0.5));
  BookHisto("Selection/hCDAvsZVertex_TotMomToBeamlineFinal",                new TH2D("CDAvsZVertex_TotMomToBeamlineFinal",                "Two-track total momentum wrt beam axis, after all selections",       100, 0., 240., 100, 0., 0.5));
  
  BookHisto("Selection/hZvertexvsBeamlineDistInitial",              new TH2D("ZvertexvsBeamlineDistInitial",              "Two track vertex wrt beam axis, before all cuts",            50, 0., 240., 100, 0., 1.));
  BookHisto("Selection/hZvertexvsBeamlineDistAfterDownstreamTrack", new TH2D("ZvertexvsBeamlineDistAfterDownstreamTrack", "Two track vertex wrt beam axis, after downstream selection", 50, 0., 240., 100, 0., 1.));
  BookHisto("Selection/hZvertexvsBeamlineDistAfterEnergyCuts",      new TH2D("ZvertexvsBeamlineDistAfterEnergyCuts",      "Two track vertex wrt beam axis, after energy cuts",          50, 0., 240., 100, 0., 1.));
  BookHisto("Selection/hZvertexvsBeamlineDistAfterVetoes",          new TH2D("ZvertexvsBeamlineDistAfterVetoes",          "Two track vertex wrt beam axis, after vetoes",               50, 0., 240., 100, 0., 1.));
  BookHisto("Selection/hZvertexvsBeamlineDistAfterGeomCuts",        new TH2D("ZvertexvsBeamlineDistAfterGeomCuts",        "Two track vertex wrt beam axis, after geometrical cuts",     50, 0., 240., 100, 0., 1.));
  BookHisto("Selection/hZvertexvsBeamlineDistFinal",                new TH2D("ZvertexvsBeamlineDistFinal",                "Two track vertex wrt beam axis, after all selections",       50, 0., 240., 100, 0., 1.));

  BookHisto("Selection/hCDAvsZVertex_TrackToBeamlineInitial",  new TH2D("CDAvsZVertex_TrackToBeamlineInitial",  "Track wrt beam axis, before all cuts",      100, 0., 240., 100, 0., 0.2));
  BookHisto("Selection/hCDAvsZVertex_TrackToBeamlineAfterCut", new TH2D("CDAvsZVertex_TrackToBeamlineAfterCut", "Track wrt beam axis, after cut",            100, 0., 240., 100, 0., 1.));
  BookHisto("Selection/hCDAvsZVertex_TrackToBeamlineFinal",    new TH2D("CDAvsZVertex_TrackToBeamlineFinal",    "Track wrt beam axis, after all selections", 100, 0., 240., 100, 0., 1.));
  BookHisto("Selection/hCDAvsZVertex_TrackToTrackInitial",     new TH2D("CDAvsZVertex_TrackToTrackInitial",     "Track1 wrt track2, before all cuts",        100, 0., 240., 100, 0., 0.2));
  BookHisto("Selection/hCDAvsZVertex_TrackToTrackAfterCut",    new TH2D("CDAvsZVertex_TrackToTrackAfterCut",    "Track1 wrt track2, after cut",              100, 0., 240., 100, 0., 1.));
  BookHisto("Selection/hCDAvsZVertex_TrackToTrackFinal",       new TH2D("CDAvsZVertex_TrackToTrackFinal",       "Track1 wrt track2, after all selections",   100, 0., 240., 100, 0., 1.));
  
  BookHisto("Selection/hBeamlineDistvsTargetDist_TotMomInitial",              new TH2D("BeamlineDistvsTargetDist_TotMomInitial",              "Two-track total momentum wrt beam axis, before all cuts",            100, 0., 1., 100, 0., 1.));
  BookHisto("Selection/hBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack", new TH2D("BeamlineDistvsTargetDist_TotMomAfterDownstreamTrack", "Two-track total momentum wrt beam axis, after downstream selection", 100, 0., 1., 100, 0., 1.));
  BookHisto("Selection/hBeamlineDistvsTargetDist_TotMomAfterEnergyCuts",      new TH2D("BeamlineDistvsTargetDist_TotMomAfterEnergyCuts",      "Two-track total momentum wrt beam axis, after energy cuts",          100, 0., 1., 100, 0., 1.));
  BookHisto("Selection/hBeamlineDistvsTargetDist_TotMomAfterVetoes",          new TH2D("BeamlineDistvsTargetDist_TotMomAfterVetoes",          "Two-track total momentum wrt beam axis, after vetoes",               100, 0., 1., 100, 0., 1.));
  BookHisto("Selection/hBeamlineDistvsTargetDist_TotMomAfterGeomCuts",        new TH2D("BeamlineDistvsTargetDist_TotMomAfterGeomCuts",        "Two-track total momentum wrt beam axis, after geometrical cuts",     100, 0., 1., 100, 0., 1.));
  BookHisto("Selection/hBeamlineDistvsTargetDist_TotMomFinal",                new TH2D("BeamlineDistvsTargetDist_TotMomFinal",                "Two-track total momentum wrt beam axis, after all selections",       100, 0., 1., 100, 0., 1.));
  
  BookHisto("Selection/hDeltaTimeFromCHOD",     new TH1D("DeltaTimeFromCHOD",     "Time difference of two tracks (CHOD candidates)",    200, -20., 20.));  
  BookHisto("Selection/hNMUV3CandAssocToTrack", new TH1D("NMUV3CandAssocToTrack", "Number of MUV3 candidates associated to each track", 4, -0.5, 3.5));
  BookHisto("Selection/hNCHODCandAssocToTrack", new TH1D("NCHODCandAssocToTrack", "Number of CHOD candidates associated to each track", 10, 0., 10.));
  
  BookHisto("Selection/hEoP", new TH1D("EoP", "E/p in LKr", 100, 0., 1.2));
  BookHisto("Selection/hEoPMuVsPi", new TH2D("EoPMuVsPi", "Muon E/p vs pion E/p in LKr", 100, 0., 1.2, 100, 0., 0.1));
  
  BookHisto("Selection/hSingleAddEnLKrHit",  new TH2D("SingleAddEnLKrHit", "Single hit energy vs time, for additional LKr hits",              250, -100., 100, 250, 0., 100.));
  BookHisto("Selection/hSingleAddEnLKrCand", new TH2D("SingleAddEnLKrCand", "Single candidate energy vs time, for additional LKr candidates", 250, -100., 100, 250, 0., 100.));
  BookHisto("Selection/hAddEnLKrHit",        new TH2D("AddEnLKrHit", "Additional energy vs time, for LKr hits",                               250, -100., 100, 250, 0., 10.));
  BookHisto("Selection/hAddEnLKrCand",       new TH2D("AddEnLKrCand", "Additional energy vs time, for LKr candidates",                        250, -100., 100, 250, 0., 10.));
  
  BookHisto("Selection/hInvMassMC",   new TH1D("InvMassMC",   "Invariant mass MC", 100, 999., 1001.));
  BookHisto("Selection/hInvMassReco", new TH1D("InvMassReco", "Invariant mass Reco", 50, 960., 1040.));

  BookHisto("Selection/hAcc",   new TH1D("Acc",   "Acceptance", 50, 1.E-6, 5.E-5));
  BookHisto("Selection/hYield", new TH1D("Yield", "Yield",      50, 1.E-17, 1.E-15));

  BookHisto("CouplingScan/hReachCoupling",    new TH2D("ReachCoupling", "Probability of N reaching the FV vs coupling",    fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("CouplingScan/hDecayCoupling",    new TH2D("DecayCoupling", "Probability of N decaying in the FV vs coupling", fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("CouplingScan/hWeightCoupling",   new TH2D("WeightCoupling", "N weight vs coupling",                           fN, fCouplingStart, fCouplingStop, 1000, 1.E-12, 1.E-9));

  fgGammaTot = new TGraph();
  fgGammaTot->SetNameTitle("CouplingScan/GammaTotCoupling", "N total decay width vs coupling");
  BookHisto(fgGammaTot);

  fgTau = new TGraph();
  fgTau->SetNameTitle("CouplingScan/TauCoupling", "N lifetime vs coupling");
  BookHisto(fgTau);

  fgAcc = new TGraph();
  fgAcc->SetNameTitle("CouplingScan/AccCoupling", "Acceptance vs coupling");
  BookHisto(fgAcc);

  fgYield = new TGraph();
  fgYield->SetNameTitle("CouplingScan/YieldCoupling", "Yield per POT vs coupling");
  BookHisto(fgYield);

  BookHisto("MassScan/hReachMass",    new TH2D("ReachMass", "Probability of N reaching the FV vs coupling",    fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("MassScan/hDecayMass",    new TH2D("DecayMass", "Probability of N decaying in the FV vs coupling", fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("MassScan/hWeightMass",   new TH2D("WeightMass", "N weight vs coupling",                           fN, fCouplingStart, fCouplingStop, 1000, 1.E-12, 1.E-9));

  fgGammaTot = new TGraph();
  fgGammaTot->SetNameTitle("MassScan/GammaTotMass", "N total decay width vs coupling");
  BookHisto(fgGammaTot);

  fgTau = new TGraph();
  fgTau->SetNameTitle("MassScan/TauMass", "N lifetime vs coupling");
  BookHisto(fgTau);

  fgAcc = new TGraph();
  fgAcc->SetNameTitle("MassScan/AccMass", "Acceptance vs coupling");
  BookHisto(fgAcc);

  fgYield = new TGraph();
  fgYield->SetNameTitle("MassScan/YieldMass", "Yield per POT vs coupling");
  BookHisto(fgYield);

  fgExclusion = new TGraph();
  fgExclusion->SetNameTitle("TotalScan/Exclusion", "Sensitivity vs N mass and coupling");
  BookHisto(fgExclusion);
}

void HeavyNeutrinoAllInclusive::Process(Int_t) {

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
  Double_t DecayFactor    = 0.;
  Double_t Weight         = 0.;
  Double_t DProdProb      = 0.;
  Int_t counter           = 0;
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
	FillHisto("Selection/hXYCHODTrue", xCHOD/1000., yCHOD/1000.);
	FillHisto("Selection/hXYMUV3True", xMUV3/1000., yMUV3/1000.);
	if (p->GetParticleName() == "pi+" || p->GetParticleName() == "pi-") {
	  mom1 = p->GetInitial4Momentum();
	  p1 = TMath::Sqrt(mom1.Px()*mom1.Px()+mom1.Py()*mom1.Py()+mom1.Pz()*mom1.Pz());
	}
	else if (p->GetParticleName() == "mu+" || p->GetParticleName() == "mu-") {
	  mom2 = p->GetInitial4Momentum();
	  p2 = TMath::Sqrt(mom2.Px()*mom2.Px()+mom2.Py()*mom2.Py()+mom2.Pz()*mom2.Pz());
	}
	if (counter == 2) {
	  FillHisto("Selection/hMomPi",  p1/1000.);
	  FillHisto("Selection/hMomMu",  p2/1000.);
	  FillHisto("Selection/hP1vsP2", p1/1000., p2/1000.);
	  FillHisto("Selection/hInvMassMC", (mom1+mom2).M2()/1000.);
	}
      }
    }
  }

  // Scan on the coupling                                                                       

  for(Double_t fCoupling = fCouplingStart; fCoupling < fCouplingStop-fCouplingStep; fCoupling += fCouplingStep) {

    fUSquared = TMath::Power(10, fCoupling);
    fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
    fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
    fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;
    fCouplings[fCoupling] = fCoupling;

    if (GetWithMC()) {
      Event *evt = GetMCEvent();
      
      // Computation of coupling-related quantities of all HNLs (good and bad)
      
      for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
	KinePart *p = evt->GetKinePart(i);      
	if (p->GetParentID() == -1 && p->GetPDGcode() == 999) {
	  fNevents++;
	  point1.SetXYZ(p->GetProdPos().X(), p->GetProdPos().Y(), p->GetProdPos().Z());
	  point2.SetXYZ(0., 0., fLInitialFV);
	  momentum1.SetXYZ(p->GetInitial4Momentum().Px(), p->GetInitial4Momentum().Py(), p->GetInitial4Momentum().Pz());
	  MN = ComputeHNLMass(p);
          fMasses[round(MN)] = round(MN);
          if (fNevents[round(MN)].count(fCoupling) == 0)
            fNevents[round(MN)][fCoupling] = 0;
          fNevents[round(MN)][fCoupling]++;
	  gammaTot = GammaTot(MN);
	  HNLTau = tauN(MN);
          fGammaTot[round(MN)][fCoupling] = gammaTot;
          fTau[round(MN)][fCoupling] = HNLTau;
	  LReach = ComputeL(point1, point2, momentum1);
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
	  
	  Weight = DProdProb*fDDecayProb*NReachProb*NDecayProb*DecayFactor*LeptonUSquared;
          fSumAll[round(MN)][fCoupling] += Weight;
	  
	  // Some more plots of KinePart quantities
	  
	  FillHisto("Selection/hZDProd", p->GetPosAtCheckPoint(0).z());
	  FillHisto("Selection/hZDDecay", p->GetProdPos().Z());
	  FillHisto("Selection/hDTheta", p->GetPosAtCheckPoint(0).x());
	  FillHisto("Selection/hDLambda", p->GetPosAtCheckPoint(0).y());
	  FillHisto("Selection/hDPath", p->GetMomAtCheckPoint(0).X());
	  FillHisto("Selection/hDMom", p->GetMomAtCheckPoint(0).Y()/1000.);
	  FillHisto("Selection/hZHNLDecay", p->GetEndPos().Z()/1000.);
	  FillHisto("Selection/hHNLGamma", p->GetInitial4Momentum().Gamma());
	  FillHisto("Selection/hHNLReachProb", NReachProb);
	  FillHisto("Selection/hHNLDecayProb", NDecayProb);
	  FillHisto("Selection/hHNLTheta", p->GetMomAtCheckPoint(0).Z());
	  FillHisto("Selection/hHNLMom", p->GetMomAtCheckPoint(0).T()/1000.);
	  FillHisto("Selection/hWeight", Weight);

          FillHisto("CouplingScan/hReachCoupling",  fCoupling, NReachProb);
          FillHisto("CouplingScan/hDecayCoupling",  fCoupling, NDecayProb);
          FillHisto("CouplingScan/hWeightCoupling", fCoupling, Weight);

          FillHisto("MassScan/hReachMass",  fCoupling, NReachProb);
          FillHisto("MassScan/hDecayMass",  fCoupling, NDecayProb);
          FillHisto("MassScan/hWeightMass", fCoupling, Weight);
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

    if (TriggerOK) {
      
      // If real data
      // Select physics triggers and L0 trigger conditions
    
      if (!GetWithMC()) {  
	Bool_t L0OK = kFALSE;
	Bool_t L1OK = kFALSE;
      
	for (UInt_t i = 0; i < fID.size(); i++) {
	  L0OK |= TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0TPData, fID[i]);
	}
      
	if (L0OK) {
	
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
	
	  if (L1OK) {
	  
	    FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
	    CutID++;
	  
	    // Select two-track events
	  
	    std::vector<DownstreamTrack> Tracks = *(std::vector<DownstreamTrack>*) GetOutput("DownstreamTrackBuilder.Output");
	  
	    FillHisto("Selection/hNtracks", Tracks.size());
	  
	    if (Tracks.size() == 2) {
	    
	      FillHisto("Selection/hN2tracks", 0.5);
	      FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
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
		    TString name = Form("Selection/hXYSpec%dReco",j);
		    FillHisto(name, x/1000., y/1000.);
		  }
		}
		Double_t x  = Cand->xAtAfterMagnet(fzCHODPlane);
		Double_t y  = Cand->yAtAfterMagnet(fzCHODPlane);
		Double_t r  = sqrt(x*x+y*y);
		Double_t r1 = frMinCHOD;
		Double_t r2 = frMaxCHOD;
		if (r > r1 && r < r2)
		  FillHisto("Selection/hXYCHODReco", x/1000., y/1000.);
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
	    
	      FillHisto("Selection/hCDAvsZVertex_TrackToBeamlineInitial",      Zvertex1 / 1000.,         CDA1 / 1000.);
	      FillHisto("Selection/hCDAvsZVertex_TrackToBeamlineInitial",      Zvertex2 / 1000.,         CDA2 / 1000.);
	      FillHisto("Selection/hCDAvsZVertex_TrackToTrackInitial",          Zvertex / 1000.,          CDA / 1000.);
	      FillHisto("Selection/hCDAvsZVertex_TotMomToBeamlineInitial",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 1
	      FillHisto("Selection/hZvertexvsBeamlineDistInitial",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 1
	      FillHisto("Selection/hBeamlineDistvsTargetDist_TotMomInitial", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 1
	    
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
		  inAcc &= true;
	      }

	      if (inAcc) {
		
		FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
		CutID++;

		// Track selection, CUT 2: Chi2 and momentum cuts
	    
		if (ChiSquare1 < 20. && ChiSquare2 < 20.) {
	      
		  FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
		  CutID++;
	      
		  if (SpectrometerCand1->GetNChambers() >= 3 && SpectrometerCand2->GetNChambers() >= 3) {
		
		    FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
		    CutID++;
		
		    // Track selection, CUT 3: Opposite-charged tracks
		
		    if (Charge1 + Charge2 == 0) {

		      FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
		      CutID++;
		  
		      // Plot the number of candidates associated to each track, for MUV3 and CHOD
		  
		      FillHisto("Selection/hNMUV3CandAssocToTrack", Tracks[0].GetNMUV3AssociationRecords());
		      FillHisto("Selection/hNMUV3CandAssocToTrack", Tracks[1].GetNMUV3AssociationRecords());
		      FillHisto("Selection/hNCHODCandAssocToTrack", Tracks[0].GetNCHODAssociationRecords());
		      FillHisto("Selection/hNCHODCandAssocToTrack", Tracks[1].GetNCHODAssociationRecords());
		  
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

		      if (inAcc) {

			FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
			CutID++;
		    
			if (CHODAssoc) {

			  FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
			  CutID++;
		      
			  // Downstream track selection, CUT 5: Extrapolation and association to LKr
		      
			  Bool_t LKrAssoc = (Tracks[0].LKrAssociationExists() && Tracks[1].LKrAssociationExists());
		      
			  if (GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kLKr) && GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kLKr)) {
			
			    FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
			    CutID++;
			
			    if (LKrAssoc) {

			      FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
			      CutID++;
			  
			      // Downstream track selection, CUT 6: Extrapolation and association to MUV3
			  
			      if (GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kMUV3) && GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kMUV3)) {
			    
				FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
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
			    
				if (Assoc == 1 || Assoc == 2) {

				  FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
				  CutID++;

				  FillHisto("Selection/hCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 2
				  FillHisto("Selection/hZvertexvsBeamlineDistAfterDownstreamTrack",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 2
				  FillHisto("Selection/hBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 2
			      
				  // Compute time of MUV3 and CHOD candidates for better resolution wrt Spectrometer tracks
			      
				  Double_t MUV3Time;
			      
				  if (Assoc == 1)
				    MUV3Time = Tracks[0].GetMUV3Time(0);
				  else if (Assoc == 2)
				    MUV3Time = Tracks[1].GetMUV3Time(0);
			      
				  Double_t CHODTime1 = Tracks[0].GetCHODTime();
				  Double_t CHODTime2 = Tracks[1].GetCHODTime();
			      
				  // Plot time difference of the two tracks
			      
				  FillHisto("Selection/hDeltaTimeFromCHOD", CHODTime1 - CHODTime2);
			      
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
			      
				  FillHisto("Selection/hEoP", EoP1);
				  FillHisto("Selection/hEoP", EoP2);
				  FillHisto("Selection/hEoPMuVsPi", PiEoP, MuEoP);
			      
				  if (MuEoP < 0.2) {

				    FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
				    CutID++;
				
				    if (PiEoP > 0.2 && PiEoP < 0.8) {

				      FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
				      CutID++;
				  
				      FillHisto("Selection/hCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 3
				      FillHisto("Selection/hZvertexvsBeamlineDistAfterEnergyCuts",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 3
				      FillHisto("Selection/hBeamlineDistvsTargetDist_TotMomAfterEnergyCuts", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 3
				  
				      // Veto cuts, CUT 8: LAV veto
				  
				      fLAVMatching->SetReferenceTime((CHODTime1 + CHODTime2) / 2);
				  
				      if (GetWithMC())
					fLAVMatching->SetTimeCuts(99999, 99999);     // LAV is not time-aligned in MC
				      else
					fLAVMatching->SetTimeCuts(-10., 10.);
				  
				      if (!fLAVMatching->LAVHasTimeMatching(LAVEvent)) {

					FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
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
				    
					if (!fSAVMatching->SAVHasTimeMatching(IRCEvent, SACEvent)) {

					  FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
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
					  
					      FillHisto("Selection/hSingleAddEnLKrCand", LKrCand->GetTime() - MUV3Time, LKrCandE/1000.);
					  
					      if (LKrPos && LKrCandE > 40. && LKrTime) {
					      LKrEnergyInTime += LKrCandE;
					      LKrWeightedTime += LKrCandE * LKrCand->GetTime();
					      }
					      }
					  
					      if (LKrEnergyInTime > 1000.) {
					      LKrWeightedTime /= LKrEnergyInTime;
					      FillHisto("Selection/hAddEnLKrCand", LKrWeightedTime - MUV3Time, LKrEnergyInTime/1000.);
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
					  
					      FillHisto("Selection/hSingleAddEnLKrHit", LKrHitTime - MUV3Time, LKrHitE/1000.);
					  
					      if (LKrPos && LKrHitE > 40. && LKrTime) {
					      LKrEnergyInTime += LKrHitE;
					      LKrWeightedTime += LKrHitE * LKrHitTime;
					      }
					      }
					  
					      if (LKrEnergyInTime > 1000.) {
					      LKrWeightedTime /= LKrEnergyInTime;    
					      FillHisto("Selection/hAddEnLKrHit", LKrWeightedTime - MUV3Time, LKrEnergyInTime/1000.);
					      return;
					      }
					  */

					  FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
					  CutID++;
				      
					  FillHisto("Selection/hCDAvsZVertex_TotMomToBeamlineAfterVetoes",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 4
					  FillHisto("Selection/hZvertexvsBeamlineDistAfterVetoes",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 4
					  FillHisto("Selection/hBeamlineDistvsTargetDist_TotMomAfterVetoes", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 4
				      
					  // Geometrical cuts, CUT 11: Cut on CDA of two tracks
				      
					  if (CDA < 10.) {

					    FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
					    CutID++;
					
					    // Geometrical cuts, CUT 12: Cut on two-track vertex wrt beamline
					
					    if (BeamlineDist > 100.) {

					      FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
					      CutID++;
					  
					      // Geometrical cuts, CUT 13: Cut on Z of two-track vertex
					  
					      if (Zvertex > 102500. && Zvertex < 180000.) {

						FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
						CutID++;
					    
						// Geometrical cuts, CUT 14: Cut on CDA of each track wrt beam axis
					    
						if (CDA1 > 50. && CDA2 > 50.) {
  
						  FillHisto("Selection/hPhysicsEventsVsCuts", CutID);
						  CutID++;
					      
						  FillHisto("Selection/hCDAvsZVertex_TrackToBeamlineAfterCut",           Zvertex1 / 1000.,         CDA1 / 1000.);
						  FillHisto("Selection/hCDAvsZVertex_TrackToBeamlineAfterCut",           Zvertex2 / 1000.,         CDA2 / 1000.);
						  FillHisto("Selection/hCDAvsZVertex_TrackToTrackAfterCut",               Zvertex / 1000.,          CDA / 1000.);
						  FillHisto("Selection/hCDAvsZVertex_TotMomToBeamlineAfterGeomCuts",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 5
						  FillHisto("Selection/hZvertexvsBeamlineDistAfterGeomCuts",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 5
						  FillHisto("Selection/hBeamlineDistvsTargetDist_TotMomAfterGeomCuts", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 5
					      
						  // Plot CDA vs Zvertex for events surviving previous selections
					      
						  FillHisto("Selection/hCDAvsZVertex_TrackToBeamlineFinal",      Zvertex1 / 1000.,         CDA1 / 1000.);
						  FillHisto("Selection/hCDAvsZVertex_TrackToBeamlineFinal",      Zvertex2 / 1000.,         CDA2 / 1000.);
						  FillHisto("Selection/hCDAvsZVertex_TrackToTrackFinal",          Zvertex / 1000.,          CDA / 1000.);
						  FillHisto("Selection/hCDAvsZVertex_TotMomToBeamlineFinal",      Zvertex / 1000.,       CDAMom / 1000.);     // Reference plot, step 6
						  FillHisto("Selection/hZvertexvsBeamlineDistFinal",              Zvertex / 1000., BeamlineDist / 1000.);     // Reference plot, step 6
						  FillHisto("Selection/hBeamlineDistvsTargetDist_TotMomFinal", TargetDist / 1000., BeamlineDist / 1000.);     // Reference plot, step 6
					      
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
					      
						  FillHisto("Selection/hInvMassReco", invMass);
					      
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
							HNLTau = tauN(MN);
							LReach = ComputeL(point1, point2, momentum1);
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
						    
							Weight = DProdProb*fDDecayProb*NReachProb*NDecayProb*DecayFactor*LeptonUSquared;
							fSumGood[round(MN)][fCoupling] += Weight;
						      }
						    }
						  }
						}
					      }
					    }
					  }
					}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void HeavyNeutrinoAllInclusive::EndOfBurstUser() {
  
  FillHisto("hNbursts", 0.5);
}

void HeavyNeutrinoAllInclusive::EndOfJobUser() {
  
  // Retrieve selection histos

  fhNk3pi    = (TH1D*) fHisto.GetTH1("Selection/hNk3pi");
  fhNbursts  = (TH1D*) fHisto.GetTH1("Selection/hNbursts");
  fhNEvents  = (TH1D*) fHisto.GetTH1("Selection/hNEvents");
  fhN2tracks = (TH1D*) fHisto.GetTH1("Selection/hN2tracks");
  fhNtracks  = (TH1D*) fHisto.GetTH1("Selection/hNtracks");

  fhZDProd       = (TH1D*) fHisto.GetTH1("Selection/hZDProd");
  fhZDDecay      = (TH1D*) fHisto.GetTH1("Selection/hZDDecay");
  fhDTheta       = (TH1D*) fHisto.GetTH1("Selection/hDTheta");
  fhDLambda      = (TH1D*) fHisto.GetTH1("Selection/hDLambda");
  fhDPath        = (TH1D*) fHisto.GetTH1("Selection/hDPath");
  fhDMom         = (TH1D*) fHisto.GetTH1("Selection/hDMom");
  fhZHNLDecay    = (TH1D*) fHisto.GetTH1("Selection/hZHNLDecay");
  fhHNLGamma     = (TH1D*) fHisto.GetTH1("Selection/hHNLGamma");
  fhHNLDecayProb = (TH1D*) fHisto.GetTH1("Selection/hHNLDecayProb");
  fhHNLReachProb = (TH1D*) fHisto.GetTH1("Selection/hHNLReachProb");
  fhHNLTheta     = (TH1D*) fHisto.GetTH1("Selection/hHNLTheta");
  fhHNLMom       = (TH1D*) fHisto.GetTH1("Selection/hHNLMom");
  fhMomPi        = (TH1D*) fHisto.GetTH1("Selection/hMomPi");
  fhMomMu        = (TH1D*) fHisto.GetTH1("Selection/hMomMu");
  fhWeight       = (TH1D*) fHisto.GetTH1("Selection/hWeight");

  fhXYSpec0Reco = (TH2D*) fHisto.GetTH2("Selection/hXYSpec0Reco");
  fhXYSpec1Reco = (TH2D*) fHisto.GetTH2("Selection/hXYSpec1Reco");
  fhXYSpec2Reco = (TH2D*) fHisto.GetTH2("Selection/hXYSpec2Reco");
  fhXYSpec3Reco = (TH2D*) fHisto.GetTH2("Selection/hXYSpec3Reco");
  fhXYCHODReco  = (TH2D*) fHisto.GetTH2("Selection/hXYCHODReco");
  fhXYCHODTrue  = (TH2D*) fHisto.GetTH2("Selection/hXYCHODTrue");
  fhXYMUV3True  = (TH2D*) fHisto.GetTH2("Selection/hXYMUV3True");
  fhP1vsP2      = (TH2D*) fHisto.GetTH2("Selection/hP1vsP2");

  fhPhysicsEventsVsCuts = (TH1D*) fHisto.GetTH1("Selection/hPhysicsEventsVsCuts");
  
  fhCDAvsZVertex_TotMomToBeamlineInitial              = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TotMomToBeamlineInitial");
  fhCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TotMomToBeamlineAfterDownstreamTrack");
  fhCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts      = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TotMomToBeamlineAfterEnergyCuts");
  fhCDAvsZVertex_TotMomToBeamlineAfterVetoes          = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TotMomToBeamlineAfterVetoes");
  fhCDAvsZVertex_TotMomToBeamlineAfterGeomCuts        = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TotMomToBeamlineAfterGeomCuts");
  fhCDAvsZVertex_TotMomToBeamlineFinal                = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TotMomToBeamlineFinal");
  
  fhZvertexvsBeamlineDistInitial              = (TH2D*) fHisto.GetTH2("Selection/hZvertexvsBeamlineDistInitial");
  fhZvertexvsBeamlineDistAfterDownstreamTrack = (TH2D*) fHisto.GetTH2("Selection/hZvertexvsBeamlineDistAfterDownstreamTrack");
  fhZvertexvsBeamlineDistAfterEnergyCuts      = (TH2D*) fHisto.GetTH2("Selection/hZvertexvsBeamlineDistAfterEnergyCuts");
  fhZvertexvsBeamlineDistAfterVetoes          = (TH2D*) fHisto.GetTH2("Selection/hZvertexvsBeamlineDistAfterVetoes");
  fhZvertexvsBeamlineDistAfterGeomCuts        = (TH2D*) fHisto.GetTH2("Selection/hZvertexvsBeamlineDistAfterGeomCuts");
  fhZvertexvsBeamlineDistFinal                = (TH2D*) fHisto.GetTH2("Selection/hZvertexvsBeamlineDistFinal");
  
  fhCDAvsZVertex_TrackToBeamlineInitial  = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TrackToBeamlineInitial");
  fhCDAvsZVertex_TrackToTrackInitial     = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TrackToTrackInitial");
  fhCDAvsZVertex_TrackToBeamlineAfterCut = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TrackToBeamlineAfterCut");
  fhCDAvsZVertex_TrackToTrackAfterCut    = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TrackToTrackAfterCut");
  fhCDAvsZVertex_TrackToBeamlineFinal    = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TrackToBeamlineFinal");
  fhCDAvsZVertex_TrackToTrackFinal       = (TH2D*) fHisto.GetTH2("Selection/hCDAvsZVertex_TrackToTrackFinal");
  
  fhBeamlineDistvsTargetDist_TotMomInitial              = (TH2D*) fHisto.GetTH2("Selection/hBeamlineDistvsTargetDist_TotMomInitial");
  fhBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack = (TH2D*) fHisto.GetTH2("Selection/hBeamlineDistvsTargetDist_TotMomAfterDownstreamTrack");
  fhBeamlineDistvsTargetDist_TotMomAfterEnergyCuts      = (TH2D*) fHisto.GetTH2("Selection/hBeamlineDistvsTargetDist_TotMomAfterEnergyCuts");
  fhBeamlineDistvsTargetDist_TotMomAfterVetoes          = (TH2D*) fHisto.GetTH2("Selection/hBeamlineDistvsTargetDist_TotMomAfterVetoes");
  fhBeamlineDistvsTargetDist_TotMomAfterGeomCuts        = (TH2D*) fHisto.GetTH2("Selection/hBeamlineDistvsTargetDist_TotMomAfterGeomCuts");
  fhBeamlineDistvsTargetDist_TotMomFinal                = (TH2D*) fHisto.GetTH2("Selection/hBeamlineDistvsTargetDist_TotMomFinal");
  
  fhDeltaTimeFromCHOD     = (TH1D*) fHisto.GetTH1("Selection/hDeltaTimeFromCHOD");
  fhNMUV3CandAssocToTrack = (TH1D*) fHisto.GetTH1("Selection/hNMUV3CandAssocToTrack");
  fhNCHODCandAssocToTrack = (TH1D*) fHisto.GetTH1("Selection/hNCHODCandAssocToTrack");
  
  fhEoP       = (TH1D*) fHisto.GetTH1("Selection/hEoP");
  fhEoPMuVsPi = (TH2D*) fHisto.GetTH1("Selection/hEoPMuVsPi");
  
  fhSingleAddEnLKrHit  = (TH2D*) fHisto.GetTH2("Selection/hSingleAddEnLKrHit");
  fhSingleAddEnLKrCand = (TH2D*) fHisto.GetTH2("Selection/hSingleAddEnLKrCand");
  fhAddEnLKrHit        = (TH2D*) fHisto.GetTH2("Selection/hAddEnLKrHit");
  fhAddEnLKrCand       = (TH2D*) fHisto.GetTH2("Selection/hAddEnLKrCand");
  
  fhInvMassMC    = (TH1D*) fHisto.GetTH1("Selection/hInvMassMC");
  fhInvMassReco  = (TH1D*) fHisto.GetTH1("Selection/hInvMassReco");

  fhAcc   = (TH1D*) fHisto.GetTH1("Selection/hAcc");
  fhYield = (TH1D*) fHisto.GetTH1("Selection/hYield");

  // Retrieve scan histos

  fhReachCoupling    = (TH2D*)fHisto.GetTH2("hReachCoupling");
  fhDecayCoupling    = (TH2D*)fHisto.GetTH2("hDecayCoupling");
  fhWeightCoupling   = (TH2D*)fHisto.GetTH2("hWeightCoupling");

  fhReachMass    = (TH2D*)fHisto.GetTH2("hReachMass");
  fhDecayMass    = (TH2D*)fHisto.GetTH2("hDecayMass");
  fhWeightMass   = (TH2D*)fHisto.GetTH2("hWeightMass");

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

  fhReachCoupling   ->GetXaxis()->SetTitle("Log of coupling");
  fhDecayCoupling   ->GetXaxis()->SetTitle("Log of coupling");
  fhWeightCoupling  ->GetXaxis()->SetTitle("Log of coupling");

  fhReachMass   ->GetXaxis()->SetTitle("N mass [GeV]");
  fhDecayMass   ->GetXaxis()->SetTitle("N mass [GeV]");
  fhWeightMass  ->GetXaxis()->SetTitle("N mass [GeV]");

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

  fhReachCoupling   ->GetYaxis()->SetTitle("Reach probability");
  fhDecayCoupling   ->GetYaxis()->SetTitle("Decay probability");
  fhWeightCoupling  ->GetYaxis()->SetTitle("Weight");

  fhReachMass   ->GetYaxis()->SetTitle("Reach probability");
  fhDecayMass   ->GetYaxis()->SetTitle("Decay probability");
  fhWeightMass  ->GetYaxis()->SetTitle("Weight");

  // Acceptance computation                                                                       

  Double_t Coupling = 0.;
  Double_t MN = 0.;
  Int_t counter = 0;

  for (auto it = fCouplings.begin(); it != fCouplings.end(); it++) {
    Coupling = it->first;
    fgGammaTot->SetPoint(counter, Coupling, fGammaTot[Coupling]);
    fgTau     ->SetPoint(counter, Coupling, fTau     [Coupling]);
    fgAcc     ->SetPoint(counter, Coupling, fAcc     [Coupling]);
    fgYield   ->SetPoint(counter, Coupling, fYield   [Coupling]);
    counter++;
  }

  counter = 0;

  for (auto it = fMasses.begin(); it != fMasses.end(); it++) {
    MN = it->first;
    for (auto it1 = fCouplings.begin(); it1 != fCouplings.end(); it1++) {
      Coupling = it1->first;
      fAcc[MN][Coupling] = fSumGood[MN][Coupling]/fSumAll[MN][Coupling];
      fProb[MN][Coupling] = fSumAll[MN][Coupling]/fNevents[MN][Coupling];
      fYield[MN][Coupling] = fAcc[MN][Coupling]*fProb[MN][Coupling];
      counter++;
    }
  }

  // Exclusion plot                                                                                    

  counter= 0;

  for (auto it = fMasses.begin(); it != fMasses.end(); it++) {
    MN = it->first;
    for (auto it1 = fCouplings.begin(); it1 != fCouplings.end(); it1++) {
      Coupling = it1->first;
      if (fYield[MN][Coupling]*1.E18 > 2.3)
        fgExclusion->SetPoint(counter, MN/1000., Coupling);
      counter++;
    }
  }

  // Titles

  fgGammaTotCoupling->GetXaxis()->SetTitle("Log of coupling");
  fgTauCoupling     ->GetXaxis()->SetTitle("Log of coupling");
  fgAccCoupling     ->GetXaxis()->SetTitle("Log of coupling");
  fgYieldCoupling   ->GetXaxis()->SetTitle("Log of coupling");
  fgGammaTotCoupling->GetYaxis()->SetTitle("Total decay width [MeV]");
  fgTauCoupling     ->GetYaxis()->SetTitle("Lifetime [ns]");
  fgAccCoupling     ->GetYaxis()->SetTitle("Acceptance");
  fgYieldCoupling   ->GetYaxis()->SetTitle("Yield");
  fgGammaTotCoupling->SetLineColor(2);
  fgTauCoupling     ->SetLineColor(2);
  fgAccCoupling     ->SetLineColor(2);
  fgYieldCoupling   ->SetLineColor(2);
  fgGammaTotCoupling->SetLineWidth(3);
  fgTauCoupling     ->SetLineWidth(3);
  fgAccCoupling     ->SetLineWidth(3);
  fgYieldCoupling   ->SetLineWidth(3);
  fgGammaTotCoupling->Draw("AC");
  fgTauCoupling     ->Draw("AC");
  fgAccCoupling     ->Draw("AC");
  fgYieldCoupling   ->Draw("AC");
  fgGammaTotCoupling->Write();
  fgTauCoupling     ->Write();
  fgAccCoupling     ->Write();
  fgYieldCoupling   ->Write();

  fgGammaTotMass->GetXaxis()->SetTitle("N mass [GeV]");
  fgTauMass     ->GetXaxis()->SetTitle("N mass [GeV]");
  fgAccMass     ->GetXaxis()->SetTitle("N mass [GeV]");
  fgYieldMass   ->GetXaxis()->SetTitle("N mass [GeV]");
  fgGammaTotMass->GetYaxis()->SetTitle("Total decay width [MeV]");
  fgTauMass     ->GetYaxis()->SetTitle("Lifetime [ns]");
  fgAccMass     ->GetYaxis()->SetTitle("Acceptance");
  fgYieldMass   ->GetYaxis()->SetTitle("Yield");
  fgGammaTotMass->SetLineColor(2);
  fgTauMass     ->SetLineColor(2);
  fgAccMass     ->SetLineColor(2);
  fgYieldMass   ->SetLineColor(2);
  fgGammaTotMass->SetLineWidth(3);
  fgTauMass     ->SetLineWidth(3);
  fgAccMass     ->SetLineWidth(3);
  fgYieldMass   ->SetLineWidth(3);
  fgGammaTotMass->Draw("AC");
  fgTauMass     ->Draw("AC");
  fgAccMass     ->Draw("AC");
  fgYieldMass   ->Draw("AC");
  fgGammaTotMass->Write();
  fgTauMass     ->Write();
  fgAccMass     ->Write();
  fgYieldMass   ->Write();

  fgExclusion->GetXaxis()->SetTitle("N mass [GeV]");
  fgExclusion->GetYaxis()->SetTitle("Log of coupling");
  fgExclusion->SetLineColor(2);
  fgExclusion->SetLineWidth(3);
  fgExclusion->Draw("AC");
  fgExclusion->Write();

  // Plot residual number of events after each cut

  const int NCuts = 25;
  const char *CutNames[NCuts]  = {"Total", "TriggerOK", "2 tracks", "Straw0 acc", "Straw1 acc", "Straw2 acc", "Straw3 acc", "Chi2", "Straw chambers", "Charge", "CHOD acc", "CHOD assoc", "LKr acc", "LKr assoc", "MUV3 acc", "MUV3 assoc", "Mu E/p", "Pi E/p", "LAV veto", "SAV veto", "LKr veto", "CDA tracks", "Beam distance", "Z vertex", "CDA beam"};

  for (Int_t i = 1; i <= NCuts; i++)
    fhPhysicsEventsVsCuts->GetXaxis()->SetBinLabel(i, CutNames[i-1]);
  
  fhPhysicsEventsVsCuts->GetXaxis()->LabelsOption("v");

  SaveAllPlots();

  return;
}

HeavyNeutrinoAllInclusive::~HeavyNeutrinoAllInclusive() {

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

  fhReachCoupling    = nullptr;
  fhDecayCoupling    = nullptr;
  fhWeightCoupling   = nullptr;
  fgAccCoupling      = nullptr;
  fgYieldCoupling    = nullptr;
  fgGammaTotCoupling = nullptr;
  fgTauCoupling      = nullptr;

  fhReachMass    = nullptr;
  fhDecayMass    = nullptr;
  fhWeightMass   = nullptr;
  fgAccMass      = nullptr;
  fgYieldMass    = nullptr;
  fgGammaTotMass = nullptr;
  fgTauMass      = nullptr;

  fgExclusion = nullptr;
}
