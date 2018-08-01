// --------------------------------------------------------------        
//                                                                                  
// History:                                                                         
//                                                                             
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018    
//                                                                                
// ---------------------------------------------------------------              
/// \class HeavyNeutrinoScan
/// \Brief                                                                       
/// Produce expected yield per POT for HNL MC samples, for single values of HNL mass and coupling, for a scan on the coupling or on the mass, or for both
/// \EndBrief                                                                  
//
/// \Detailed                                                                        
/// After retrieving all HNL weights using the tool called HNLWeight, and checking if an event passed the selection implemented in the analyzer HeavyNeutrino, acceptance and yield per POT are computed. 
/// Different sets of plots are produced, all together at once:
/// for one value of HNL mass and coupling (that can be set as external parameters; default values are MassForSingleValue = 1000 MeV and CouplingForSingleValue = -6 = Log(10^-6)), quantities related to the D meson, the HNL and its product are produced;acceptance and yield are also computed; the yield is computed as a function of the HNL momentum as well;
/// a scan is enabled, to produce plots as a function of the coupling;
/// if the MC sample contains several masses, plots are produced as a function of the HNL mass; in this case, an expected exclusion plot is produced as a function of the HNL mass and coupling. Thus, the analyzer is able to run either on MC samples of just one HNL mass or on samples containing HNLs of different masses;
/// for HNL masses of 1 GeV, several plots are produced for toy-MC comparison: momentum and transverse momentum of the D meson, HNL and muon daughter (both for all HNLs and only for the good HNL), for two decay chains, as produced in the MC routine: DS->Nmu; N->pimu, and D0->KNmu; N->pimu.
///                  
/// This analyzer makes use of two ToolsLib, called HNLFunctions and HNLWeight.  
/// The values of the ratios between specific-flavour couplings can be either set as external parameters from command line or taken as default values.               
/// For example, if the user sets USquared = 1.E-10, UeSquaredRatio = 5., UmuSquaredRatio = 1., UtauSquaredRatio = 3.5, the specific-flavour coupling values will be: UeSquared = 5.25E-11, UmuSquared = 1.05E-11, UtauSquared = 3.68E-11.   
/// The values of the initial and final coupling and the scan step can be either set as external parameters from command line or taken as default values.               
/// For example, if the user sets CouplingStart = -3, CouplingStop = -2. and CouplingStep = 0.5, the scan will be performed for Log(U2) = [-3., -2.5, -2.], this meaning U2 = [1.E-3., 1.E-2.5, 1.E-2.].       
/// The user must also set the initial and final mass and the mass step, as generated in the MC simulation, in GeV.
/// The values of the initial and final momentum and its step can be either set as external parameters from command line or taken as default values. These are used to produce plots of the acceptance and yield per POT as a function of the HNL momentum (default values are MomStart = 0. GeV, MomStop = 150. GeV and MomStep = 5. GeV).
/// Other parameters to be set from command line are:
/// length of the beginning Z of the fiducial volume (default InitialFV = 102500. mm); length of the FV (default LFV = 77500. mm); decay mode of the HNLs, as set in the MC macro (default Mode = 0, which is pimu mode).                               
///
/// \author Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch)                           
/// \EndDetailed                                           

#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include <TGraphAsymmErrors.h>
#include "MCSimple.hh"
#include "MCInfo.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "BeamParameters.hh"
#include "GeometricAcceptance.hh"
#include "HNLFunctions.hh"
#include "HNLWeight.hh"
#include "HeavyNeutrinoScan.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

/// \class HeavyNeutrinoScan

HeavyNeutrinoScan::HeavyNeutrinoScan(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrinoScan") {

  fReadingData = GetIsTree();

  RequestAllMCTrees();

  AddParam("USquared", &fUSquared, 1.E-6);
  AddParam("InitialUeSquaredRatio", &fInitialUeSquaredRatio, 1.);
  AddParam("InitialUmuSquaredRatio", &fInitialUmuSquaredRatio, 16.);
  AddParam("InitialUtauSquaredRatio", &fInitialUtauSquaredRatio, 3.8);
  AddParam("CouplingStart", &fCouplingStart, -8.);
  AddParam("CouplingStop", &fCouplingStop, -2.);
  AddParam("CouplingStep", &fCouplingStep, 1.);
  AddParam("MassStart", &fMassStart, 0.3);
  AddParam("MassStop", &fMassStop, 1.8);
  AddParam("MassStep", &fMassStep, 0.1);
  AddParam("InitialFV", &fInitialFV, 102500.);
  AddParam("LFV", &fLFV, 77500.);
  AddParam("Mode", &fMode, 0);
  AddParam("MomStop", &fMomStop, 150.);
  AddParam("MomStart", &fMomStart, 0.);
  AddParam("MomStep", &fMomStep, 5.);
  AddParam("MassForSingleValue", &fMassForSingleValue, 1.);
  AddParam("CouplingForSingleValue", &fCouplingForSingleValue, -6.);

  fUeSquaredRatio = fInitialUeSquaredRatio;
  fUmuSquaredRatio = fInitialUmuSquaredRatio;
  fUtauSquaredRatio = fInitialUtauSquaredRatio;
  fNMom = round((std::abs(fMomStop - fMomStart))/fMomStep);
  fN = round((std::abs(fCouplingStop-fCouplingStart))/fCouplingStep);
  fNMass = round((std::abs(fMassStop-fMassStart))/fMassStep);
  fUeSquared = fUSquared/(fInitialUeSquaredRatio + fInitialUmuSquaredRatio + fInitialUtauSquaredRatio)*fInitialUeSquaredRatio;
  fUmuSquared = fUSquared/(fInitialUeSquaredRatio + fInitialUmuSquaredRatio + fInitialUtauSquaredRatio)*fInitialUmuSquaredRatio;
  fUtauSquared = fUSquared/(fInitialUeSquaredRatio + fInitialUmuSquaredRatio + fInitialUtauSquaredRatio)*fInitialUtauSquaredRatio;

  fCDAcomp = new TwoLinesCDA();
  fDistcomp = new PointLineDistance();

  if (fReadingData) {  
    fGammaTotFile.open("GammaTotFile.txt", ios::out);
    fTauFile.open("TauFile.txt", ios::out);
    fNEventsFile.open("NEventsFile.txt", ios::out);
    fSumGoodFile.open("SumGoodFile.txt", ios::out);
  }
  else {
    fGammaTotFile.open("GammaTotFile.txt", ios::in);
    fTauFile.open("TauFile.txt", ios::in);
    fNEventsFile.open("NEventsFile.txt", ios::in);
    fSumGoodFile.open("SumGoodFile.txt", ios::in);
  }
}

void HeavyNeutrinoScan::InitHist() {

  if (fReadingData) {
    
    // One value of mass and coupling

    BookHisto("SingleValue/hZMotherProd",   new TH1D("ZMotherProd", "N mother production point", 100000, -1., 27.));
    BookHisto("SingleValue/hZDProd",        new TH1D("ZDProd", "D meson production point", 100000, -1., 27.));
    BookHisto("SingleValue/hZTauProd",      new TH1D("ZTauProd", "#tau production point", 100000, -1., 27.));
    BookHisto("SingleValue/hZDDecay",       new TH1D("ZDDecay", "N production point", 100000, -1., 27.));
    BookHisto("SingleValue/hDTheta",        new TH1D("DTheta", "N mother polar angle", 100,  0., 0.3));
    BookHisto("SingleValue/hDLambda",       new TH1D("DLambda", "N mother decay length", 100, -1., 40.));
    BookHisto("SingleValue/hDPath",         new TH1D("DPath", "N mother path along Z", 100, -1., 50.));
    BookHisto("SingleValue/hDMom",          new TH1D("DMom", "N mother momentum", 100, -1., 170.));
    BookHisto("SingleValue/hZHNLDecay",     new TH1D("ZHNLDecay", "N decay point", 100., 90., 190.));
    BookHisto("SingleValue/hHNLGamma",      new TH1D("HNLGamma", "N Lorentz gamma", 50., 0., 170.));
    BookHisto("SingleValue/hHNLDecayProb",  new TH1D("HNLDecayProb", "N decay probability", 100., 0., 0.006));
    BookHisto("SingleValue/hHNLReachProb",  new TH1D("HNLReachProb", "N reach probability", 100., 0.99, 1.));
    BookHisto("SingleValue/hHNLTheta",      new TH1D("HNLTheta", "N polar angle", 100., 0., 0.5));
    BookHisto("SingleValue/hHNLMom",        new TH1D("HNLMom", "N momentum", 100., -0.5, 200.));
    BookHisto("SingleValue/hWeight",        new TH1D("Weight", "N weight", 1000, 1.E-20, 1.E-7));
    BookHisto("SingleValue/hCoupling",      new TH1D("Coupling", "Coupling", 100, -10., 0.));
    BookHisto("SingleValue/hMass",          new TH1D("Mass", "Mass", 10, 0.3, 1.2));    
    BookHisto("SingleValue/hTransverse",    new TH1D("Transverse", "N transverse position at decay point", 100, 0., 2.));
    BookHisto("SingleValue/hDThetaMom",     new TH2D("DThetaMom", "D meson polar angle vs momentum", 100, 0., 150., 50, 0., 0.5));
    BookHisto("SingleValue/hXYDecay",       new TH2D("XYDecay", "X, Y of N at decay point", 100, -1., 1., 100, -1., 1.));
    BookHisto("SingleValue/hXYMuon",        new TH2D("XYMuon", "X, Y of muon daughter at CH1", 100, -1., 1., 100, -1., 1.));
    BookHisto("SingleValue/hXYPion",        new TH2D("XYPion", "X, Y of pion daughter at CH1", 100, -1., 1., 100, -1., 1.));
    BookHisto("SingleValue/hBeamvsTarTrue", new TH2D("BeamvsTarTrue", "N true trajectory (target)", 20, 0., 0.2, 100, 0., 1.)); // target extrapolation for all events
    BookHisto("SingleValue/hBeamvsTAXTrue", new TH2D("BeamvsTAXTrue", "N true trajectory (TAX)", 20, 0., 0.2, 100, 0., 1.)); // TAX extrapolation for all events

    BookHisto("SingleValue/hG", new TH1D("sG", "", fNMom, fMomStart, fMomStop));
    BookHisto("SingleValue/hR", new TH1D("sR", "", fNMom, fMomStart, fMomStop));
    BookHisto("SingleValue/hA", new TH1D("sA", "", fNMom, fMomStart, fMomStop));
    BookHisto("SingleValue/hN", new TH1D("sN", "", fNMom, fMomStart, fMomStop));

    // Coupling scan 

    BookHisto("CouplingScan/hReachCoupling",  new TH2D("ReachCoupling", "Probability of N reaching the FV vs coupling", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2., 1000, -0.1, 1.1));
    BookHisto("CouplingScan/hDecayCoupling",  new TH2D("DecayCoupling", "Probability of N decaying in the FV vs coupling", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2., 1000, -0.1, 1.1));
    BookHisto("CouplingScan/hProbCoupling",   new TH2D("ProbCoupling", "Probability of N reaching and decaying in the FV vs coupling", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2., 1000, -0.1, 1.1));
    BookHisto("CouplingScan/hWeightCoupling", new TH2D("WeightCoupling", "N weight vs coupling", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2., 100000, 1.E-12, 1.E-7));
    BookHisto("CouplingScan/hEnergyCoupling", new TH2D("EnergyCoupling", "N energy vs coupling", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2., 100, 0., 100.));

    BookHisto("CouplingScan/hG",  new TH1D("cG", "",  fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("CouplingScan/hA",  new TH1D("cA", "",  fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("CouplingScan/hR",  new TH1D("cR", "",  fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("CouplingScan/hG1", new TH1D("cG1", "", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("CouplingScan/hA1", new TH1D("cA1", "", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("CouplingScan/hR1", new TH1D("cR1", "", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("CouplingScan/hG2", new TH1D("cG2", "", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("CouplingScan/hA2", new TH1D("cA2", "", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("CouplingScan/hR2", new TH1D("cR2", "", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("CouplingScan/hN",  new TH1D("cN",  "", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("CouplingScan/hN1", new TH1D("cN1", "", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("CouplingScan/hN2", new TH1D("cN2", "", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
   
    // Mass scan

    BookHisto("MassScan/hReachMass",  new TH2D("ReachMass", "Probability of N reaching the FV vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 1000, -0.1, 1.1));
    BookHisto("MassScan/hDecayMass",  new TH2D("DecayMass", "Probability of N decaying in the FV vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 1000, -0.1, 1.1));
    BookHisto("MassScan/hProbMass",   new TH2D("ProbMass", "Probability of N reaching and decaying in the FV vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 1000, -0.1, 1.1));
    BookHisto("MassScan/hWeightMass", new TH2D("WeightMass", "N weight vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 10000, 1.E-13, 1.E-9));
    BookHisto("MassScan/hEnergyMass", new TH2D("EnergyMass", "N energy vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 100, 0., 100.));

    BookHisto("MassScan/hG",  new TH1D("mG", "",  fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    BookHisto("MassScan/hR",  new TH1D("mR", "",  fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    BookHisto("MassScan/hA",  new TH1D("mA", "",  fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    BookHisto("MassScan/hG1", new TH1D("mG1", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    BookHisto("MassScan/hA1", new TH1D("mA1", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    BookHisto("MassScan/hR1", new TH1D("mR1", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    BookHisto("MassScan/hG2", new TH1D("mG2", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    BookHisto("MassScan/hA2", new TH1D("mA2", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    BookHisto("MassScan/hR2", new TH1D("mR2", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    BookHisto("MassScan/hN",  new TH1D("mN",  "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    BookHisto("MassScan/hN1", new TH1D("mN1", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    BookHisto("MassScan/hN2", new TH1D("mN2", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
    
    // Toy-MC comparison: DS

    BookHisto("ToyMC/DS/hpDS",   new TH1D("pDS", "D_{S} momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 40, 0., 200.));
    BookHisto("ToyMC/DS/hpDS1",  new TH1D("pDS1", "D_{S} momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 20, 40., 200.));
    BookHisto("ToyMC/DS/hptDS",  new TH1D("ptDS", "D_{S} transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 10, 0., 5.));
    BookHisto("ToyMC/DS/hptDS1", new TH1D("ptDS1", "D_{S} transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 6, 0., 3.));

    BookHisto("ToyMC/DS/hpNS",   new TH1D("pNS", "N momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 50, 0., 300.));
    BookHisto("ToyMC/DS/hpNS1",  new TH1D("pNS1", "N momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 50, 0., 300.));
    BookHisto("ToyMC/DS/hptNS",  new TH1D("ptNS", "N transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 20, 0., 2.));
    BookHisto("ToyMC/DS/hptNS1", new TH1D("ptNS1", "N transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 20, 0., 2.));

    BookHisto("ToyMC/DS/hpmuS",   new TH1D("pmuS", "Muon daughter momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 50, 0., 300.));
    BookHisto("ToyMC/DS/hpmuS1",  new TH1D("pmuS1", "Muon daughter momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 50, 0., 300.));
    BookHisto("ToyMC/DS/hptmuS",  new TH1D("ptmuS", "Muon daughter transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 20, 0., 2.));
    BookHisto("ToyMC/DS/hptmuS1", new TH1D("ptmuS1", "Muon daughter transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 20, 0., 2.));

    // Toy-MC comparison: D0

    BookHisto("ToyMC/D0/hpD0",   new TH1D("pD0", "D^{0} momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 40, 0., 200.));
    BookHisto("ToyMC/D0/hpD01",  new TH1D("pD01", "D^{0} momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 40, 60., 300.));
    BookHisto("ToyMC/D0/hptD0",  new TH1D("ptD0", "D^{0} transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 10, 0., 5.));
    BookHisto("ToyMC/D0/hptD01", new TH1D("ptD01", "D^{0} transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 6, 0., 3.));

    BookHisto("ToyMC/D0/hpN0",   new TH1D("pN0", "N momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 50, 0., 300.));
    BookHisto("ToyMC/D0/hpN01",  new TH1D("pN01", "N momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 50, 0., 300.));
    BookHisto("ToyMC/D0/hptN0",  new TH1D("ptN0", "N transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 20, 0., 2.));
    BookHisto("ToyMC/D0/hptN01", new TH1D("ptN01", "N transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 20, 0., 2.));

    BookHisto("ToyMC/D0/hpmu0",   new TH1D("pmu0", "Muon daughter momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 50, 0., 300.));
    BookHisto("ToyMC/D0/hpmu01",  new TH1D("pmu01", "Muon daughter momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 50, 0., 300.));
    BookHisto("ToyMC/D0/hptmu0",  new TH1D("ptmu0", "Muon daughter transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 20, 0., 2.));
    BookHisto("ToyMC/D0/hptmu01", new TH1D("ptmu01", "Muon daughter transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 20, 0., 2.));
    
    // Needed for error graphs
    
    fHisto.GetTH1("SingleValue/hG")->Sumw2();
    fHisto.GetTH1("SingleValue/hR")->Sumw2();
    fHisto.GetTH1("SingleValue/hA")->Sumw2();
    fHisto.GetTH1("SingleValue/hN")->Sumw2();
    fHisto.GetTH1("SingleValue/hTransverse")->Sumw2();
    fHisto.GetTH1("CouplingScan/hG")->Sumw2();
    fHisto.GetTH1("CouplingScan/hG1")->Sumw2();
    fHisto.GetTH1("CouplingScan/hG2")->Sumw2();
    fHisto.GetTH1("CouplingScan/hR")->Sumw2();
    fHisto.GetTH1("CouplingScan/hR1")->Sumw2();
    fHisto.GetTH1("CouplingScan/hR2")->Sumw2();
    fHisto.GetTH1("CouplingScan/hA")->Sumw2();
    fHisto.GetTH1("CouplingScan/hA1")->Sumw2();
    fHisto.GetTH1("CouplingScan/hA2")->Sumw2();
    fHisto.GetTH1("CouplingScan/hN")->Sumw2();
    fHisto.GetTH1("CouplingScan/hN1")->Sumw2();
    fHisto.GetTH1("CouplingScan/hN2")->Sumw2();
    fHisto.GetTH1("MassScan/hG")->Sumw2();
    fHisto.GetTH1("MassScan/hG1")->Sumw2();
    fHisto.GetTH1("MassScan/hG2")->Sumw2();
    fHisto.GetTH1("MassScan/hR")->Sumw2();
    fHisto.GetTH1("MassScan/hR1")->Sumw2();
    fHisto.GetTH1("MassScan/hR2")->Sumw2();
    fHisto.GetTH1("MassScan/hA")->Sumw2();
    fHisto.GetTH1("MassScan/hA1")->Sumw2();
    fHisto.GetTH1("MassScan/hA2")->Sumw2();
    fHisto.GetTH1("MassScan/hN")->Sumw2();
    fHisto.GetTH1("MassScan/hN1")->Sumw2();
    fHisto.GetTH1("MassScan/hN2")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hpDS")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hpDS1")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hptDS")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hptDS1")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hpNS")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hpNS1")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hptNS")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hptNS1")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hpmuS")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hpmuS1")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hptmuS")->Sumw2();
    fHisto.GetTH1("ToyMC/DS/hptmuS1")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hpD0")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hpD01")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hptD0")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hptD01")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hpN0")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hpN01")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hptN0")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hptN01")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hpmu0")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hpmu01")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hptmu0")->Sumw2();
    fHisto.GetTH1("ToyMC/D0/hptmu01")->Sumw2();
  }
  else {

    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/ZMotherProd", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/ZDProd", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/ZTauProd", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/ZDDecay", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/DTheta", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/DLambda", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/DPath", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/DMom", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/ZHNLDecay", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/HNLGamma", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/HNLDecayProb", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/HNLReachProb", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/HNLTheta", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/HNLMom", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/Weight", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/Coupling", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/Mass", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/Transverse", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "SingleValue/DThetaMom", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "SingleValue/XYDecay", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "SingleValue/XYMuon", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "SingleValue/XYPion", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "SingleValue/BeamvsTarTrue", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "SingleValue/BeamvsTAXTrue", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "CouplingScan/ReachCoupling", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "CouplingScan/DecayCoupling", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "CouplingScan/ProbCoupling", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "CouplingScan/WeightCoupling", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "CouplingScan/EnergyCoupling", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "MassScan/ReachMass", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "MassScan/DecayMass", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "MassScan/ProbMass", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "MassScan/WeightMass", true))->Write();
    ((TH2D*)RequestHistogram(fAnalyzerName, "MassScan/EnergyMass", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/pDS", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/pDS1", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/ptDS", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/ptDS1", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/pNS", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/pNS1", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/ptNS", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/ptNS1", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/pmuS", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/pmuS1", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/ptmuS", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/DS/ptmuS1", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/pD0", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/pD01", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/ptD0", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/ptD01", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/pN0", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/pN01", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/ptN0", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/ptN01", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/pmu0", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/pmu01", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/ptmu0", true))->Write();
    ((TH1D*)RequestHistogram(fAnalyzerName, "ToyMC/D0/ptmu01", true))->Write();

    // One value

    BookHisto("SingleValue/ErrorAccMomSel", new TGraphAsymmErrors());
    fHisto.GetTGraph("SingleValue/ErrorAccMomSel")->SetNameTitle("SingleValue/ErrorAccMomSel", "Acceptance vs N momentum, selection");
    BookHisto("SingleValue/ErrorAccMomReg", new TGraphAsymmErrors());
    fHisto.GetTGraph("SingleValue/ErrorAccMomReg")->SetNameTitle("SingleValue/ErrorAccMomReg", "Acceptance vs N momentum, regeneration");
    BookHisto("SingleValue/ErrorAccMomFV", new TGraphAsymmErrors());
    fHisto.GetTGraph("SingleValue/ErrorAccMomFV")->SetNameTitle("SingleValue/ErrorAccMomFV", "Acceptance vs N momentum, FV");

    BookHisto("SingleValue/ErrorYieldMom", new TGraphAsymmErrors());
    fHisto.GetTGraph("SingleValue/ErrorYieldMom")->SetNameTitle("SingleValue/ErrorYieldMom", "Yield per POT vs N momentum");

    // Coupling

    BookHisto("CouplingScan/MeanCoupling", new TGraph());
    fHisto.GetTGraph("CouplingScan/MeanCoupling")->SetNameTitle("CouplingScan/MeanCoupling", "Mean probability of N reaching and decaying in the FV vs coupling");
    BookHisto("CouplingScan/GammaTotCoupling", new TGraph());
    fHisto.GetTGraph("CouplingScan/GammaTotCoupling")->SetNameTitle("CouplingScan/GammaTotCoupling", "N total decay width vs coupling");
    BookHisto("CouplingScan/TauCoupling", new TGraph());
    fHisto.GetTGraph("CouplingScan/TauCoupling")->SetNameTitle("CouplingScan/TauCoupling", "N mean lifetime vs coupling");

    BookHisto("CouplingScan/ErrorAccCouplingSel", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSel")->SetNameTitle("CouplingScan/ErrorAccCouplingSel", "Selection");
    BookHisto("CouplingScan/ErrorAccCouplingReg", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingReg")->SetNameTitle("CouplingScan/ErrorAccCouplingReg", "Regeneration");
    BookHisto("CouplingScan/ErrorAccCouplingFV", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFV")->SetNameTitle("CouplingScan/ErrorAccCouplingFV", "FV");
    BookHisto("CouplingScan/ErrorYieldCoupling", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorYieldCoupling")->SetNameTitle("CouplingScan/ErrorYieldCoupling", "All");

    BookHisto("CouplingScan/ErrorAccCouplingSelTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSelTarget")->SetNameTitle("CouplingScan/ErrorAccCouplingSelTarget", "Selection Target");
    BookHisto("CouplingScan/ErrorAccCouplingRegTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingRegTarget")->SetNameTitle("CouplingScan/ErrorAccCouplingRegTarget", "Regeneration Target");
    BookHisto("CouplingScan/ErrorAccCouplingFVTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFVTarget")->SetNameTitle("CouplingScan/ErrorAccCouplingFVTarget", "FV Target");
    BookHisto("CouplingScan/ErrorYieldCouplingTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTarget")->SetNameTitle("CouplingScan/ErrorYieldCouplingTarget", "Target");

    BookHisto("CouplingScan/ErrorAccCouplingSelTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSelTAX")->SetNameTitle("CouplingScan/ErrorAccCouplingSelTAX", "Selection TAX");
    BookHisto("CouplingScan/ErrorAccCouplingRegTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingRegTAX")->SetNameTitle("CouplingScan/ErrorAccCouplingRegTAX", "Regeneration TAX");
    BookHisto("CouplingScan/ErrorAccCouplingFVTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFVTAX")->SetNameTitle("CouplingScan/ErrorAccCouplingFVTAX", "FV TAX");
    BookHisto("CouplingScan/ErrorYieldCouplingTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTAX")->SetNameTitle("CouplingScan/ErrorYieldCouplingTAX", "TAX");

    // Mass

    BookHisto("MassScan/MeanMass", new TGraph());
    fHisto.GetTGraph("MassScan/MeanMass")->SetNameTitle("MassScan/MeanMass", "Mean probability of N reaching and decaying in the FV vs N mass");
    BookHisto("MassScan/GammaTotMass", new TGraph());
    fHisto.GetTGraph("MassScan/GammaTotMass")->SetNameTitle("MassScan/GammaTotMass", "N total decay width vs N mass");
    BookHisto("MassScan/TauMass", new TGraph());
    fHisto.GetTGraph("MassScan/TauMass")->SetNameTitle("MassScan/TauMass", "N mean lifetime vs N mass");

    BookHisto("MassScan/ErrorAccMassSel", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassSel")->SetNameTitle("MassScan/ErrorAccMassSel", "Selection");
    BookHisto("MassScan/ErrorAccMassReg", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassReg")->SetNameTitle("MassScan/ErrorAccMassReg", "Regeneration");
    BookHisto("MassScan/ErrorAccMassFV", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassFV")->SetNameTitle("MassScan/ErrorAccMassFV", "FV");
    BookHisto("MassScan/ErrorYieldMass", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorYieldMass")->SetNameTitle("MassScan/ErrorYieldMass", "All");

    BookHisto("MassScan/ErrorAccMassSelTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassSelTarget")->SetNameTitle("MassScan/ErrorAccMassSelTarget", "Selection Target");
    BookHisto("MassScan/ErrorAccMassRegTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassRegTarget")->SetNameTitle("MassScan/ErrorAccMassRegTarget", "Regeneration Target");
    BookHisto("MassScan/ErrorAccMassFVTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassFVTarget")->SetNameTitle("MassScan/ErrorAccMassFVTarget", "FV Target");
    BookHisto("MassScan/ErrorYieldMassTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorYieldMassTarget")->SetNameTitle("MassScan/ErrorYieldMassTarget", "Target");

    BookHisto("MassScan/ErrorAccMassSelTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassSelTAX")->SetNameTitle("MassScan/ErrorAccMassSelTAX", "Selection TAX");
    BookHisto("MassScan/ErrorAccMassRegTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassRegTAX")->SetNameTitle("MassScan/ErrorAccMassRegTAX", "Regeneration TAX");
    BookHisto("MassScan/ErrorAccMassFVTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassFVTAX")->SetNameTitle("MassScan/ErrorAccMassFVTAX", "FV TAX");
    BookHisto("MassScan/ErrorYieldMassTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorYieldMassTAX")->SetNameTitle("MassScan/ErrorYieldMassTAX", "TAX");

    // Total

    BookHisto("TotalScan/hExclusion", new TH2D("Exclusion", "Sensitivity as a function of N mass and coupling", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));
    BookHisto("TotalScan/Contours", new TGraph());
    fHisto.GetTGraph("TotalScan/Contours")->SetNameTitle("TotalScan/Contours", "Contours for sensitivity");
  }
}

void HeavyNeutrinoScan::Process(Int_t) {

  if (fReadingData) {

    Double_t fCoupling = -999.;  
    Double_t MN = 0.;
    Double_t HNLTau = 0.;
    Double_t gammaTot = 0.;
    Double_t momN = 0.;
    Double_t momBin = 0.;
    Double_t NDecayProb = 0.;
    Double_t NReachProb = 0.;
    Double_t Weight = 0.;
    Double_t scale = 1.E10;
    Bool_t isGood = false;

    Bool_t IsHNLGood = *(Bool_t*)GetOutput("HeavyNeutrino.Output");

    // Scan on the coupling                          
  
    for(Int_t couplingIndex  = fCouplingStart*10; couplingIndex <= fCouplingStop*10; couplingIndex += fCouplingStep*10) {
      fCoupling = couplingIndex/10.;
      fUSquared = TMath::Power(10., fCoupling);
      fUeSquaredRatio = fInitialUeSquaredRatio;
      fUmuSquaredRatio = fInitialUmuSquaredRatio;
      fUtauSquaredRatio = fInitialUtauSquaredRatio;
      fUeSquared = fUSquared/(fInitialUeSquaredRatio + fInitialUmuSquaredRatio + fInitialUtauSquaredRatio)*fInitialUeSquaredRatio;
      fUmuSquared = fUSquared/(fInitialUeSquaredRatio + fInitialUmuSquaredRatio + fInitialUtauSquaredRatio)*fInitialUmuSquaredRatio;
      fUtauSquared = fUSquared/(fInitialUeSquaredRatio + fInitialUmuSquaredRatio + fInitialUtauSquaredRatio)*fInitialUtauSquaredRatio;
      fCouplings[fCoupling] = fCoupling;

      // Retrieve HNL weights

      if (GetWithMC()) {
	Event *evt = GetMCEvent();

	std::vector<std::map<std::string, Double_t>> Weights = ComputeWeight(evt, fUSquared, fInitialUeSquaredRatio, fInitialUmuSquaredRatio, fInitialUtauSquaredRatio, fLInitialFV, fLFV, fMode);

	for (UInt_t i = 0; i < Weights.size(); i++) {
	  MN = Weights[i]["Mass"];
	  HNLTau = Weights[i]["Lifetime"];
	  momN = Weights[i]["Momentum"];
	  gammaTot = GammaTot(MN);
	  NReachProb = Weights[i]["ReachProb"];
	  NDecayProb = Weights[i]["DecayProb"];
	  Weight = Weights[i]["Weight"];    
	  isGood = Weights[i]["IsGood"];
	  fGammaTot[round(MN)][fCoupling] = gammaTot;
	  fTau[round(MN)][fCoupling] = HNLTau;
	  fMasses[round(MN)] = round(MN);

	  // Maps for exclusion plot

	  if (fNEvents[round(MN)].count(fCoupling) == 0)
	    fNEvents[round(MN)][fCoupling] = 0;
	  fNEvents[round(MN)][fCoupling]++;
	
	  if (IsHNLGood == true && isGood == true)
	    fSumGood[round(MN)][fCoupling] += Weight;

	  // Histos for acceptance and yield: coupling

	  if (round(MN)/1000. == fMassForSingleValue) {
	    FillHisto("CouplingScan/hA", fCoupling, Weight*scale);
	    FillHisto("CouplingScan/hN", fCoupling, scale);

	    if (Weights[i]["ProdProb"] == fDBeProdProb) {
	      FillHisto("CouplingScan/hA1", fCoupling, Weight*scale);
	      FillHisto("CouplingScan/hN1", fCoupling, scale);
	    }
	    else if (Weights[i]["ProdProb"] == fDCuProdProb) {
	      FillHisto("CouplingScan/hA2", fCoupling, Weight*scale);
	      FillHisto("CouplingScan/hN2", fCoupling, scale);
	    }

	    if (IsHNLGood == true && isGood == true) {
	      FillHisto("CouplingScan/hG", fCoupling, Weight*scale);

	      if (Weights[i]["ProdProb"] == fDBeProdProb) { 
		FillHisto("CouplingScan/hG1", fCoupling, Weight*scale);
	      }
	      else if (Weights[i]["ProdProb"] == fDCuProdProb) {
		FillHisto("CouplingScan/hG2", fCoupling, Weight*scale);
	      }
	    }

	    if (isGood == true) {
	      FillHisto("CouplingScan/hR", fCoupling, Weight*scale);
	      if (Weights[i]["ProdProb"] == fDBeProdProb) {
                FillHisto("CouplingScan/hR1", fCoupling, Weight*scale);
              }
              else if (Weights[i]["ProdProb"] == fDCuProdProb) {
                FillHisto("CouplingScan/hR2", fCoupling, Weight*scale);
              }
	    }

	    // Histos for weight components: coupling

	    FillHisto("CouplingScan/hReachCoupling", fCoupling, NReachProb);
	    FillHisto("CouplingScan/hDecayCoupling", fCoupling, NDecayProb);
	    FillHisto("CouplingScan/hProbCoupling", fCoupling, NReachProb*NDecayProb);
	    FillHisto("CouplingScan/hWeightCoupling", fCoupling, Weight);	
	    FillHisto("CouplingScan/hEnergyCoupling", fCoupling, momN/1000.*momN/1000. + MN/1000.*MN/1000.);
	  }

	  // Histos for acceptance and yield: mass

	  if (fCoupling == fCouplingForSingleValue) {
	    FillHisto("MassScan/hA", round(MN)/1000., Weight*scale);
	    FillHisto("MassScan/hN", round(MN)/1000., scale);
	  
	    if (Weights[i]["ProdProb"] == fDBeProdProb) {
	      FillHisto("MassScan/hA1", round(MN)/1000., Weight*scale);
	      FillHisto("MassScan/hN1", round(MN)/1000., scale);
	    }
	    else if (Weights[i]["ProdProb"] == fDCuProdProb) {
	      FillHisto("MassScan/hA2", round(MN)/1000., Weight*scale);
	      FillHisto("MassScan/hN2", round(MN)/1000., scale);
	    }
	  
	    if (IsHNLGood == true && isGood == true) {
	      FillHisto("MassScan/hG", round(MN)/1000., Weight*scale);
	      if (Weights[i]["ProdProb"] == fDBeProdProb) {
		FillHisto("MassScan/hG1", round(MN)/1000., Weight*scale);
	      }
	      else if (Weights[i]["ProdProb"] == fDCuProdProb) {
		FillHisto("MassScan/hG2", round(MN)/1000., Weight*scale);
	      }
	    }

	    if (isGood == true) {
	      FillHisto("MassScan/hR", round(MN)/1000., Weight*scale);
	      if (Weights[i]["ProdProb"] == fDBeProdProb) {
                FillHisto("MassScan/hR1", round(MN)/1000., Weight*scale);
              }
              else if (Weights[i]["ProdProb"] == fDCuProdProb) {
                FillHisto("MassScan/hR2", round(MN)/1000., Weight*scale);
              }
	    }

	    // Histos for weight components: mass

	    FillHisto("MassScan/hReachMass", round(MN)/1000., NReachProb);
	    FillHisto("MassScan/hDecayMass", round(MN)/1000., NDecayProb);
	    FillHisto("MassScan/hProbMass", round(MN)/1000., NReachProb*NDecayProb);
	    FillHisto("MassScan/hWeightMass", round(MN)/1000., Weight);
	    FillHisto("MassScan/hEnergyMass", round(MN)/1000., momN/1000.*momN/1000. + MN/1000.*MN/1000.);
	  }

	  // Histos for weight components: one value

	  if (fCoupling == fCouplingForSingleValue && round(MN)/1000. == fMassForSingleValue) {
	    FillHisto("SingleValue/hCoupling",     fCoupling);
	    FillHisto("SingleValue/hMass",         round(MN)/1000.);
	    FillHisto("SingleValue/hHNLReachProb", NReachProb);
	    FillHisto("SingleValue/hHNLDecayProb", NDecayProb);
	    FillHisto("SingleValue/hWeight",       Weight);
	  }
	}
      }
    }
  
    // Scan on the N momentum
  
    if (GetWithMC()) {
      Event *evt = GetMCEvent();

      fUSquared = TMath::Power(10., fCouplingForSingleValue);
      std::vector<std::map<std::string, Double_t>> Weights = ComputeWeight(evt, fUSquared, fInitialUeSquaredRatio, fInitialUmuSquaredRatio, fInitialUtauSquaredRatio, fLInitialFV, fLFV, fMode);

      for (UInt_t i = 0; i < Weights.size(); i++) {
	MN =  round(Weights[i]["Mass"]);
	isGood = Weights[i]["IsGood"];

	if (MN/1000. == fMassForSingleValue) {
	  momN = Weights[i]["Momentum"]/1000.;
	  Weight = Weights[i]["Weight"];

	  if (momN >= fMomStart && momN <= fMomStop) {
	    momBin = fMomStep*trunc(momN/fMomStep);
	    FillHisto("SingleValue/hA", momBin, Weight*scale);
	    FillHisto("SingleValue/hN", momBin, scale);

	    if (IsHNLGood == true && isGood == true) {
	      FillHisto("SingleValue/hG", momBin, Weight*scale);
	    }

	    if (isGood == true) {
	      FillHisto("SingleValue/hR", momBin, Weight*scale);
	    }
	  }
	}
      }
    }
  
    // Some plots from MC

    Int_t kineCounter = 0;
    Double_t goodWeight = 0.;
    Double_t WeightTom = 0.;
    Double_t WeightGaia = 0.;
    Double_t goodWeightGaia = 0.;

    if (GetWithMC()) {
      Event *evt = GetMCEvent();

      // Weights for comparison

      std::vector<std::map<std::string, Double_t>> Weights = ComputeWeight(evt, fUSquared, fInitialUeSquaredRatio, fInitialUmuSquaredRatio, fInitialUtauSquaredRatio, fLInitialFV, fLFV, fMode);  
      std::vector<std::map<std::string, Double_t>> WeightsTom = ComputeWeight(evt, 1.E-6, 1., 1., 0., fLInitialFV, fLFV, fMode);
      std::vector<std::map<std::string, Double_t>> WeightsGaia = ComputeWeight(evt, 1.E-6, 0., 1., 0., fLInitialFV, fLFV, fMode);

      for (UInt_t i = 0; i < Weights.size(); i++) {
        if (Weights[i]["IsGood"] == 1)
	  goodWeight = Weights[i]["Weight"];
      }

      for (UInt_t i = 0; i < WeightsGaia.size(); i++) {
	if (WeightsGaia[i]["IsGood"] == 1)
          goodWeightGaia = WeightsGaia[i]["Weight"];
      }

      for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
	KinePart *p = evt->GetKinePart(i);

	if (p->GetParentID() == 0) {
	  Double_t Z = GeometricAcceptance::GetInstance()->GetZStraw(0);
	  Double_t X = p->GetProdPos().X() + p->GetInitial4Momentum().X()/p->GetInitial4Momentum().Z()*(Z - p->GetProdPos().Z());
	  Double_t Y = p->GetProdPos().Y() + p->GetInitial4Momentum().Y()/p->GetInitial4Momentum().Z()*(Z - p->GetProdPos().Z());
	  
	if (p->GetPDGcode() == 13 || p->GetPDGcode() == -13)
	  FillHisto("SingleValue/hXYMuon", X/1000., Y/1000., Weight);
	else if (p->GetPDGcode() == 211 || p->GetPDGcode() == -211)
	  FillHisto("SingleValue/hXYPion", X/1000., Y/1000., Weight);
	}

	if (p->GetParentID() == -1 && p->GetPDGcode() == 999) {
	  Weight = Weights[kineCounter]["Weight"];
	  WeightTom = WeightsTom[kineCounter]["Weight"];
	  WeightGaia = WeightsGaia[kineCounter]["Weight"];

	  FillHisto("SingleValue/hZMotherProd", p->GetPosAtCheckPoint(0).z()/1000.);

	  if (p->GetPosAtCheckPoint(1).z() != 0.)
	    FillHisto("SingleValue/hZDProd", p->GetPosAtCheckPoint(1).z()/1000.);
	  if (p->GetPosAtCheckPoint(2).z() != 0.)
	    FillHisto("SingleValue/hZTauProd", p->GetPosAtCheckPoint(2).z()/1000.);

	  FillHisto("SingleValue/hZDDecay", p->GetProdPos().Z()/1000.);
	  FillHisto("SingleValue/hDTheta", p->GetPosAtCheckPoint(1).x());
	  FillHisto("SingleValue/hDLambda", p->GetPosAtCheckPoint(1).y());
	  FillHisto("SingleValue/hDPath", p->GetMomAtCheckPoint(1).X());
	  FillHisto("SingleValue/hDMom", p->GetMomAtCheckPoint(1).Y()/1000.);
	  FillHisto("SingleValue/hZHNLDecay", p->GetEndPos().Z()/1000.);
	  FillHisto("SingleValue/hHNLGamma", p->GetInitial4Momentum().Gamma());
	  FillHisto("SingleValue/hHNLTheta", p->GetMomAtCheckPoint(0).Z());
	  FillHisto("SingleValue/hHNLMom", p->GetMomAtCheckPoint(0).T()/1000.);

	  if (p->GetEndProcessName() == "good") {

	    fDistcomp->SetLineDir(TVector3(p->GetInitial4Momentum().X(), p->GetInitial4Momentum().Y(), p->GetInitial4Momentum().Z()));
	    fDistcomp->SetLinePoint1(TVector3(p->GetEndPos().X(), p->GetEndPos().Y(), p->GetEndPos().Z()));
	    fDistcomp->SetPoint(0., 0., -26.5); // Z = mean of Z of N prod point             
	    fDistcomp->ComputeDistance();

	    Double_t TargetDist = fDistcomp->GetDistance();

	    fDistcomp->SetLineDir(TVector3(p->GetInitial4Momentum().X(), p->GetInitial4Momentum().Y(), p->GetInitial4Momentum().Z()));
	    fDistcomp->SetLinePoint1(TVector3(p->GetEndPos().X(), p->GetEndPos().Y(), p->GetEndPos().Z()));
	    fDistcomp->SetPoint(0., -22., 23230.);
	    fDistcomp->ComputeDistance();

	    Double_t TAXDist = fDistcomp->GetDistance();

	    fDistcomp->SetLinePoint1(0., 0., 102000.);
	    fDistcomp->SetLineDir(0., 0., 1.);
	    fDistcomp->SetPoint(TVector3(p->GetEndPos().X(), p->GetEndPos().Y(), p->GetEndPos().Z()));

	    fDistcomp->ComputeDistance();

	    Double_t BeamlineDist = fDistcomp->GetDistance();


	    FillHisto("SingleValue/hBeamvsTarTrue", TargetDist/1000., BeamlineDist/1000., goodWeight);
	    FillHisto("SingleValue/hBeamvsTAXTrue", TAXDist/1000., BeamlineDist/1000., goodWeight);
	  }

	  // Tommaso comparison

	  if (TMath::Sqrt(p->GetInitial4Momentum().E()*p->GetInitial4Momentum().E() - p->GetInitial4Momentum().P()*p->GetInitial4Momentum().P())/1000. == 0.8) {
	    if (p->GetParticleName() == "D+->e+N" || p->GetParticleName() == "D-->e-N" || p->GetParticleName() == "D+->mu+N" || p->GetParticleName() == "D-->mu-N" || p->GetParticleName() == "DS+->e+N" || p->GetParticleName() == "DS-->e-N"|| p->GetParticleName() == "DS+->mu+N" || p->GetParticleName() == "DS-->mu-N") {
	      if (p->GetProdPos().Z() <= 400. && p->GetProdPos().Z() >= -400.) {
		FillHisto("SingleValue/hDThetaMom", p->GetMomAtCheckPoint(1).Y()/1000., p->GetPosAtCheckPoint(1).x(), WeightTom);
		FillHisto("SingleValue/hXYDecay", p->GetEndPos().X()/1000., p->GetEndPos().Y()/1000., WeightTom);
		FillHisto("SingleValue/hTransverse", TMath::Sqrt(p->GetEndPos().X()/1000.*p->GetEndPos().X()/1000. + p->GetEndPos().Y()/1000.*p->GetEndPos().Y()/1000.), WeightTom);
	      }
	    }
	  }
	  
	  // Toy-MC comparison: all
	  
	  if ((TMath::Sqrt(p->GetInitial4Momentum().E()*p->GetInitial4Momentum().E() - p->GetInitial4Momentum().P()*p->GetInitial4Momentum().P()))/1000. == fMassForSingleValue) {
	    if (p->GetMomAtCheckPoint(3).X() != 0. && (p->GetParticleName() == "DS+->mu+N" || p->GetParticleName() == "DS-->mu-N")) {
	      FillHisto("ToyMC/DS/hpDS", p->GetMomAtCheckPoint(3).X()/1000., WeightGaia);
	      FillHisto("ToyMC/DS/hptDS", p->GetMomAtCheckPoint(3).Y()/1000., WeightGaia);
	      FillHisto("ToyMC/DS/hpNS", p->GetMomAtCheckPoint(3).Z()/1000., WeightGaia);
	      FillHisto("ToyMC/DS/hptNS", p->GetMomAtCheckPoint(3).T()/1000., WeightGaia);
	      FillHisto("ToyMC/DS/hpmuS", p->GetMomAtCheckPoint(4).X()/1000., WeightGaia);
	      FillHisto("ToyMC/DS/hptmuS", p->GetMomAtCheckPoint(4).Y()/1000., WeightGaia);
	    }
	    
	    if (p->GetMomAtCheckPoint(5).X() != 0. && (p->GetParticleName() == "D0->mu+K-N" || p->GetParticleName() == "D0->mu-K+N" || p->GetParticleName() == "D0Bar->mu+K-N" || p->GetParticleName() == "D0Bar->mu-K+N")) {
	      FillHisto("ToyMC/D0/hpD0", p->GetMomAtCheckPoint(5).X()/1000., WeightGaia);
	      FillHisto("ToyMC/D0/hptD0", p->GetMomAtCheckPoint(5).Y()/1000., WeightGaia);
	      FillHisto("ToyMC/D0/hpN0", p->GetMomAtCheckPoint(5).Z()/1000., WeightGaia);
	      FillHisto("ToyMC/D0/hptN0", p->GetMomAtCheckPoint(5).T()/1000., WeightGaia);
	      FillHisto("ToyMC/D0/hpmu0", p->GetMomAtCheckPoint(6).X()/1000., WeightGaia);
	      FillHisto("ToyMC/D0/hptmu0", p->GetMomAtCheckPoint(6).Y()/1000., WeightGaia);
	    }

	    // Toy-MC comparison: good

	    if (p->GetEndProcessName() == "good") {
	      if (p->GetMomAtCheckPoint(3).X() != 0.) {
		FillHisto("ToyMC/DS/hpDS1", p->GetMomAtCheckPoint(3).X()/1000., goodWeightGaia);
		FillHisto("ToyMC/DS/hptDS1", p->GetMomAtCheckPoint(3).Y()/1000., goodWeightGaia);
		FillHisto("ToyMC/DS/hpNS1", p->GetMomAtCheckPoint(3).Z()/1000., goodWeightGaia);
		FillHisto("ToyMC/DS/hptNS1", p->GetMomAtCheckPoint(3).T()/1000., goodWeightGaia);
		FillHisto("ToyMC/DS/hpmuS1", p->GetMomAtCheckPoint(4).X()/1000., goodWeightGaia);
		FillHisto("ToyMC/DS/hptmuS1", p->GetMomAtCheckPoint(4).Y()/1000., goodWeightGaia);
	      }

	      if (p->GetMomAtCheckPoint(5).X() != 0.) {
		FillHisto("ToyMC/D0/hpD01", p->GetMomAtCheckPoint(5).X()/1000., goodWeightGaia);
		FillHisto("ToyMC/D0/hptD01", p->GetMomAtCheckPoint(5).Y()/1000., goodWeightGaia);
		FillHisto("ToyMC/D0/hpN01", p->GetMomAtCheckPoint(5).Z()/1000., goodWeightGaia);
		FillHisto("ToyMC/D0/hptN01", p->GetMomAtCheckPoint(5).T()/1000., goodWeightGaia);
		FillHisto("ToyMC/D0/hpmu01", p->GetMomAtCheckPoint(6).X()/1000., goodWeightGaia);
		FillHisto("ToyMC/D0/hptmu01", p->GetMomAtCheckPoint(6).Y()/1000., goodWeightGaia);
	      }
	    }
	  }
	}
	kineCounter++;
      }
    }
  }
}

void HeavyNeutrinoScan::EndOfJobUser() {
 
  if (fReadingData) {

    // X axis
 
    fHisto.GetTH1("SingleValue/hZMotherProd")->GetXaxis()->SetTitle("Position along Z [m]");
    fHisto.GetTH1("SingleValue/hZDProd")->GetXaxis()->SetTitle("Position along Z [m]");
    fHisto.GetTH1("SingleValue/hZTauProd")->GetXaxis()->SetTitle("Position along Z [m]");
    fHisto.GetTH1("SingleValue/hZDDecay")->GetXaxis()->SetTitle("Position along Z [m]");
    fHisto.GetTH1("SingleValue/hDTheta")->GetXaxis()->SetTitle("Polar angle [rad]");
    fHisto.GetTH1("SingleValue/hDLambda")->GetXaxis()->SetTitle("Decay length [mm]");
    fHisto.GetTH1("SingleValue/hDPath")->GetXaxis()->SetTitle("Z [mm]");
    fHisto.GetTH1("SingleValue/hDMom")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("SingleValue/hZHNLDecay")->GetXaxis()->SetTitle("Position along Z [m]");
    fHisto.GetTH1("SingleValue/hHNLGamma")->GetXaxis()->SetTitle("Lorentz factor");
    fHisto.GetTH1("SingleValue/hHNLDecayProb")->GetXaxis()->SetTitle("Decay probability");
    fHisto.GetTH1("SingleValue/hHNLReachProb")->GetXaxis()->SetTitle("Reach probability");
    fHisto.GetTH1("SingleValue/hHNLTheta")->GetXaxis()->SetTitle("Polar angle [rad]");
    fHisto.GetTH1("SingleValue/hHNLMom")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("SingleValue/hWeight")->GetXaxis()->SetTitle("Weight");
    fHisto.GetTH1("SingleValue/hCoupling")->GetXaxis()->SetTitle("Coupling");
    fHisto.GetTH1("SingleValue/hMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");  
    fHisto.GetTH1("SingleValue/hTransverse")->GetXaxis()->SetTitle("Distance [m]");
    fHisto.GetTH2("SingleValue/hDThetaMom")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH2("SingleValue/hXYDecay")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("SingleValue/hXYMuon")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("SingleValue/hXYPion")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("SingleValue/hBeamvsTarTrue")->GetXaxis()->SetTitle("Impact parameter [m]");
    fHisto.GetTH2("SingleValue/hBeamvsTAXTrue")->GetXaxis()->SetTitle("Impact parameter [m]");
    fHisto.GetTH2("CouplingScan/hReachCoupling")->GetXaxis()->SetTitle("Log(U^{2})");
    fHisto.GetTH2("CouplingScan/hDecayCoupling")->GetXaxis()->SetTitle("Log(U^{2})");
    fHisto.GetTH2("CouplingScan/hProbCoupling")->GetXaxis()->SetTitle("Log(U^{2})");
    fHisto.GetTH2("CouplingScan/hWeightCoupling")->GetXaxis()->SetTitle("Log(U^{2})");
    fHisto.GetTH2("CouplingScan/hEnergyCoupling")->GetXaxis()->SetTitle("Log(U^{2})");
    fHisto.GetTH2("MassScan/hReachMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("MassScan/hDecayMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("MassScan/hProbMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("MassScan/hWeightMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("MassScan/hEnergyMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH1("ToyMC/DS/hpDS")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/DS/hpDS1")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/DS/hptDS")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/DS/hptDS1")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/DS/hpNS")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/DS/hpNS1")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/DS/hptNS")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/DS/hptNS1")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/DS/hpmuS")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/DS/hpmuS1")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/DS/hptmuS")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/DS/hptmuS1")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hpD0")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hpD01")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hptD0")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hptD01")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hpN0")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hpN01")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hptN0")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hptN01")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hpmu0")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hpmu01")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hptmu0")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("ToyMC/D0/hptmu01")->GetXaxis()->SetTitle("P [GeV/c]");

    // Y axis

    fHisto.GetTH2("SingleValue/hDThetaMom")->GetYaxis()->SetTitle("Polar angle [rad]");
    fHisto.GetTH2("SingleValue/hXYDecay")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("SingleValue/hXYMuon")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("SingleValue/hXYPion")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("SingleValue/hBeamvsTarTrue")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
    fHisto.GetTH2("SingleValue/hBeamvsTAXTrue")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
    fHisto.GetTH2("CouplingScan/hReachCoupling")->GetYaxis()->SetTitle("Reach probability");
    fHisto.GetTH2("CouplingScan/hDecayCoupling")->GetYaxis()->SetTitle("Decay probability");
    fHisto.GetTH2("CouplingScan/hProbCoupling")->GetYaxis()->SetTitle("Reach and decay probability");
    fHisto.GetTH2("CouplingScan/hWeightCoupling")->GetYaxis()->SetTitle("Weight");
    fHisto.GetTH2("CouplingScan/hEnergyCoupling")->GetYaxis()->SetTitle("N energy [GeV]");
    fHisto.GetTH2("MassScan/hReachMass")->GetYaxis()->SetTitle("Reach probability");
    fHisto.GetTH2("MassScan/hDecayMass")->GetYaxis()->SetTitle("Decay probability");
    fHisto.GetTH2("MassScan/hProbMass")->GetYaxis()->SetTitle("Reach and decay probability");
    fHisto.GetTH2("MassScan/hWeightMass")->GetYaxis()->SetTitle("Weight");
    fHisto.GetTH2("MassScan/hEnergyMass")->GetYaxis()->SetTitle("N energy [GeV]");
    fHisto.GetTH1("ToyMC/DS/hpDS")->GetYaxis()->SetTitle("Entries/5 GeV/c");
    fHisto.GetTH1("ToyMC/DS/hpDS1")->GetYaxis()->SetTitle("Entries/8 GeV/c");
    fHisto.GetTH1("ToyMC/DS/hptDS")->GetYaxis()->SetTitle("Entries/0.5 GeV/c");
    fHisto.GetTH1("ToyMC/DS/hptDS1")->GetYaxis()->SetTitle("Entries/0.5 GeV/c");
    fHisto.GetTH1("ToyMC/DS/hpNS")->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fHisto.GetTH1("ToyMC/DS/hpNS1")->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fHisto.GetTH1("ToyMC/DS/hptNS")->GetYaxis()->SetTitle("Entries/0.1 GeV/c");
    fHisto.GetTH1("ToyMC/DS/hptNS1")->GetYaxis()->SetTitle("Entries/0.1 GeV/c");
    fHisto.GetTH1("ToyMC/DS/hpmuS")->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fHisto.GetTH1("ToyMC/DS/hpmuS1")->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fHisto.GetTH1("ToyMC/DS/hptmuS")->GetYaxis()->SetTitle("Entries/0.1 GeV/c");
    fHisto.GetTH1("ToyMC/DS/hptmuS1")->GetYaxis()->SetTitle("Entries/0.1 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hpD0")->GetYaxis()->SetTitle("Entries/5 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hpD01")->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hptD0")->GetYaxis()->SetTitle("Entries/0.5 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hptD01")->GetYaxis()->SetTitle("Entries/0.5 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hpN0")->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hpN01")->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hptN0")->GetYaxis()->SetTitle("Entries/0.1 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hptN01")->GetYaxis()->SetTitle("Entries/0.1 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hpmu0")->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hpmu01")->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hptmu0")->GetYaxis()->SetTitle("Entries/0.1 GeV/c");
    fHisto.GetTH1("ToyMC/D0/hptmu01")->GetYaxis()->SetTitle("Entries/0.1 GeV/c");

    Double_t Coupling = 0.;
    Double_t MN = 0.;

    // Write to files for exclusion plot

    for (auto it = fMasses.begin(); it != fMasses.end(); it++) {
      MN = it->first;
      for (auto it1 = fCouplings.begin(); it1 != fCouplings.end(); it1++) {
	Coupling = it1->first;
	fGammaTotFile << MN << " " << Coupling << " " << fGammaTot[MN][Coupling] << "\n";
	fTauFile << MN << " " << Coupling << " " << fTau[MN][Coupling] << "\n";
	fNEventsFile << MN << " " << Coupling << " " << fNEvents[MN][Coupling] << "\n";
	fSumGoodFile << MN << " " << Coupling << " " << fSumGood[MN][Coupling] << "\n";
      }
    }

    fGammaTotFile.close();
    fTauFile.close();
    fNEventsFile.close();
    fSumGoodFile.close();
  }
  else {

    // Acceptance and yield computation: momentum

    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomSel"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sG", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sR", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomReg"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sR", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sA", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomFV"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sA", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sN", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorYieldMom"))->Divide((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sG", true), (TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sN", true), "cl=0.683 mode");
  
    // Acceptance and yield computation: coupling

    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSelTarget"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cG1", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cR1", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingRegTarget"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cR1", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cA1", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFVTarget"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cA1", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cN1", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTarget"))->Divide((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cG1", true), (TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cN1", true), "cl=0.683 mode");

    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSelTAX"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cG2", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cR2", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingRegTAX"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cR2", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cA2", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFVTAX"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cA2", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cN2", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTAX"))->Divide((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cG2", true), (TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cN2", true), "cl=0.683 mode");

    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSel"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cG", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cR", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingReg"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cR", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cA", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFV"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cA", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cN", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCoupling"))->Divide((TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cG", true), (TH1D*)RequestHistogram(fAnalyzerName, "CouplingScan/cN", true), "cl=0.683 mode");


    // Acceptance and yield computation: mass

    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassSel"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mG", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassReg"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassFV"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mN", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMass"))->Divide((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mG", true), (TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mN", true), "cl=0.683 mode");

    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassSelTarget"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mG1", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR1", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassRegTarget"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR1", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA1", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassFVTarget"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA1", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mN1", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTarget"))->Divide((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mG1", true), (TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mN1", true), "cl=0.683 mode");

    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassSelTAX"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mG2", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR2", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassRegTAX"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR2", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA2", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassFVTAX"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA2", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mN2", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTAX"))->Divide((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mG2", true), (TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mN2", true), "cl=0.683 mode");

    // Retrieve file info for exclusion plot
    
    Double_t Coupling = 0.;
    Double_t MN = 0.;
    Double_t POT = 1.E18;
    Int_t couplingCounter = 0;
    Int_t massCounter = 0;
    string line;

    while (getline(fGammaTotFile, line)) {
      stringstream ss(line);
      Double_t value;
      ss >> MN >> Coupling >> value;
      fGammaTot[MN][Coupling] = value;
    }

    while (getline(fTauFile, line)) {
      stringstream ss(line);
      Double_t value;
      ss >> MN >> Coupling >> value;
      fTau[MN][Coupling] = value;
    }

    while (getline(fNEventsFile, line)) {
      stringstream ss(line);
      Double_t value;
      ss >> MN >> Coupling >> value;
      fNEvents[MN][Coupling] = value;
    }

    while (getline(fSumGoodFile, line)) {
      stringstream ss(line);
      Double_t value;
      ss >> MN >> Coupling >> value;
      fSumGood[MN][Coupling] = value;
    }

    // Lifetime and exclusion plots

    for (Int_t i = fMassStart*1000.; i <= fMassStop*1000.; i += fMassStep*1000) {
      MN = i;
      for (Int_t j  = fCouplingStart*10; j <= fCouplingStop*10; j += fCouplingStep*10) {
	Coupling = j/10.;

	if (MN/1000. == fMassForSingleValue) {
	  fHisto.GetTGraph("CouplingScan/GammaTotCoupling")->SetPoint(couplingCounter, Coupling, fGammaTot[MN][Coupling]);
	  fHisto.GetTGraph("CouplingScan/TauCoupling")->SetPoint(couplingCounter, Coupling, fTau[MN][Coupling]);
	  couplingCounter++;
	}

	fYield[MN][Coupling] = fSumGood[MN][Coupling]/fNEvents[MN][Coupling];
	FillHisto("TotalScan/hExclusion", MN/1000., Coupling, fYield[MN][Coupling]*POT);
      }
    
      Coupling = fCouplingForSingleValue;
      fHisto.GetTGraph("MassScan/GammaTotMass")->SetPoint(massCounter, MN/1000., fGammaTot[MN][Coupling]);
      fHisto.GetTGraph("MassScan/TauMass")->SetPoint(massCounter, MN/1000., fTau[MN][Coupling]);
      massCounter++;
    }

    // Sensitivity

    for (int i=1; i<=fHisto.GetTH2("TotalScan/hExclusion")->GetNbinsX(); i++) {
      for (int j=1; j<=fHisto.GetTH2("TotalScan/hExclusion")->GetNbinsY(); j++) {
	cout<<fHisto.GetTH2("TotalScan/hExclusion")->GetXaxis()->GetBinCenter(i)<<" "<<fHisto.GetTH2("TotalScan/hExclusion")->GetYaxis()->GetBinCenter(j)<<" "<<fHisto.GetTH2("TotalScan/hExclusion")->GetBinContent(i,j)<<endl;
      }
    }
    cout<<"--------------------------"<<endl;
    EvaluateUL(fHisto.GetTH2("TotalScan/hExclusion"), fHisto.GetTGraph("TotalScan/Contours"));

    // Mean probability vs mass and coupling

    for (Int_t i = 0; i < static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "CouplingScan/ProbCoupling", true))->GetNbinsX()+1; i++) {
      TH1D *h = static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "CouplingScan/ProbCoupling", true))->ProjectionY("", i, i+1, "");
      fHisto.GetTGraph("CouplingScan/MeanCoupling")->SetPoint(i, static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "CouplingScan/ProbCoupling", true))->GetXaxis()->GetBinCenter(i), h->GetMean());
    }
  
    for (Int_t i = 0; i < static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "MassScan/ProbMass", true))->GetNbinsX()+1; i++) {      
      TH1D *h = static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "MassScan/ProbMass", true))->ProjectionY("", i, i+1, "");           
      fHisto.GetTGraph("MassScan/MeanMass")->SetPoint(i, static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "MassScan/ProbMass", true))->GetXaxis()->GetBinCenter(i), h->GetMean());
    }  
    
    // Cosmetics
  
    CosmeticsGraph(fHisto.GetTGraph("CouplingScan/GammaTotCoupling"), "Log(U^{2})", "Total decay width [MeV]", 2);
    CosmeticsGraph(fHisto.GetTGraph("CouplingScan/TauCoupling"), "Log(U^{2})", "Lifetime [ns]", 2);
    CosmeticsGraph(fHisto.GetTGraph("CouplingScan/MeanCoupling"), "Log(U^{2})", "Mean probability", 2);
    CosmeticsGraph(fHisto.GetTGraph("MassScan/GammaTotMass"), "N mass [GeV/c^{2}]", "Total decay width [MeV]", 2);
    CosmeticsGraph(fHisto.GetTGraph("MassScan/TauMass"), "N mass [GeV/c^{2}]", "Lifetime [ns]", 2);
    CosmeticsGraph(fHisto.GetTGraph("MassScan/MeanMass"), "N mass [GeV/c^{2}]", "Mean probability", 2);
    CosmeticsGraph(fHisto.GetTGraph("TotalScan/Contours"), "N mass [GeV/c^{2}]", "Log(U^{2})", 2);  

    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSel")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingReg")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFV")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCoupling")), "Log(U^{2})", "Yield per POT", 2);

    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSelTarget")), "Log(U^{2})", "Acceptance", 9);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingRegTarget")), "Log(U^{2})", "Acceptance", 9);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFVTarget")), "Log(U^{2})", "Acceptance", 9);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTarget")), "Log(U^{2})", "Yield per POT", 9);

    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSelTAX")), "Log(U^{2})", "Acceptance", 8);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingRegTAX")), "Log(U^{2})", "Acceptance", 8);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFVTAX")), "Log(U^{2})", "Acceptance", 8);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTAX")), "Log(U^{2})", "Yield per POT", 8);

    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassSel")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassReg")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassFV")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMass")), "Log(U^{2})", "Yield per POT", 2);

    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassSelTarget")), "Log(U^{2})", "Acceptance", 9);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassRegTarget")), "Log(U^{2})", "Acceptance", 9);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassFVTarget")), "Log(U^{2})", "Acceptance", 9);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTarget")), "Log(U^{2})", "Yield per POT", 9);

    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassSelTAX")), "Log(U^{2})", "Acceptance", 8);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassRegTAX")), "Log(U^{2})", "Acceptance", 8);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassFVTAX")), "Log(U^{2})", "Acceptance", 8);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTAX")), "Log(U^{2})", "Yield per POT", 8);

    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomSel")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomReg")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomFV")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorYieldMom")), "Log(U^{2})", "Yield per POT", 2);

    fHisto.GetTH2("TotalScan/hExclusion")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("TotalScan/hExclusion")->GetYaxis()->SetTitle("Log(U^{2})");
  }
  
  //remove("GammaTotFile.txt"); 
  //remove("TauFile.txt"); 
  //remove("NEventsFile.txt"); 
  //remove("SumGoodFile.txt"); 

  SaveAllPlots();
  
  return;
}

void HeavyNeutrinoScan::CosmeticsGraph(TGraph* g, const char* x, const char* y, Int_t col) {

  g->GetXaxis()->SetTitle(x);
  g->GetYaxis()->SetTitle(y);
  g->SetLineColor(col);
  g->SetLineWidth(3);
  g->Draw("AC");
  
  return;
}

void HeavyNeutrinoScan::CosmeticsGraph(TGraphAsymmErrors* g, const char* x, const char* y, Int_t col) {

  g->GetXaxis()->SetTitle(x);
  g->GetYaxis()->SetTitle(y);
  g->SetLineColor(1);
  g->SetLineWidth(1);
  g->SetMarkerStyle(8);
  g->SetMarkerColor(col);
  g->Draw("AP");
  
  return;
}

void HeavyNeutrinoScan::EvaluateUL(TH2* h, TGraph* gr) {

  Int_t nContours = 6;
  Double_t* CLs = new Double_t[nContours];

  CLs[0] = 0.6826895; // 1-sigma contour                              
  CLs[1] = 0.90; // 90% CL exclusion contours  
  CLs[2] = 0.9544997; // 2-sigma contour       
  CLs[3] = 0.9973002; // 3-sigma contour                
  CLs[4] = 0.9999366; // 4-sigma contour                
  CLs[5] = 0.9999994; // 5-sigma contour       

  Double_t* contours = new Double_t[nContours];

  for (Int_t i = 0; i < nContours; i++) {
    Double_t sigLevel = -TMath::Log(1-CLs[i]); // maximum number of events, given the CL, for Nobs=0 and Bkg=0.                     

    contours[i] = sigLevel;
  }

  h->SetContour(nContours, contours);

  std::vector<TGraph*> contoursHisto = ExtractContours(h);

  if (contoursHisto.size() == 0) {
    cout<<"Graf for UL is empty"<<endl;
    return;
  }

  Int_t nn = contoursHisto[1]->GetN();

  TCanvas *c = new TCanvas();
  contoursHisto[1]->Draw();
  c->Print("cont.pdf");

  Double_t* xx = new Double_t[nn];
  Double_t* yy = new Double_t[nn];

  xx = contoursHisto[1]->GetX();
  yy = contoursHisto[1]->GetY();

  for (Int_t i = 0; i < nn; i++) {
    gr->SetPoint(i, xx[i], yy[i]);
  }
  
  gr->SetLineColor(kRed);
  gr->SetLineWidth(3);
  gr->Draw("AP*");
  
  return;
}
  
std::vector<TGraph*> HeavyNeutrinoScan::ExtractContours(TH2* h) {
  
  // Draw contours as filled regions and save points 
  
  TCanvas *c = new TCanvas();
  h->Draw("cont z list");
  c->Print("minc.pdf");

  TObjArray* conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
  TList* contLevel = NULL;
  Int_t TotalConts = 0;
  std::vector<TGraph*> graf;

  if (conts == NULL){
    TotalConts = 0;
    return graf;
  } 
  else {
    TotalConts = conts->GetSize();
  }

  for (Int_t i = 0; i < TotalConts; i++) {
    contLevel = (TList*)conts->At(i);
    if (contLevel) {
      graf.push_back(new TGraph());
      graf[i]->Merge(contLevel);
      graf[i]->SetName(Form("c%d_i%d", h->GetUniqueID(), i));
    }
  }

  //remove("minc.pdf");
  //remove("cont.pdf");

  return graf;
}
