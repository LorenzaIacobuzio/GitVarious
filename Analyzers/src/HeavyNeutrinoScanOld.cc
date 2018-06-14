// ---------------------------------------------------------------                                      
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
/// \Detailed                                                                                           
/// After retrieving all HNL weights using the tool called HNLWeight, and checking if an event passed the selection implemented in the analyzer HeavyNeutrino, acceptance and yield per POT are computed. 
/// Different sets of plots are produced, all together at once:
/// for one value of HNL mass and coupling (that can be set as external parameters; default values are MassForSingleValue = 1000 MeV and CouplingForSingleValue = -6 = Log(10^-6)), quantities related to the D meson, the HNL and its product are produced;acceptance and yield are also computed; the yield is computed as a function of the HNL momentum as well;
/// a scan is enabled, to produce plots as a function of the coupling;
/// if the MC sample contains several masses, plots are produced as a function of the HNL mass; in this case, an expected exclusion plot is produced as a function of the HNL mass and coupling. Thus, the analyzer is able to run either on MC samples of just one HNL mass or on samples containing HNLs of different masses;
/// For HNL masses of 1 GeV, several plots are produced for toy-MC comparison: momentum and transverse momentum of the D meson, HNL and muon daughter (both for all HNLs and only for the good HNL), for two decay chains, as produced in the MC routine: DS->Nmu; N->pimu, and D0->KNmu; N->pimu.
///                  
/// This analyzer makes use of two ToolsLib, called HNLFunctions and HNLWeight.  
/// The values of the ratios between specific-flavour couplings can be either set as external parameters from command line or taken as default values.               
/// For example, if the user sets USquared = 1.E-10, UeSquaredRatio = 5., UmuSquaredRatio = 1., UtauSquaredRatio = 3.5, the specific-flavour coupling values will be: UeSquared = 5.25E-11, UmuSquared = 1.05E-11, UtauSquared = 3.68E-11.   
/// The values of the initial and final coupling and the scan step can be either set as external parameters from command line or taken as default values.               
/// For example, if the user sets CouplingStart = -3, CouplingStop = -2. and CouplingStep = 0.5, the scan will be performed for Log(U2) = [-3., -2.5, -2.], this meaning U2 = [1.E-3., 1.E-2.5, 1.E-2.].                                                          
///
/// This analyzer produces .txt files needed for computing error bars for plots generated during the first step, and it can be run in --histo mode to produce these plots:
/// ./myExecutable -i outFile.root -o outFile1.root --histo
/// The files are deleted at the end of the second step.
/// The values of the initial and final momentum and its step can be either set as external parameters from command line or taken as default values. These are used to produce plots of the acceptance and yield per POT as a function of the HNL momentum (default values are MomStart = 0. GeV, MomStop = 150. GeV and MomStep = 5. GeV).
/// Other parameters to be set from command line are:
/// length of the beginning Z of the fiducial volume (default InitialFV = 102500. mm); length of the FV (default LFV = 77500. mm); decay mode of the HNLs, as set in the MC macro (default Mode = 0, which is pimu mode); number of events in each chunk the sample is split into, to compute error bars in --histo mode.                               
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
#include "MCInfo.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "BeamParameters.hh"
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
  AddParam("UeSquaredRatio", &fUeSquaredRatio, 1.);
  AddParam("UmuSquaredRatio", &fUmuSquaredRatio, 16.);
  AddParam("UtauSquaredRatio", &fUtauSquaredRatio, 3.8);
  AddParam("CouplingStart", &fCouplingStart, -10.);
  AddParam("CouplingStop", &fCouplingStop, 0.);
  AddParam("CouplingStep", &fCouplingStep, 1.);
  AddParam("InitialFV", &fInitialFV, 102500.);
  AddParam("LFV", &fLFV, 77500.);
  AddParam("Mode", &fMode, 0);
  AddParam("MomStop", &fMomStop, 150.);
  AddParam("MomStart", &fMomStart, 0.);
  AddParam("MomStep", &fMomStep, 5.);
  AddParam("MassForSingleValue", &fMassForSingleValue, 1000.);
  AddParam("CouplingForSingleValue", &fCouplingForSingleValue, -6.);
  AddParam("SplitStep", &fSplitStep, 100);

  fNMom = round((std::abs(fMomStop - fMomStart))/fMomStep);
  fN = round((std::abs(fCouplingStop-fCouplingStart))/fCouplingStep);
  fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
  fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
  fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;

  if (fReadingData) {
    
    fErrorFile.      open("ErrorBars.txt",       fstream::out | fstream::in | fstream::trunc);
    fErrorFileTarget.open("ErrorBarsTarget.txt", fstream::out | fstream::in | fstream::trunc);
    fErrorFileTAX.   open("ErrorBarsTAX.txt",    fstream::out | fstream::in | fstream::trunc);
    fErrorFileMom.   open("ErrorBarsMom.txt",    fstream::out | fstream::in | fstream::trunc);

    if (!fErrorFile.is_open() || !fErrorFileTarget.is_open() || !fErrorFileTAX.is_open() || !fErrorFileMom.is_open()) {
      cout << "Cannot open error bar files" << endl;
      _exit(1);
    }
   
    // One value histos
    
    fhZMotherProd  = nullptr;
    fhZDProd       = nullptr;
    fhZTauProd     = nullptr;
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
    fhCoupling     = nullptr;
    fhMass         = nullptr;
    fhAcc          = nullptr;
    fhAccTarget    = nullptr;
    fhAccTAX       = nullptr;
    fhYield        = nullptr;
    fhYieldTarget  = nullptr;
    fhYieldTAX     = nullptr;
    fhTransverse   = nullptr;
    fhDThetaMom    = nullptr;
    fhXYDecay      = nullptr;

    // Scan histos
    
    fhReachCoupling       = nullptr;
    fhDecayCoupling       = nullptr;
    fhProbCoupling        = nullptr;
    fhWeightCoupling      = nullptr;
    fhEnergyCoupling      = nullptr;
    fgAccCoupling         = nullptr;
    fgAccCouplingTarget   = nullptr;
    fgAccCouplingTAX      = nullptr;
    fgYieldCoupling       = nullptr;
    fgYieldCouplingTarget = nullptr;
    fgYieldCouplingTAX    = nullptr;
    fgGammaTotCoupling    = nullptr;
    fgTauCoupling         = nullptr;
    fhReachMass           = nullptr;
    fhDecayMass           = nullptr;
    fhProbMass            = nullptr;
    fhWeightMass          = nullptr;
    fhEnergyMass          = nullptr;
    fgAccMass             = nullptr;
    fgAccMassTarget       = nullptr;
    fgAccMassTAX          = nullptr;
    fgYieldMass           = nullptr;
    fgYieldMassTarget     = nullptr;
    fgYieldMassTAX        = nullptr;
    fgGammaTotMass        = nullptr;
    fgTauMass             = nullptr;
    fgAccMom              = nullptr;
    fgYieldMom            = nullptr;
    fgExclusion           = nullptr;

    // Toy-MC comparison

    fhpDS    = nullptr;
    fhpDS1   = nullptr;
    fhptDS   = nullptr;
    fhptDS1  = nullptr;
    fhpNS    = nullptr;
    fhpNS1   = nullptr;
    fhptNS   = nullptr;
    fhptNS1  = nullptr;
    fhpmuS   = nullptr;
    fhpmuS1  = nullptr;
    fhptmuS  = nullptr;
    fhptmuS1 = nullptr;
    fhpD0    = nullptr;
    fhpD01   = nullptr;
    fhptD0   = nullptr;
    fhptD01  = nullptr;
    fhpN0    = nullptr;
    fhpN01   = nullptr;
    fhptN0   = nullptr;
    fhptN01  = nullptr;
    fhpmu0   = nullptr;
    fhpmu01  = nullptr;
    fhptmu0  = nullptr;
    fhptmu01 = nullptr;
  }
  else {

    fErrorFile.      open("ErrorBars.txt",       fstream::out | fstream::in);
    fErrorFileTarget.open("ErrorBarsTarget.txt", fstream::out | fstream::in);
    fErrorFileTAX.   open("ErrorBarsTAX.txt",    fstream::out | fstream::in);
    fErrorFileMom.   open("ErrorBarsMom.txt",    fstream::out | fstream::in);

    if (!fErrorFile.is_open() || !fErrorFileTarget.is_open() || !fErrorFileTAX.is_open() || !fErrorFileMom.is_open()) {
      cout << "Cannot open error bar files" << endl;
      _exit(1);
    }

    fgErrorAccCoupling         = nullptr; 
    fgErrorAccCouplingTarget   = nullptr;
    fgErrorAccCouplingTAX      = nullptr;
    fgErrorYieldCoupling       = nullptr;
    fgErrorYieldCouplingTarget = nullptr;
    fgErrorYieldCouplingTAX    = nullptr;
    fgErrorAccMass             = nullptr;
    fgErrorAccMassTarget       = nullptr;
    fgErrorAccMassTAX          = nullptr;
    fgErrorYieldMass           = nullptr;
    fgErrorYieldMassTarget     = nullptr;
    fgErrorYieldMassTAX        = nullptr;
    fgErrorAccMom              = nullptr;
    fgErrorYieldMom            = nullptr;
    fgErrorExclusion           = nullptr;
  }
}

void HeavyNeutrinoScan::InitHist() {
  
  if (fReadingData) {
    
    // One value of mass and coupling
    
    BookHisto("SingleValue/hZMotherProd",  new TH1D("ZMotherProd", "N mother production point", 200, 0., 27.));
    BookHisto("SingleValue/hZDProd",       new TH1D("ZDProd",      "D meson production point",  200, 0., 27.));
    BookHisto("SingleValue/hZTauProd",     new TH1D("ZTauProd",    "#tau production point",     200, 0., 27.));
    BookHisto("SingleValue/hZDDecay",      new TH1D("ZDDecay",     "N production point",        200, 0., 27.));
    BookHisto("SingleValue/hDTheta",       new TH1D("DTheta",  "N mother polar angle",        100,  0., 0.3));
    BookHisto("SingleValue/hDLambda",      new TH1D("DLambda", "N mother decay length", 100, -1., 40.));
    BookHisto("SingleValue/hDPath",        new TH1D("DPath",   "N mother path along Z", 100, -1., 50.));
    BookHisto("SingleValue/hDMom",         new TH1D("DMom",    "N mother momentum",     100, -1., 170.));
    BookHisto("SingleValue/hZHNLDecay",    new TH1D("ZHNLDecay", "N decay point", 100., 90., 190.));
    BookHisto("SingleValue/hHNLGamma",     new TH1D("HNLGamma", "N Lorentz gamma", 50., 0., 170.));
    BookHisto("SingleValue/hHNLDecayProb", new TH1D("HNLDecayProb", "N decay probability", 100., 0., 0.006));
    BookHisto("SingleValue/hHNLReachProb", new TH1D("HNLReachProb", "N reach probability", 100., 0.99, 1.));
    BookHisto("SingleValue/hHNLTheta",     new TH1D("HNLTheta", "N polar angle", 100., 0., 0.5));
    BookHisto("SingleValue/hHNLMom",       new TH1D("HNLMom", "N momentum", 100., -0.5, 200.));
    BookHisto("SingleValue/hWeight",       new TH1D("Weight", "N weight", 1000, 1.E-20, 1.E-7));
    BookHisto("SingleValue/hCoupling",     new TH1D("Coupling", "Coupling", 100, -10., 0.));
    BookHisto("SingleValue/hMass",         new TH1D("Mass", "Mass", 10, 0.3, 1.2));
    BookHisto("SingleValue/hAcc",          new TH1D("Acc",       "Acceptance for one value of N mass and coupling",                             1, 1.E-30, 1.));
    BookHisto("SingleValue/hAccTarget",    new TH1D("AccTarget", "Acceptance for one value of N mass and coupling, for target-produced events", 1, 1.E-30, 1.));
    BookHisto("SingleValue/hAccTAX",       new TH1D("AccTAX",    "Acceptance for one value of N mass and coupling, for TAX-produced events",    1, 1.E-30, 1.));
    BookHisto("SingleValue/hYield",        new TH1D("Yield",       "Yield per POT for one value of N mass and coupling",                                1, 1.E-50, 1.E-10));
    BookHisto("SingleValue/hYieldTarget",  new TH1D("YieldTarget", "Yield per POT for one value of N mass and coupling, for target-produced events",    1, 1.E-50, 1.E-10));
    BookHisto("SingleValue/hYieldTAX",     new TH1D("YieldTAX",    "Yield per POT for one value of N mass and coupling, for TAX-produced events",       1, 1.E-50, 1.E-10));
    BookHisto("SingleValue/hTransverse",   new TH1D("Transverse",  "N transverse position at decay point",                                              100, 0., 2.));
    BookHisto("SingleValue/hDThetaMom",    new TH2D("DThetaMom",   "D meson polar angle vs momentum", 100, 0., 150., 50, 0., 0.5));
    BookHisto("SingleValue/hXYDecay",      new TH2D("XYDecay",     "X, Y of N at decay point", 100, 0., 1., 100, 0., 1.));

    fgAccMom = new TGraphErrors();
    fgAccMom->SetNameTitle("SingleValue/AccMom", "Acceptance vs N momentum");
    BookHisto(fgAccMom);
    
    fgYieldMom = new TGraphErrors();
    fgYieldMom->SetNameTitle("SingleValue/YieldMom", "Yield per POT vs N momentum");
    BookHisto(fgYieldMom);

    // Coupling scan 
    
    BookHisto("CouplingScan/hReachCoupling",  new TH2D("ReachCoupling",  "Probability of N reaching the FV vs coupling",                 fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
    BookHisto("CouplingScan/hDecayCoupling",  new TH2D("DecayCoupling",  "Probability of N decaying in the FV vs coupling",              fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
    BookHisto("CouplingScan/hProbCoupling",   new TH2D("ProbCoupling",   "Probability of N reaching and decaying in the FV vs coupling", fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
    BookHisto("CouplingScan/hWeightCoupling", new TH2D("WeightCoupling", "N weight vs coupling",                                         fN, fCouplingStart, fCouplingStop, 1000, 1.E-13, 1.E-6));
    BookHisto("CouplingScan/hEnergyCoupling", new TH2D("EnergyCoupling", "N energy vs coupling",                                         fN, fCouplingStart, fCouplingStop, 100, 0., 100.));
    
    fgGammaTotCoupling = new TGraphErrors();
    fgGammaTotCoupling->SetNameTitle("CouplingScan/GammaTotCoupling", "N total decay width vs coupling");
    BookHisto(fgGammaTotCoupling);
    
    fgTauCoupling = new TGraphErrors();
    fgTauCoupling->SetNameTitle("CouplingScan/TauCoupling", "N lifetime vs coupling");
    BookHisto(fgTauCoupling);
    
    fgAccCoupling = new TGraphErrors();
    fgAccCoupling->SetNameTitle("CouplingScan/AccCoupling", "All events");
    BookHisto(fgAccCoupling);
    
    fgAccCouplingTarget = new TGraphErrors();
    fgAccCouplingTarget->SetNameTitle("CouplingScan/AccCouplingTarget", "Target-produced events");
    BookHisto(fgAccCouplingTarget);
    
    fgAccCouplingTAX = new TGraphErrors();
    fgAccCouplingTAX->SetNameTitle("CouplingScan/AccCouplingTAX", "TAX-produced events");
    BookHisto(fgAccCouplingTAX);
    
    fgYieldCoupling = new TGraphErrors();
    fgYieldCoupling->SetNameTitle("CouplingScan/YieldCoupling", "All events");
    BookHisto(fgYieldCoupling);
    
    fgYieldCouplingTarget = new TGraphErrors();
    fgYieldCouplingTarget->SetNameTitle("CouplingScan/YieldCouplingTarget", "Target-produced events");
    BookHisto(fgYieldCouplingTarget);
    
    fgYieldCouplingTAX = new TGraphErrors();
    fgYieldCouplingTAX->SetNameTitle("CouplingScan/YieldCouplingTAX", "TAX-produced events");
    BookHisto(fgYieldCouplingTAX);
    
    // Mass scan
    
    BookHisto("MassScan/hReachMass",  new TH2D("ReachMass",  "Probability of N reaching the FV vs N mass",                 10, 0.25, 1.25, 1000, -0.1, 1.1));
    BookHisto("MassScan/hDecayMass",  new TH2D("DecayMass",  "Probability of N decaying in the FV vs N mass",              10, 0.25, 1.25, 1000, -0.1, 1.1));
    BookHisto("MassScan/hProbMass",   new TH2D("ProbMass",   "Probability of N reaching and decaying in the FV vs N mass", 10, 0.25, 1.25, 1000, -0.1, 1.1));
    BookHisto("MassScan/hWeightMass", new TH2D("WeightMass", "N weight vs N mass",                                         10, 0.25, 1.25, 1000, 1.E-14, 1.E-10));
    BookHisto("MassScan/hEnergyMass", new TH2D("EnergyMass", "N energy vs N mass",                                         10, 0.25, 1.25, 100, 0., 100.));

    fgGammaTotMass = new TGraphErrors();
    fgGammaTotMass->SetNameTitle("MassScan/GammaTotMass", "N total decay width vs N mass");
    BookHisto(fgGammaTotMass);
    
    fgTauMass = new TGraphErrors();
    fgTauMass->SetNameTitle("MassScan/TauMass", "N lifetime vs N mass");
    BookHisto(fgTauMass);
    
    fgAccMass = new TGraphErrors();
    fgAccMass->SetNameTitle("MassScan/AccMass", "All events");
    BookHisto(fgAccMass);
    
    fgAccMassTarget = new TGraphErrors();
    fgAccMassTarget->SetNameTitle("MassScan/AccMassTarget", "Target-produced events");
    BookHisto(fgAccMassTarget);
    
    fgAccMassTAX = new TGraphErrors();
    fgAccMassTAX->SetNameTitle("MassScan/AccMassTAX", "TAX-produced events");
    BookHisto(fgAccMassTAX);
    
    fgYieldMass = new TGraphErrors();
    fgYieldMass->SetNameTitle("MassScan/YieldMass", "All events");
    BookHisto(fgYieldMass);
    
    fgYieldMassTarget = new TGraphErrors();
    fgYieldMassTarget->SetNameTitle("MassScan/YieldMassTarget", "Target-produced events");
    BookHisto(fgYieldMassTarget);
    
    fgYieldMassTAX = new TGraphErrors();
    fgYieldMassTAX->SetNameTitle("MassScan/YieldMassTAX", "TAX-produced events");
    BookHisto(fgYieldMassTAX);
    
    // Total scan
    
    fgExclusion = new TGraphErrors();
    fgExclusion->SetNameTitle("TotalScan/Exclusion", "Sensitivity as a function of coupling and N mass");
    BookHisto(fgExclusion);

    // Toy-MC comparison

    BookHisto("ToyMC/DS/hpDS",   new TH1D("pDS",   "D_{S} momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs",                  40, 0., 200.));
    BookHisto("ToyMC/DS/hpDS1",  new TH1D("pDS1",  "D_{S} momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only",            20, 40., 200.));
    BookHisto("ToyMC/DS/hptDS",  new TH1D("ptDS",  "D_{S} transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs",       10, 0., 5.));
    BookHisto("ToyMC/DS/hptDS1", new TH1D("ptDS1", "D_{S} transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only",  6, 0., 3.));

    BookHisto("ToyMC/DS/hpNS",   new TH1D("pNS",   "N momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs",                  50, 0., 300.));
    BookHisto("ToyMC/DS/hpNS1",  new TH1D("pNS1",  "N momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only",            20, 40., 300.));
    BookHisto("ToyMC/DS/hptNS",  new TH1D("ptNS",  "N transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs",       20, 0., 2.));
    BookHisto("ToyMC/DS/hptNS1", new TH1D("ptNS1", "N transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 10, 0., 2.));

    BookHisto("ToyMC/DS/hpmuS",   new TH1D("pmuS",   "Muon daughter momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs",                  50, 0., 300.));
    BookHisto("ToyMC/DS/hpmuS1",  new TH1D("pmuS1",  "Muon daughter momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only",            20, 0., 120.));
    BookHisto("ToyMC/DS/hptmuS",  new TH1D("ptmuS",  "Muon daughter transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs",       20, 0., 2.));
    BookHisto("ToyMC/DS/hptmuS1", new TH1D("ptmuS1", "Muon daughter transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 10, 0., 2.));

    BookHisto("ToyMC/D0/hpD0",   new TH1D("pD0",   "D^{0} momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs",                  40, 0., 200.));
    BookHisto("ToyMC/D0/hpD01",  new TH1D("pD01",  "D^{0} momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only",            40, 60., 300.));
    BookHisto("ToyMC/D0/hptD0",  new TH1D("ptD0",  "D^{0} transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs",       10, 0., 5.));
    BookHisto("ToyMC/D0/hptD01", new TH1D("ptD01", "D^{0} transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only",  6, 0., 3.));

    BookHisto("ToyMC/D0/hpN0",   new TH1D("pN0",   "N momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs",                  50, 0., 300.));
    BookHisto("ToyMC/D0/hpN01",  new TH1D("pN01",  "N momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only",            20, 40., 300.));
    BookHisto("ToyMC/D0/hptN0",  new TH1D("ptN0",  "N transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs",       20, 0., 2.));
    BookHisto("ToyMC/D0/hptN01", new TH1D("ptN01", "N transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 10, 0., 2.));

    BookHisto("ToyMC/D0/hpmu0",   new TH1D("pmu0",   "Muon daughter momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs",                  50, 0., 300.));
    BookHisto("ToyMC/D0/hpmu01",  new TH1D("pmu01",  "Muon daughter momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only",            20, 0., 120.));
    BookHisto("ToyMC/D0/hptmu0",  new TH1D("ptmu0",  "Muon daughter transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs",       20, 0., 2.));
    BookHisto("ToyMC/D0/hptmu01", new TH1D("ptmu01", "Muon daughter transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 10, 0., 2.));
  }
  else {

    // Errorgraphs for mass and coupling scans

    fgErrorAccCoupling         = (TGraphErrors*)RequestHistogram(fAnalyzerName, "CouplingScan/AccCoupling",         true);
    fgErrorAccCouplingTarget   = (TGraphErrors*)RequestHistogram(fAnalyzerName, "CouplingScan/AccCouplingTarget",   true);
    fgErrorAccCouplingTAX      = (TGraphErrors*)RequestHistogram(fAnalyzerName, "CouplingScan/AccCouplingTAX",      true);
    fgErrorYieldCoupling       = (TGraphErrors*)RequestHistogram(fAnalyzerName, "CouplingScan/YieldCoupling",       true);
    fgErrorYieldCouplingTarget = (TGraphErrors*)RequestHistogram(fAnalyzerName, "CouplingScan/YieldCouplingTarget", true);
    fgErrorYieldCouplingTAX    = (TGraphErrors*)RequestHistogram(fAnalyzerName, "CouplingScan/YieldCouplingTAX",    true);
    fgErrorAccMass             = (TGraphErrors*)RequestHistogram(fAnalyzerName, "MassScan/AccMass",                 true);
    fgErrorAccMassTarget       = (TGraphErrors*)RequestHistogram(fAnalyzerName, "MassScan/AccMassTarget",           true);
    fgErrorAccMassTAX          = (TGraphErrors*)RequestHistogram(fAnalyzerName, "MassScan/AccMassTAX",              true);
    fgErrorYieldMass           = (TGraphErrors*)RequestHistogram(fAnalyzerName, "MassScan/YieldMass",               true);
    fgErrorYieldMassTarget     = (TGraphErrors*)RequestHistogram(fAnalyzerName, "MassScan/YieldMassTarget",         true);
    fgErrorYieldMassTAX        = (TGraphErrors*)RequestHistogram(fAnalyzerName, "MassScan/YieldMassTAX",            true);
    fgErrorAccMom              = (TGraphErrors*)RequestHistogram(fAnalyzerName, "SingleValue/AccMom",               true);
    fgErrorYieldMom            = (TGraphErrors*)RequestHistogram(fAnalyzerName, "SingleValue/YieldMom",             true);
    fgErrorExclusion           = (TGraphErrors*)RequestHistogram(fAnalyzerName, "TotalScan/Exclusion",              true);
  }
}

void HeavyNeutrinoScan::Process(Int_t) {

  if (!fReadingData) return;

  Double_t fCoupling          = -999.;  
  Double_t MN                 = 0.;
  Double_t HNLTau             = 0.;
  Double_t gammaTot           = 0.;
  Double_t momN               = 0.;
  Double_t momBin             = 0.;
  Double_t NDecayProb         = 0.;
  Double_t NReachProb         = 0.;
  Double_t Weight             = 0.;
  Bool_t isGood               = false;

  Bool_t IsHNLGood = *(Bool_t*)GetOutput("HeavyNeutrino.Output");

  // Scan on the coupling                                                                       
  
  for(Int_t couplingIndex  = fCouplingStart*10; couplingIndex <= fCouplingStop*10; couplingIndex += fCouplingStep*10) {
    fCoupling = couplingIndex/10.;
    fUSquared = TMath::Power(10., fCoupling);
    fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
    fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
    fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;
    fCouplings[fCoupling] = fCoupling;

    if (GetWithMC()) {
      Event *evt = GetMCEvent();

      std::vector<std::map<std::string, Double_t>> Weights = ComputeWeight(evt, fUSquared, fUeSquaredRatio, fUmuSquaredRatio, fUtauSquaredRatio, fLInitialFV, fLFV, fMode);

      for (UInt_t i = 0; i < Weights.size(); i++) {
	MN =  Weights[i]["Mass"];
	HNLTau = Weights[i]["Lifetime"];
	momN = Weights[i]["Momentum"];
	gammaTot = GammaTot(MN);
	NReachProb = Weights[i]["ReachProb"];
	NDecayProb = Weights[i]["DecayProb"];
	Weight = Weights[i]["Weight"];    
	isGood = Weights[i]["IsGood"];
	fMasses[round(MN)] = round(MN);
	fGammaTot[round(MN)][fCoupling] = gammaTot;
	fTau[round(MN)][fCoupling] = HNLTau;

	if (fNevents[round(MN)].count(fCoupling) == 0)
	  fNevents[round(MN)][fCoupling] = 0;
	fNevents[round(MN)][fCoupling]++;
	fSumAll[round(MN)][fCoupling] += Weight;

	if (Weights[i]["ProdProb"] == fDBeProdProb) {
	  if (fNeventsTarget[round(MN)].count(fCoupling) == 0)
	    fNeventsTarget[round(MN)][fCoupling] = 0;
	  fNeventsTarget[round(MN)][fCoupling]++;
	  fSumAllTarget[round(MN)][fCoupling] += Weight;
	}
	else if (Weights[i]["ProdProb"] == fDCuProdProb) {
	  if (fNeventsTAX[round(MN)].count(fCoupling) == 0)
	    fNeventsTAX[round(MN)][fCoupling] = 0;
	  fNeventsTAX[round(MN)][fCoupling]++;
	  fSumAllTAX[round(MN)][fCoupling] += Weight;
	}

	if (IsHNLGood == true && isGood == true) {
	  fSumGood[round(MN)][fCoupling] += Weight;
	  if (Weights[i]["ProdProb"] == fDBeProdProb) { 
	    fSumGoodTarget[round(MN)][fCoupling] += Weight;
	  }
	  else if (Weights[i]["ProdProb"] == fDCuProdProb) {
	    fSumGoodTAX[round(MN)][fCoupling] += Weight;
	  }
	}
	
	if (i == 0) {
	  if (fEvtCounter[round(MN)].count(fCoupling) == 0)
	    fEvtCounter[round(MN)][fCoupling] = 0;
	  fEvtCounter[round(MN)][fCoupling]++;
	}

	if (fEvtCounter[round(MN)][fCoupling] % fSplitStep == 0 && i == Weights.size() - 1) { 
	  if (fSumAll[round(MN)][fCoupling] != 0. && fNevents[round(MN)][fCoupling] != 0)
	    fErrorFile << round(MN) << "\t" << fCoupling << "\t" << fSumGood[round(MN)][fCoupling]/fSumAll[round(MN)][fCoupling] << "\t" << fSumGood[round(MN)][fCoupling]/fNevents[round(MN)][fCoupling] << endl;
	  else
	    fErrorFile << round(MN) << "\t" << fCoupling << "\t" << "0." << "\t" << "0." << endl;
	  if (fSumAllTarget[round(MN)][fCoupling] != 0. && fNeventsTarget[round(MN)][fCoupling] != 0)
	    fErrorFileTarget << round(MN) << "\t" << fCoupling << "\t" << fSumGoodTarget[round(MN)][fCoupling]/fSumAllTarget[round(MN)][fCoupling] << "\t" << fSumGoodTarget[round(MN)][fCoupling]/fNeventsTarget[round(MN)][fCoupling] << endl;
	  else
	    fErrorFileTarget << round(MN) << "\t" << fCoupling << "\t" << "0." << "\t" << "0." << endl;
	  if (fSumAllTAX[round(MN)][fCoupling] != 0. && fNeventsTAX[round(MN)][fCoupling] != 0)
	    fErrorFileTAX << round(MN) << "\t" << fCoupling << "\t" << fSumGoodTAX[round(MN)][fCoupling]/fSumAllTAX[round(MN)][fCoupling] << "\t" << fSumGoodTAX[round(MN)][fCoupling]/fNeventsTAX[round(MN)][fCoupling] << endl;
	  else
	    fErrorFileTAX << round(MN) << "\t" << fCoupling << "\t" << "0." << "\t" << "0." << endl;
	}
	
	if (round(MN) == fMassForSingleValue) {
	  FillHisto("CouplingScan/hReachCoupling",  fCoupling, NReachProb);
	  FillHisto("CouplingScan/hDecayCoupling",  fCoupling, NDecayProb);
	  FillHisto("CouplingScan/hProbCoupling",   fCoupling, NReachProb*NDecayProb);
	  FillHisto("CouplingScan/hWeightCoupling", fCoupling, Weight);	
	  FillHisto("CouplingScan/hEnergyCoupling", fCoupling, momN/1000.*momN/1000. + MN/1000.*MN/1000.);	
	}

	if (fCoupling == fCouplingForSingleValue) {
	  FillHisto("MassScan/hReachMass",  round(MN)/1000., NReachProb);
	  FillHisto("MassScan/hDecayMass",  round(MN)/1000., NDecayProb);
	  FillHisto("MassScan/hProbMass",   round(MN)/1000., NReachProb*NDecayProb);
	  FillHisto("MassScan/hWeightMass", round(MN)/1000., Weight);
	  FillHisto("MassScan/hEnergyMass", round(MN)/1000., momN/1000.*momN/1000. + MN/1000.*MN/1000.);	
	}

	if (fCoupling == fCouplingForSingleValue && round(MN) == fMassForSingleValue) {
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
    std::vector<std::map<std::string, Double_t>> Weights = ComputeWeight(evt, fUSquared, fUeSquaredRatio, fUmuSquaredRatio, fUtauSquaredRatio, fLInitialFV, fLFV, fMode);

    for (UInt_t i = 0; i < Weights.size(); i++) {
      MN =  round(Weights[i]["Mass"]);
      isGood = Weights[i]["IsGood"];
      if (MN == fMassForSingleValue) {
	momN = Weights[i]["Momentum"]/1000.;
	Weight = Weights[i]["Weight"];

	if (momN >= fMomStart && momN <= fMomStop) {
	  momBin = fMomStep*trunc(momN/fMomStep);
	  fMomenta[momBin] = momBin;

	  if (fNeventsMom.count(momBin) == 0)
	    fNeventsMom[momBin] = 0;
	  fNeventsMom[momBin]++;
	  fSumAllMom[momBin] += Weight;
	  
	  if (IsHNLGood == true && isGood == true) {
	    fSumGoodMom[momBin] += Weight;
	  }

	  if (i == 0) {
	    if (fEvtCounterMom.count(momBin) == 0)
	      fEvtCounterMom[momBin] = 0;
	    fEvtCounterMom[momBin]++;
	  }

	  if (fEvtCounterMom[momBin] % (fSplitStep*2) == 0 && i == Weights.size() - 1) {
	    if (fSumAllMom[momBin] != 0. && fNeventsMom[momBin] != 0) {
	      fErrorFileMom << momBin << "\t" << fSumGoodMom[momBin]/fSumAllMom[momBin] << "\t" << fSumGoodMom[momBin]/fNeventsMom[momBin] << endl;
	    }
	    else
	      fErrorFileMom << momBin << "\t" << "0." << "\t" << "0." << endl;
	  }
	}
      }
    }
  }

  // Some plots
  
  if (GetWithMC()) {
    Event *evt = GetMCEvent();
    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);
      if (p->GetParentID() == -1 && p->GetPDGcode() == 999) {
	FillHisto("SingleValue/hZMotherProd", p->GetPosAtCheckPoint(0).z()/1000.);

	if (p->GetPosAtCheckPoint(1).z() != 0.)
	  FillHisto("SingleValue/hZDProd", p->GetPosAtCheckPoint(1).z()/1000.);
	if (p->GetPosAtCheckPoint(2).z() != 0.)
	  FillHisto("SingleValue/hZTauProd", p->GetPosAtCheckPoint(2).z()/1000.);

	FillHisto("SingleValue/hZDDecay",    p->GetProdPos().Z()/1000.);
	FillHisto("SingleValue/hDTheta",     p->GetPosAtCheckPoint(1).x());
	FillHisto("SingleValue/hDLambda",    p->GetPosAtCheckPoint(1).y());
	FillHisto("SingleValue/hDPath",      p->GetMomAtCheckPoint(1).X());
	FillHisto("SingleValue/hDMom",       p->GetMomAtCheckPoint(1).Y()/1000.);
	FillHisto("SingleValue/hZHNLDecay",  p->GetEndPos().Z()/1000.);
	FillHisto("SingleValue/hHNLGamma",   p->GetInitial4Momentum().Gamma());
	FillHisto("SingleValue/hHNLTheta",   p->GetMomAtCheckPoint(0).Z());
	FillHisto("SingleValue/hHNLMom",     p->GetMomAtCheckPoint(0).T()/1000.);
	FillHisto("SingleValue/hDThetaMom",  p->GetMomAtCheckPoint(1).Y()/1000., p->GetPosAtCheckPoint(1).x());
	FillHisto("SingleValue/hXYDecay",    p->GetEndPos().X()/1000., p->GetEndPos().Y()/1000.);
	FillHisto("SingleValue/hTransverse", TMath::Sqrt(p->GetEndPos().X()/1000.*p->GetEndPos().X()/1000. + p->GetEndPos().Y()/1000.*p->GetEndPos().Y()/1000.));

	if (TMath::Sqrt(p->GetInitial4Momentum().E()*p->GetInitial4Momentum().E() - p->GetInitial4Momentum().P()*p->GetInitial4Momentum().P()) == fMassForSingleValue) {
	  if (p->GetMomAtCheckPoint(3).X() != 0.) {
	    FillHisto("ToyMC/DS/hpDS",   p->GetMomAtCheckPoint(3).X()/1000.);
	    FillHisto("ToyMC/DS/hptDS",  p->GetMomAtCheckPoint(3).Y()/1000.);
	    FillHisto("ToyMC/DS/hpNS",   p->GetMomAtCheckPoint(3).Z()/1000.);
	    FillHisto("ToyMC/DS/hptNS",  p->GetMomAtCheckPoint(3).T()/1000.);
	    FillHisto("ToyMC/DS/hpmuS",  p->GetMomAtCheckPoint(4).X()/1000.);
	    FillHisto("ToyMC/DS/hptmuS", p->GetMomAtCheckPoint(4).Y()/1000.);
	  }
	  
	  if (p->GetMomAtCheckPoint(5).X() != 0.) {
	    FillHisto("ToyMC/D0/hpD0",   p->GetMomAtCheckPoint(5).X()/1000.);
	    FillHisto("ToyMC/D0/hptD0",  p->GetMomAtCheckPoint(5).Y()/1000.);
	    FillHisto("ToyMC/D0/hpN0",   p->GetMomAtCheckPoint(5).Z()/1000.);
	    FillHisto("ToyMC/D0/hptN0",  p->GetMomAtCheckPoint(5).T()/1000.);
	    FillHisto("ToyMC/D0/hpmu0",  p->GetMomAtCheckPoint(6).X()/1000.);
	    FillHisto("ToyMC/D0/hptmu0", p->GetMomAtCheckPoint(6).Y()/1000.);
	  }
	  
	  if (p->GetEndProcessName() == "good") {
	    if (p->GetMomAtCheckPoint(3).X() != 0.) {
	      FillHisto("ToyMC/DS/hpDS1",   p->GetMomAtCheckPoint(3).X()/1000.);
	      FillHisto("ToyMC/DS/hptDS1",  p->GetMomAtCheckPoint(3).Y()/1000.);
	      FillHisto("ToyMC/DS/hpNS1",   p->GetMomAtCheckPoint(3).Z()/1000.);
	      FillHisto("ToyMC/DS/hptNS1",  p->GetMomAtCheckPoint(3).T()/1000.);
	      FillHisto("ToyMC/DS/hpmuS1",  p->GetMomAtCheckPoint(4).X()/1000.);
	      FillHisto("ToyMC/DS/hptmuS1", p->GetMomAtCheckPoint(4).Y()/1000.);
	    }
	    
	    if (p->GetMomAtCheckPoint(5).X() != 0.) {
	      FillHisto("ToyMC/D0/hpD01",   p->GetMomAtCheckPoint(5).X()/1000.);
	      FillHisto("ToyMC/D0/hptD01",  p->GetMomAtCheckPoint(5).Y()/1000.);
	      FillHisto("ToyMC/D0/hpN01",   p->GetMomAtCheckPoint(5).Z()/1000.);
	      FillHisto("ToyMC/D0/hptN01",  p->GetMomAtCheckPoint(5).T()/1000.);
	      FillHisto("ToyMC/D0/hpmu01",  p->GetMomAtCheckPoint(6).X()/1000.);
	      FillHisto("ToyMC/D0/hptmu01", p->GetMomAtCheckPoint(6).Y()/1000.);
	    }
	  }
	}
      }
    }
  }
}

void HeavyNeutrinoScan::EndOfJobUser() {

  if (fReadingData) {
  
    // Retrieve histos

    fhZMotherProd  = (TH1D*)fHisto.GetTH1("SingleValue/hZMotherProd");
    fhZDProd       = (TH1D*)fHisto.GetTH1("SingleValue/hZDProd");
    fhZTauProd     = (TH1D*)fHisto.GetTH1("SingleValue/hZTauProd");
    fhZDDecay      = (TH1D*)fHisto.GetTH1("SingleValue/hZDDecay");
    fhDTheta       = (TH1D*)fHisto.GetTH1("SingleValue/hDTheta");
    fhDLambda      = (TH1D*)fHisto.GetTH1("SingleValue/hDLambda");
    fhDPath        = (TH1D*)fHisto.GetTH1("SingleValue/hDPath");
    fhDMom         = (TH1D*)fHisto.GetTH1("SingleValue/hDMom");
    fhZHNLDecay    = (TH1D*)fHisto.GetTH1("SingleValue/hZHNLDecay");
    fhHNLGamma     = (TH1D*)fHisto.GetTH1("SingleValue/hHNLGamma");
    fhHNLDecayProb = (TH1D*)fHisto.GetTH1("SingleValue/hHNLDecayProb");
    fhHNLReachProb = (TH1D*)fHisto.GetTH1("SingleValue/hHNLReachProb");
    fhHNLTheta     = (TH1D*)fHisto.GetTH1("SingleValue/hHNLTheta");
    fhHNLMom       = (TH1D*)fHisto.GetTH1("SingleValue/hHNLMom");
    fhWeight       = (TH1D*)fHisto.GetTH1("SingleValue/hWeight");
    fhCoupling     = (TH1D*)fHisto.GetTH1("SingleValue/hCoupling");
    fhMass         = (TH1D*)fHisto.GetTH1("SingleValue/hMass");
    fhAcc          = (TH1D*)fHisto.GetTH1("SingleValue/hAcc");
    fhAccTarget    = (TH1D*)fHisto.GetTH1("SingleValue/hAccTarget");
    fhAccTAX       = (TH1D*)fHisto.GetTH1("SingleValue/hAccTAX");
    fhYield        = (TH1D*)fHisto.GetTH1("SingleValue/hYield");
    fhYieldTarget  = (TH1D*)fHisto.GetTH1("SingleValue/hYieldTarget");
    fhYieldTAX     = (TH1D*)fHisto.GetTH1("SingleValue/hYieldTAX");
    fhTransverse   = (TH1D*)fHisto.GetTH1("SingleValue/hTransverse");
    fhDThetaMom    = (TH2D*)fHisto.GetTH2("SingleValue/hDThetaMom");
    fhXYDecay      = (TH2D*)fHisto.GetTH2("SingleValue/hXYDecay");

    fhReachCoupling  = (TH2D*)fHisto.GetTH2("CouplingScan/hReachCoupling");
    fhDecayCoupling  = (TH2D*)fHisto.GetTH2("CouplingScan/hDecayCoupling");
    fhProbCoupling   = (TH2D*)fHisto.GetTH2("CouplingScan/hProbCoupling");
    fhWeightCoupling = (TH2D*)fHisto.GetTH2("CouplingScan/hWeightCoupling");
    fhEnergyCoupling = (TH2D*)fHisto.GetTH2("CouplingScan/hEnergyCoupling");

    fhReachMass  = (TH2D*)fHisto.GetTH2("MassScan/hReachMass");
    fhDecayMass  = (TH2D*)fHisto.GetTH2("MassScan/hDecayMass");
    fhProbMass   = (TH2D*)fHisto.GetTH2("MassScan/hProbMass");
    fhWeightMass = (TH2D*)fHisto.GetTH2("MassScan/hWeightMass");
    fhEnergyMass = (TH2D*)fHisto.GetTH2("MassScan/hEnergyMass");

    fhpDS    = (TH1D*)fHisto.GetTH1("ToyMC/DS/hpDS");
    fhpDS1   = (TH1D*)fHisto.GetTH1("ToyMC/DS/hpDS1");
    fhptDS   = (TH1D*)fHisto.GetTH1("ToyMC/DS/hptDS");
    fhptDS1  = (TH1D*)fHisto.GetTH1("ToyMC/DS/hptDS1");
    fhpNS    = (TH1D*)fHisto.GetTH1("ToyMC/DS/hpNS");
    fhpNS1   = (TH1D*)fHisto.GetTH1("ToyMC/DS/hpNS1");
    fhptNS   = (TH1D*)fHisto.GetTH1("ToyMC/DS/hptNS");
    fhptNS1  = (TH1D*)fHisto.GetTH1("ToyMC/DS/hptNS1");
    fhpmuS   = (TH1D*)fHisto.GetTH1("ToyMC/DS/hpmuS");
    fhpmuS1  = (TH1D*)fHisto.GetTH1("ToyMC/DS/hpmuS1");
    fhptmuS  = (TH1D*)fHisto.GetTH1("ToyMC/DS/hptmuS");
    fhptmuS1 = (TH1D*)fHisto.GetTH1("ToyMC/DS/hptmuS1");
    fhpD0    = (TH1D*)fHisto.GetTH1("ToyMC/D0/hpD0");
    fhpD01   = (TH1D*)fHisto.GetTH1("ToyMC/D0/hpD01");
    fhptD0   = (TH1D*)fHisto.GetTH1("ToyMC/D0/hptD0");
    fhptD01  = (TH1D*)fHisto.GetTH1("ToyMC/D0/hptD01");
    fhpN0    = (TH1D*)fHisto.GetTH1("ToyMC/D0/hpN0");
    fhpN01   = (TH1D*)fHisto.GetTH1("ToyMC/D0/hpN01");
    fhptN0   = (TH1D*)fHisto.GetTH1("ToyMC/D0/hptN0");
    fhptN01  = (TH1D*)fHisto.GetTH1("ToyMC/D0/hptN01");
    fhpmu0   = (TH1D*)fHisto.GetTH1("ToyMC/D0/hpmu0");
    fhpmu01  = (TH1D*)fHisto.GetTH1("ToyMC/D0/hpmu01");
    fhptmu0  = (TH1D*)fHisto.GetTH1("ToyMC/D0/hptmu0");
    fhptmu01 = (TH1D*)fHisto.GetTH1("ToyMC/D0/hptmu01");
    
    // X axis title

    fhZMotherProd ->GetXaxis()->SetTitle("Position along Z [m]");
    fhZDProd      ->GetXaxis()->SetTitle("Position along Z [m]");
    fhZTauProd    ->GetXaxis()->SetTitle("Position along Z [m]");
    fhZDDecay     ->GetXaxis()->SetTitle("Position along Z [m]");
    fhDTheta      ->GetXaxis()->SetTitle("Theta [rad]");
    fhDLambda     ->GetXaxis()->SetTitle("Decay length [mm]");
    fhDPath       ->GetXaxis()->SetTitle("Z [mm]");
    fhDMom        ->GetXaxis()->SetTitle("P [GeV/c]");
    fhZHNLDecay   ->GetXaxis()->SetTitle("Position along Z [m]");
    fhHNLGamma    ->GetXaxis()->SetTitle("Lorentz gamma");
    fhHNLDecayProb->GetXaxis()->SetTitle("Decay probability");
    fhHNLReachProb->GetXaxis()->SetTitle("Reach probability");
    fhHNLTheta    ->GetXaxis()->SetTitle("Theta [rad]");
    fhHNLMom      ->GetXaxis()->SetTitle("P [GeV/c]");
    fhWeight      ->GetXaxis()->SetTitle("Weight");
    fhCoupling    ->GetXaxis()->SetTitle("Coupling");
    fhMass        ->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fhAcc         ->GetXaxis()->SetTitle("Acceptance");
    fhAccTarget   ->GetXaxis()->SetTitle("Acceptance");
    fhAccTAX      ->GetXaxis()->SetTitle("Acceptance");
    fhYield       ->GetXaxis()->SetTitle("Yield per POT");
    fhYieldTarget ->GetXaxis()->SetTitle("Yield per POT");
    fhYieldTAX    ->GetXaxis()->SetTitle("Yield per POT");
    fhTransverse  ->GetXaxis()->SetTitle("Distance [m]");
    fhDThetaMom   ->GetXaxis()->SetTitle("P [GeV/c]");
    fhXYDecay     ->GetXaxis()->SetTitle("X [m]");

    fhReachCoupling ->GetXaxis()->SetTitle("Log(U^{2})");
    fhDecayCoupling ->GetXaxis()->SetTitle("Log(U^{2})");
    fhProbCoupling  ->GetXaxis()->SetTitle("Log(U^{2})");
    fhWeightCoupling->GetXaxis()->SetTitle("Log(U^{2})");
    fhEnergyCoupling->GetXaxis()->SetTitle("Log(U^{2})");

    fhReachMass ->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fhDecayMass ->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fhProbMass  ->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fhWeightMass->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fhEnergyMass->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");

    fhpDS   ->GetXaxis()->SetTitle("P [GeV/c]");
    fhpDS1  ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptDS  ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptDS1 ->GetXaxis()->SetTitle("P [GeV/c]");
    fhpNS   ->GetXaxis()->SetTitle("P [GeV/c]");
    fhpNS1  ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptNS  ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptNS1 ->GetXaxis()->SetTitle("P [GeV/c]");
    fhpmuS  ->GetXaxis()->SetTitle("P [GeV/c]");
    fhpmuS1 ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptmuS ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptmuS1->GetXaxis()->SetTitle("P [GeV/c]");
    fhpD0   ->GetXaxis()->SetTitle("P [GeV/c]");
    fhpD01  ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptD0  ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptD01 ->GetXaxis()->SetTitle("P [GeV/c]");
    fhpN0   ->GetXaxis()->SetTitle("P [GeV/c]");
    fhpN01  ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptN0  ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptN01 ->GetXaxis()->SetTitle("P [GeV/c]");
    fhpmu0  ->GetXaxis()->SetTitle("P [GeV/c]");
    fhpmu01 ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptmu0 ->GetXaxis()->SetTitle("P [GeV/c]");
    fhptmu01->GetXaxis()->SetTitle("P [GeV/c]");   

    // Y axis title

    fhDThetaMom     ->GetYaxis()->SetTitle("Polar angle [rad]");
    fhXYDecay       ->GetYaxis()->SetTitle("Y [m]");
    fhReachCoupling ->GetYaxis()->SetTitle("Reach probability");
    fhDecayCoupling ->GetYaxis()->SetTitle("Decay probability");
    fhProbCoupling  ->GetYaxis()->SetTitle("Reach and decay probability");
    fhWeightCoupling->GetYaxis()->SetTitle("Weight");
    fhEnergyCoupling->GetYaxis()->SetTitle("N energy [GeV]");
    fhReachMass     ->GetYaxis()->SetTitle("Reach probability");
    fhDecayMass     ->GetYaxis()->SetTitle("Decay probability");
    fhProbMass      ->GetYaxis()->SetTitle("Reach and decay probability");
    fhWeightMass    ->GetYaxis()->SetTitle("Weight");
    fhEnergyMass    ->GetYaxis()->SetTitle("N energy [GeV]");

    fhpDS   ->GetYaxis()->SetTitle("Entries/5 GeV/c");
    fhpDS1  ->GetYaxis()->SetTitle("Entries/8 GeV/c");
    fhptDS  ->GetYaxis()->SetTitle("Entries/0.5 GeV/c");
    fhptDS1 ->GetYaxis()->SetTitle("Entries/0.5 GeV/c");
    fhpNS   ->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fhpNS1  ->GetYaxis()->SetTitle("Entries/13 GeV/c");
    fhptNS  ->GetYaxis()->SetTitle("Entries/0.1 GeV/c");
    fhptNS1 ->GetYaxis()->SetTitle("Entries/0.2 GeV/c");
    fhpmuS  ->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fhpmuS1 ->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fhptmuS ->GetYaxis()->SetTitle("Entries/0.1 GeV/c");
    fhptmuS1->GetYaxis()->SetTitle("Entries/0.2 GeV/c");
    fhpD0   ->GetYaxis()->SetTitle("Entries/5 GeV/c");
    fhpD01  ->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fhptD0  ->GetYaxis()->SetTitle("Entries/0.5 GeV/c");
    fhptD01 ->GetYaxis()->SetTitle("Entries/0.5 GeV/c");
    fhpN0   ->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fhpN01  ->GetYaxis()->SetTitle("Entries/13 GeV/c");
    fhptN0  ->GetYaxis()->SetTitle("Entries/0.1 GeV/c");
    fhptN01 ->GetYaxis()->SetTitle("Entries/0.2 GeV/c");
    fhpmu0  ->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fhpmu01 ->GetYaxis()->SetTitle("Entries/6 GeV/c");
    fhptmu0 ->GetYaxis()->SetTitle("Entries/0.1 GeV/c");
    fhptmu01->GetYaxis()->SetTitle("Entries/0.2 GeV/c");   

    // Acceptance computation                                                                       

    Double_t Coupling          = 0.;
    Double_t MN                = 0.;
    Double_t Mom               = 0.;
    Int_t counter              = 0;
    Int_t couplingCounter      = 0;
    Int_t massCounter          = 0;
    Int_t momCounter           = 0;

    for (auto it = fMasses.begin(); it != fMasses.end(); it++) {
      MN = it->first;
      for (auto it1 = fCouplings.begin(); it1 != fCouplings.end(); it1++) {
	Coupling = it1->first;

	fSumAll       [MN][Coupling] != 0 ? fAcc        [MN][Coupling] =       fSumGood[MN][Coupling]/       fSumAll[MN][Coupling] : fAcc        [MN][Coupling] = 0;
	fSumAllTarget [MN][Coupling] != 0 ? fAccTarget  [MN][Coupling] = fSumGoodTarget[MN][Coupling]/ fSumAllTarget[MN][Coupling] : fAccTarget  [MN][Coupling] = 0;
	fSumAllTAX    [MN][Coupling] != 0 ? fAccTAX     [MN][Coupling] =    fSumGoodTAX[MN][Coupling]/    fSumAllTAX[MN][Coupling] : fAccTAX     [MN][Coupling] = 0;
	fNevents      [MN][Coupling] != 0 ? fYield      [MN][Coupling] =       fSumGood[MN][Coupling]/      fNevents[MN][Coupling] : fYield      [MN][Coupling] = 0;
	fNeventsTarget[MN][Coupling] != 0 ? fYieldTarget[MN][Coupling] = fSumGoodTarget[MN][Coupling]/fNeventsTarget[MN][Coupling] : fYieldTarget[MN][Coupling] = 0;
	fNeventsTAX   [MN][Coupling] != 0 ? fYieldTAX   [MN][Coupling] =    fSumGoodTAX[MN][Coupling]/   fNeventsTAX[MN][Coupling] : fYieldTAX   [MN][Coupling] = 0;
	fNevents      [MN][Coupling] != 0 ? fProb       [MN][Coupling] =        fSumAll[MN][Coupling]/      fNevents[MN][Coupling] : fProb       [MN][Coupling] = 0;

	if (Coupling == fCouplingForSingleValue && MN == fMassForSingleValue) {
	  FillHisto("SingleValue/hAcc",                 fAcc[MN][Coupling]);
	  FillHisto("SingleValue/hAccTarget",     fAccTarget[MN][Coupling]);
	  FillHisto("SingleValue/hAccTAX",           fAccTAX[MN][Coupling]);
	  FillHisto("SingleValue/hYield",             fYield[MN][Coupling]);
	  FillHisto("SingleValue/hYieldTarget", fYieldTarget[MN][Coupling]);
	  FillHisto("SingleValue/hYieldTAX",       fYieldTAX[MN][Coupling]);
	}

	if (MN == fMassForSingleValue) {
	  fgGammaTotCoupling   ->SetPoint(couplingCounter, Coupling, fGammaTot   [MN][Coupling]);
	  fgTauCoupling        ->SetPoint(couplingCounter, Coupling, fTau        [MN][Coupling]);
	  fgAccCoupling        ->SetPoint(couplingCounter, Coupling, fAcc        [MN][Coupling]);
	  fgAccCouplingTarget  ->SetPoint(couplingCounter, Coupling, fAccTarget  [MN][Coupling]);
	  fgAccCouplingTAX     ->SetPoint(couplingCounter, Coupling, fAccTAX     [MN][Coupling]);
	  fgYieldCoupling      ->SetPoint(couplingCounter, Coupling, fYield      [MN][Coupling]);
	  fgYieldCouplingTarget->SetPoint(couplingCounter, Coupling, fYieldTarget[MN][Coupling]);
	  fgYieldCouplingTAX   ->SetPoint(couplingCounter, Coupling, fYieldTAX   [MN][Coupling]);
	  couplingCounter++;
	}
 
	if (fYield[MN][Coupling]*1.E18 > 2.3) {
	  fgExclusion->SetPoint(counter, MN/1000., Coupling);
	  counter++;
	}
      }
    
      Coupling = fCouplingForSingleValue;
      fgGammaTotMass       ->SetPoint(massCounter, MN/1000., fGammaTot   [MN][Coupling]);
      fgTauMass            ->SetPoint(massCounter, MN/1000., fTau        [MN][Coupling]);
      fgAccMass            ->SetPoint(massCounter, MN/1000., fAcc        [MN][Coupling]);
      fgAccMassTarget      ->SetPoint(massCounter, MN/1000., fAccTarget  [MN][Coupling]);
      fgAccMassTAX         ->SetPoint(massCounter, MN/1000., fAccTAX     [MN][Coupling]);
      fgYieldMass          ->SetPoint(massCounter, MN/1000., fYield      [MN][Coupling]);
      fgYieldMassTarget    ->SetPoint(massCounter, MN/1000., fYieldTarget[MN][Coupling]);
      fgYieldMassTAX       ->SetPoint(massCounter, MN/1000., fYieldTAX   [MN][Coupling]);
      massCounter++;
    }
  
    for (auto it = fMomenta.begin(); it != fMomenta.end(); it++) {
      Mom = it->first;
      fSumAllMom [Mom] != 0 ? fAccMom  [Mom] = fSumGoodMom[Mom]/ fSumAllMom[Mom] : fAccMom  [Mom] = 0.;
      fNeventsMom[Mom] != 0 ? fYieldMom[Mom] = fSumGoodMom[Mom]/fNeventsMom[Mom] : fYieldMom[Mom] = 0.;
      fgAccMom  ->SetPoint(momCounter, Mom, fAccMom  [Mom]);
      fgYieldMom->SetPoint(momCounter, Mom, fYieldMom[Mom]);
      momCounter++;
    }

    // Cosmetics

    CosmeticsGraph(fgGammaTotCoupling,    "Log(U^{2})",  "Total decay width [MeV]", 2);
    CosmeticsGraph(fgTauCoupling,         "Log(U^{2})",  "Lifetime [ns]",           2);
    CosmeticsGraph(fgAccCoupling,         "Log(U^{2})",  "Acceptance",              2);
    CosmeticsGraph(fgAccCouplingTarget,   "Log(U^{2})",  "Acceptance",              8);
    CosmeticsGraph(fgAccCouplingTAX,      "Log(U^{2})",  "Acceptance",              9);
    CosmeticsGraph(fgYieldCoupling,       "Log(U^{2})",  "Yield per POT",           2);
    CosmeticsGraph(fgYieldCouplingTarget, "Log(U^{2})",  "Yield per POT",           8);
    CosmeticsGraph(fgYieldCouplingTAX,    "Log(U^{2})",  "Yield per POT",           9);
    CosmeticsGraph(fgGammaTotMass,        "N mass [GeV]",     "Total decay width [MeV]", 2);
    CosmeticsGraph(fgTauMass,             "N mass [GeV]",     "Lifetime [ns]",           2);
    CosmeticsGraph(fgAccMass,             "N mass [GeV]",     "Acceptance",              2);
    CosmeticsGraph(fgAccMassTarget,       "N mass [GeV]",     "Acceptance",              8);
    CosmeticsGraph(fgAccMassTAX,          "N mass [GeV]",     "Acceptance",              9);
    CosmeticsGraph(fgYieldMass,           "N mass [GeV]",     "Yield per POT",           2);
    CosmeticsGraph(fgYieldMassTarget,     "N mass [GeV]",     "Yield per POT",           8);
    CosmeticsGraph(fgYieldMassTAX,        "N mass [GeV]",     "Yield per POT",           9);
    CosmeticsGraph(fgAccMom,              "N momentum [GeV]", "Acceptance",              2);
    CosmeticsGraph(fgYieldMom,            "N momentum [GeV]", "Yield per POT",           2);
    CosmeticsGraph(fgExclusion,           "N mass [GeV]",     "Log(U^{2})",         2);
  }
  else {

    PlotErrorBars(fgErrorAccCoupling, fgErrorAccMass, fgErrorYieldCoupling, fgErrorYieldMass, "ErrorBars.txt");
    PlotErrorBars(fgErrorAccCouplingTarget, fgErrorAccMassTarget, fgErrorYieldCouplingTarget, fgErrorYieldMassTarget, "ErrorBarsTarget.txt");
    PlotErrorBars(fgErrorAccCouplingTAX, fgErrorAccMassTAX, fgErrorYieldCouplingTAX, fgErrorYieldMassTAX, "ErrorBarsTAX.txt");
    PlotErrorBarsMom(fgErrorAccMom, fgErrorYieldMom, "ErrorBarsMom.txt");
  }

  SaveAllPlots();

  return;
}

// Cosmetics for plots

void HeavyNeutrinoScan::CosmeticsGraph(TGraphErrors* g, const char* x, const char* y, Int_t col) {

  g->GetXaxis()->SetTitle(x);
  g->GetYaxis()->SetTitle(y);
  g->SetLineColor(col);
  g->SetLineWidth(3);
  g->Draw("AC");
  g->Write();
  
  return;
}

// Compute error bars for graphs

void HeavyNeutrinoScan::PlotErrorBars(TGraphErrors* g, TGraphErrors* g1, TGraphErrors* g2, TGraphErrors* g3, std::string fileName) {

  fstream file;
  std::string line;
  std::string MUpair;
  std::string massStr, couplingStr;
  Double_t Mass, Coupling, Acc, Yield;
  std::map<std::string, std::vector<Double_t>> errorListAcc;
  std::map<std::string, std::vector<Double_t>> errorListYield;
  std::map<std::string, Double_t> errorResAcc;
  std::map<std::string, Double_t> errorResYield;

  file.open(fileName);

  if (!file.is_open()) {
    cout << "Cannot open file " << fileName << " in --histo mode" << endl;
    _exit(1);
  }

  while (!file.eof()) {
    while (getline(file, line)) {
      stringstream s(line);
      s >> Mass >> Coupling >> Acc >> Yield;
      massStr = std::to_string(Mass);
      couplingStr = std::to_string(Coupling);
      massStr.erase(massStr.find_last_not_of('0') + 1, std::string::npos);
      couplingStr.erase(couplingStr.find_last_not_of('0') + 1, std::string::npos);
      MUpair = massStr + " " + couplingStr;
      errorListAcc[MUpair].push_back(Acc);
      errorListYield[MUpair].push_back(Yield);
    }
  }
  
  for (auto it = errorListAcc.begin(); it != errorListAcc.end(); it++) {
    errorResAcc[it->first] = ComputeRMS(it->second)[1];
  }

  for (auto it = errorListYield.begin(); it != errorListYield.end(); it++)
    errorResYield[it->first] = ComputeRMS(it->second)[1];

  for (auto it = errorResAcc.begin(); it != errorResAcc.end(); it++) {                        
    stringstream ss(it->first);                
    ss >> Mass >> Coupling;                                                                 
    if (Mass == fMassForSingleValue) {
      for (Int_t p = 0; p < g->GetN(); p++) {
	Double_t x = 0.;                                                                
	Double_t y = 0.;
	g->GetPoint(p, x, y);
	if (x == Coupling)
	  g->SetPointError(p, 0., errorResAcc[it->first]);
      }
    }
    if (Coupling == (fCouplingStart - fCouplingStop)/2.) {
      for (Int_t p = 0; p < g1->GetN(); p++) {
        Double_t x = 0.;
        Double_t y = 0.;
	g1->GetPoint(p, x, y);
	if (x*1000. == Mass) {
	  g1->SetPointError(p, 0., errorResAcc[it->first]);
	}
      }
    }
  }

  for (auto it = errorResYield.begin(); it != errorResYield.end(); it++) {                        
    stringstream ss(it->first);                                                               
    ss >> Mass >> Coupling;                                                                 
    if (Mass == fMassForSingleValue) {
      for (Int_t p = 0; p < g2->GetN(); p++) {
	Double_t x = 0.;                                                                
	Double_t y = 0.;
	g2->GetPoint(p, x, y);
	if (x == Coupling) 
	  g2->SetPointError(p, 0., errorResYield[it->first]);
      }
    }
    if (Coupling == (fCouplingStart - fCouplingStop)/2.) {
      for (Int_t p = 0; p < g3->GetN(); p++) {
        Double_t x = 0.;
        Double_t y = 0.;
        g3->GetPoint(p, x, y);
        if (x*1000. == Mass)
          g3->SetPointError(p, 0., errorResYield[it->first]);
      }
    }
  }

  return;
}

// Compute error bars for momentum graphs

void HeavyNeutrinoScan::PlotErrorBarsMom(TGraphErrors* g, TGraphErrors* g1, std::string fileName) {
  
  fstream file;
  std::string line;
  std::string momStr;
  Double_t Mom, Acc, Yield;
  std::map<std::string, std::vector<Double_t>> errorListAcc;
  std::map<std::string, std::vector<Double_t>> errorListYield;
  std::map<std::string, Double_t> errorResAcc;
  std::map<std::string, Double_t> errorResYield;

  file.open(fileName);
  
  if (!file.is_open()) {
    cout << "Cannot open file " << fileName << " in --histo mode" << endl;
    _exit(1);
  }
  
  while (!file.eof()) {
    while (getline(file, line)) {
      stringstream s(line);
      s >> Mom >> Acc >> Yield;
      momStr = std::to_string(Mom);
      momStr.erase(momStr.find_last_not_of('0') + 1, std::string::npos);
      errorListAcc[momStr].push_back(Acc);
      errorListYield[momStr].push_back(Yield);
    }
  }

  for (auto it = errorListAcc.begin(); it != errorListAcc.end(); it++)
    errorResAcc[it->first] = ComputeRMS(it->second)[1];

  for (auto it = errorListYield.begin(); it != errorListYield.end(); it++)
    errorResYield[it->first] = ComputeRMS(it->second)[1];

  for (auto it = errorResAcc.begin(); it != errorResAcc.end(); it++) {                        
    stringstream ss(it->first);                                                    
    ss >> Mom;

    for (Int_t p = 0; p < g->GetN(); p++) {
      Double_t x = 0.;                                                                
      Double_t y = 0.;
      g->GetPoint(p, x, y);
      if (x == Mom)
	g->SetPointError(p, 0., errorResAcc[it->first]);
    }
  }
  
  for (auto it = errorResYield.begin(); it != errorResYield.end(); it++) {                        
    stringstream ss(it->first);                                                               
    ss >> Mom;

    for (Int_t p = 0; p < g1->GetN(); p++) {
      Double_t x = 0.;                                                                
      Double_t y = 0.;
      g1->GetPoint(p, x, y);
      if (x == Mom)
	g1->SetPointError(p, 0., errorResYield[it->first]);
    }
  }
  
  return;
}

// Compute mean and RMS of a sample

std::vector<Double_t> HeavyNeutrinoScan::ComputeRMS(std::vector<Double_t> v) {

  Double_t Mean = 0.;
  Double_t RMS  = 0.;
  Double_t Sum  = 0.;
  Double_t Diff = 0.;
  std::vector<Double_t> res;

  for (auto it = v.begin(); it != v.end(); it++) {
    Sum += *it;
  }

  Mean = Sum/v.size();

  for (auto it = v.begin(); it != v.end(); it++)
    Diff += (*it - Mean)*(*it - Mean);
  
  if (v.size() == 1)
    RMS = 0.;
  else 
    RMS = TMath::Sqrt(Diff/(v.size() - 1));

  res.push_back(Mean);
  res.push_back(RMS);

  return res;
}

HeavyNeutrinoScan::~HeavyNeutrinoScan() {
  
  fErrorFile.      close();
  fErrorFileTarget.close();
  fErrorFileTAX.   close();
  fErrorFileMom.   close();

  if (fReadingData) {
    
    fhZMotherProd         = nullptr;
    fhZDProd              = nullptr;
    fhZTauProd            = nullptr;
    fhZDDecay             = nullptr;
    fhDTheta              = nullptr;
    fhDLambda             = nullptr;
    fhDPath               = nullptr;
    fhDMom                = nullptr;
    fhZHNLDecay           = nullptr;
    fhHNLGamma            = nullptr;
    fhHNLDecayProb        = nullptr;
    fhHNLReachProb        = nullptr;
    fhHNLTheta            = nullptr;
    fhHNLMom              = nullptr;
    fhWeight              = nullptr;
    fhCoupling            = nullptr;
    fhMass                = nullptr;
    fhAcc                 = nullptr;
    fhAccTarget           = nullptr;
    fhAccTAX              = nullptr;
    fhYield               = nullptr;
    fhYieldTarget         = nullptr;
    fhYieldTAX            = nullptr;
    fhTransverse          = nullptr;
    fhDThetaMom           = nullptr;
    fhXYDecay             = nullptr;

    fhReachCoupling       = nullptr;
    fhDecayCoupling       = nullptr;
    fhProbCoupling        = nullptr;
    fhWeightCoupling      = nullptr;
    fhEnergyCoupling      = nullptr;
    fgAccCoupling         = nullptr;
    fgAccCouplingTarget   = nullptr;
    fgAccCouplingTAX      = nullptr;
    fgYieldCoupling       = nullptr;
    fgYieldCouplingTarget = nullptr;
    fgYieldCouplingTAX    = nullptr;
    fgGammaTotCoupling    = nullptr;
    fgTauCoupling         = nullptr;

    fhReachMass           = nullptr;
    fhDecayMass           = nullptr;
    fhProbMass            = nullptr;
    fhWeightMass          = nullptr;
    fhEnergyMass          = nullptr;
    fgAccMass             = nullptr;
    fgAccMassTarget       = nullptr;
    fgAccMassTAX          = nullptr;
    fgYieldMass           = nullptr;
    fgYieldMassTarget     = nullptr;
    fgYieldMassTAX        = nullptr;
    fgGammaTotMass        = nullptr;
    fgTauMass             = nullptr;

    fgAccMom              = nullptr;
    fgYieldMom            = nullptr;
    fgExclusion           = nullptr;

    fhpDS                 = nullptr;
    fhpDS1                = nullptr;
    fhptDS                = nullptr;
    fhptDS1               = nullptr;
    fhpNS                 = nullptr;
    fhpNS1                = nullptr;
    fhptNS                = nullptr;
    fhptNS1               = nullptr;
    fhpmuS                = nullptr;
    fhpmuS1               = nullptr;
    fhptmuS               = nullptr;
    fhptmuS1              = nullptr;
    fhpD0                 = nullptr;
    fhpD01                = nullptr;
    fhptD0                = nullptr;
    fhptD01               = nullptr;
    fhpN0                 = nullptr;
    fhpN01                = nullptr;
    fhptN0                = nullptr;
    fhptN01               = nullptr;
    fhpmu0                = nullptr;
    fhpmu01               = nullptr;
    fhptmu0               = nullptr;
    fhptmu01              = nullptr;
  }
  else {
    /*
    remove("ErrorBars.txt");
    remove("ErrorBarsTarget.txt");
    remove("ErrorBarsTAX.txt");
    remove("ErrorBarsMom.txt");
    */
    fgErrorAccCoupling         = nullptr; 
    fgErrorAccCouplingTarget   = nullptr;
    fgErrorAccCouplingTAX      = nullptr;
    fgErrorYieldCoupling       = nullptr;
    fgErrorYieldCouplingTarget = nullptr;
    fgErrorYieldCouplingTAX    = nullptr;
    fgErrorAccMass             = nullptr;
    fgErrorAccMassTarget       = nullptr;
    fgErrorAccMassTAX          = nullptr;
    fgErrorYieldMass           = nullptr;
    fgErrorYieldMassTarget     = nullptr;
    fgErrorYieldMassTAX        = nullptr;
    fgErrorAccMom              = nullptr;
    fgErrorYieldMom            = nullptr;
    fgErrorExclusion           = nullptr;
  }
}
