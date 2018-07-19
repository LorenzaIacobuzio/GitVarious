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
#include "HNLFunctions.hh"
#include "HNLWeight.hh"
#include "HeavyNeutrinoScan.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

/// \class HeavyNeutrinoScan

HeavyNeutrinoScan::HeavyNeutrinoScan(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrinoScan") {

  RequestAllMCTrees();

  AddParam("USquared", &fUSquared, 1.E-6);
  AddParam("UeSquaredRatio", &fUeSquaredRatio, 1.);
  AddParam("UmuSquaredRatio", &fUmuSquaredRatio, 16.);
  AddParam("UtauSquaredRatio", &fUtauSquaredRatio, 3.8);
  AddParam("CouplingStart", &fCouplingStart, -7.);
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

  fNMom = round((std::abs(fMomStop - fMomStart))/fMomStep);
  fN = round((std::abs(fCouplingStop-fCouplingStart))/fCouplingStep);
  fNMass = round((std::abs(fMassStop-fMassStart))/fMassStep);
  fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
  fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
  fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;
}

void HeavyNeutrinoScan::InitHist() {
    
  // One value of mass and coupling

  BookHisto("SingleValue/hZMotherProd",  new TH1D("ZMotherProd", "N mother production point", 100000, -1., 27.));
  BookHisto("SingleValue/hZDProd",       new TH1D("ZDProd",      "D meson production point",  100000, -1., 27.));
  BookHisto("SingleValue/hZTauProd",     new TH1D("ZTauProd",    "#tau production point",     100000, -1., 27.));
  BookHisto("SingleValue/hZDDecay",      new TH1D("ZDDecay",     "N production point",        100000, -1., 27.));
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
  BookHisto("SingleValue/hTransverse",   new TH1D("Transverse",  "N transverse position at decay point", 100, 0., 2.));
  BookHisto("SingleValue/hDThetaMom",    new TH2D("DThetaMom",   "D meson polar angle vs momentum", 100, 0., 150., 50, 0., 0.5));
  BookHisto("SingleValue/hXYDecay",      new TH2D("XYDecay",     "X, Y of N at decay point", 100, -1., 1., 100, -1., 1.));

  BookHisto("SingleValue/hG",            new TH1D("sG",  "", fNMom, fMomStart, fMomStop));
  BookHisto("SingleValue/hA",            new TH1D("sA",  "", fNMom, fMomStart, fMomStop));
  BookHisto("SingleValue/hN",            new TH1D("sN",  "", fNMom, fMomStart, fMomStop));

  BookHisto("SingleValue/ErrorAccMom",   new TGraphAsymmErrors());
  fHisto.GetTGraph("SingleValue/ErrorAccMom")->SetNameTitle("SingleValue/ErrorAccMom", "Acceptance vs N momentum");

  BookHisto("SingleValue/ErrorYieldMom", new TGraphAsymmErrors());
  fHisto.GetTGraph("SingleValue/ErrorYieldMom")->SetNameTitle("SingleValue/ErrorYieldMom", "Yield per POT vs N momentum");

  // Coupling scan 

  BookHisto("CouplingScan/hReachCoupling",  new TH2D("ReachCoupling",  "Probability of N reaching the FV vs coupling",                      fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("CouplingScan/hDecayCoupling",  new TH2D("DecayCoupling",  "Probability of N decaying in the FV vs coupling",                   fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("CouplingScan/hProbCoupling",   new TH2D("ProbCoupling",   "Probability of N reaching and decaying in the FV vs coupling",      fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("CouplingScan/hWeightCoupling", new TH2D("WeightCoupling", "N weight vs coupling",                                              fN, fCouplingStart, fCouplingStop, 10000, 1.E-11, 1.E-7));
  BookHisto("CouplingScan/hEnergyCoupling", new TH2D("EnergyCoupling", "N energy vs coupling",                                              fN, fCouplingStart, fCouplingStop, 100, 0., 100.));

  BookHisto("CouplingScan/hG",  new TH1D("cG",  "", fN, fCouplingStart, fCouplingStop));
  BookHisto("CouplingScan/hA",  new TH1D("cA",  "", fN, fCouplingStart, fCouplingStop));
  BookHisto("CouplingScan/hG1", new TH1D("cG1", "", fN, fCouplingStart, fCouplingStop));
  BookHisto("CouplingScan/hA1", new TH1D("cA1", "", fN, fCouplingStart, fCouplingStop));
  BookHisto("CouplingScan/hG2", new TH1D("cG2", "", fN, fCouplingStart, fCouplingStop));
  BookHisto("CouplingScan/hA2", new TH1D("cA2", "", fN, fCouplingStart, fCouplingStop));
  BookHisto("CouplingScan/hN",  new TH1D("cN",  "", fN, fCouplingStart, fCouplingStop));
  BookHisto("CouplingScan/hN1", new TH1D("cN1", "", fN, fCouplingStart, fCouplingStop));
  BookHisto("CouplingScan/hN2", new TH1D("cN2", "", fN, fCouplingStart, fCouplingStop));

  BookHisto("CouplingScan/MeanCoupling", new TGraph());
  fHisto.GetTGraph("CouplingScan/MeanCoupling")->SetNameTitle("CouplingScan/MeanCoupling", "Mean probability of N reaching and decaying in the FV vs coupling");

  BookHisto("CouplingScan/GammaTotCoupling", new TGraph());
  fHisto.GetTGraph("CouplingScan/GammaTotCoupling")->SetNameTitle("CouplingScan/GammaTotCoupling", "N total decay width vs coupling");

  BookHisto("CouplingScan/TauCoupling", new TGraph());
  fHisto.GetTGraph("CouplingScan/TauCoupling")->SetNameTitle("CouplingScan/TauCoupling", "N mean lifetime vs coupling");

  BookHisto("CouplingScan/ErrorAccCoupling", new TGraphAsymmErrors());
  fHisto.GetTGraph("CouplingScan/ErrorAccCoupling")->SetNameTitle("CouplingScan/ErrorAccCoupling", "All");
    
  BookHisto("CouplingScan/ErrorAccCouplingTarget", new TGraphAsymmErrors());
  fHisto.GetTGraph("CouplingScan/ErrorAccCouplingTarget")->SetNameTitle("CouplingScan/ErrorAccCouplingTarget", "Target");

  BookHisto("CouplingScan/ErrorAccCouplingTAX", new TGraphAsymmErrors());
  fHisto.GetTGraph("CouplingScan/ErrorAccCouplingTAX")->SetNameTitle("CouplingScan/ErrorAccCouplingTAX", "TAX");

  BookHisto("CouplingScan/ErrorYieldCoupling", new TGraphAsymmErrors());
  fHisto.GetTGraph("CouplingScan/ErrorYieldCoupling")->SetNameTitle("CouplingScan/ErrorYieldCoupling", "All");
    
  BookHisto("CouplingScan/ErrorYieldCouplingTarget", new TGraphAsymmErrors());
  fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTarget")->SetNameTitle("CouplingScan/ErrorYieldCouplingTarget", "Target");

  BookHisto("CouplingScan/ErrorYieldCouplingTAX", new TGraphAsymmErrors());
  fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTAX")->SetNameTitle("CouplingScan/ErrorYieldCouplingTAX", "TAX");
   
  // Mass scan

  BookHisto("MassScan/hReachMass",  new TH2D("ReachMass",  "Probability of N reaching the FV vs N mass",                 fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 1000, -0.1, 1.1));
  BookHisto("MassScan/hDecayMass",  new TH2D("DecayMass",  "Probability of N decaying in the FV vs N mass",              fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 1000, -0.1, 1.1));
  BookHisto("MassScan/hProbMass",   new TH2D("ProbMass",   "Probability of N reaching and decaying in the FV vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 1000, -0.1, 1.1));
  BookHisto("MassScan/hWeightMass", new TH2D("WeightMass", "N weight vs N mass",                                         fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 10000, 1.E-13, 1.E-9));
  BookHisto("MassScan/hEnergyMass", new TH2D("EnergyMass", "N energy vs N mass",                                         fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 100, 0., 100.));

  BookHisto("MassScan/hG",  new TH1D("mG",  "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
  BookHisto("MassScan/hA",  new TH1D("mA",  "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
  BookHisto("MassScan/hG1", new TH1D("mG1", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
  BookHisto("MassScan/hA1", new TH1D("mA1", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
  BookHisto("MassScan/hG2", new TH1D("mG2", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
  BookHisto("MassScan/hA2", new TH1D("mA2", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
  BookHisto("MassScan/hN",  new TH1D("mN",  "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
  BookHisto("MassScan/hN1", new TH1D("mN1", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));
  BookHisto("MassScan/hN2", new TH1D("mN2", "", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2.));

  BookHisto("MassScan/MeanMass", new TGraph());
  fHisto.GetTGraph("MassScan/MeanMass")->SetNameTitle("MassScan/MeanMass", "Mean probability of N reaching and decaying in the FV vs N mass");

  BookHisto("MassScan/GammaTotMass", new TGraph());
  fHisto.GetTGraph("MassScan/GammaTotMass")->SetNameTitle("MassScan/GammaTotMass", "N total decay width vs N mass");

  BookHisto("MassScan/TauMass", new TGraph());
  fHisto.GetTGraph("MassScan/TauMass")->SetNameTitle("MassScan/TauMass", "N mean lifetime vs N mass");

  BookHisto("MassScan/ErrorAccMass", new TGraphAsymmErrors());
  fHisto.GetTGraph("MassScan/ErrorAccMass")->SetNameTitle("MassScan/ErrorAccMass", "All");
    
  BookHisto("MassScan/ErrorAccMassTarget", new TGraphAsymmErrors());
  fHisto.GetTGraph("MassScan/ErrorAccMassTarget")->SetNameTitle("MassScan/ErrorAccMassTarget", "Target");

  BookHisto("MassScan/ErrorAccMassTAX", new TGraphAsymmErrors());
  fHisto.GetTGraph("MassScan/ErrorAccMassTAX")->SetNameTitle("MassScan/ErrorAccMassTAX", "TAX");

  BookHisto("MassScan/ErrorYieldMass", new TGraphAsymmErrors());
  fHisto.GetTGraph("MassScan/ErrorYieldMass")->SetNameTitle("MassScan/ErrorYieldMass", "All");
    
  BookHisto("MassScan/ErrorYieldMassTarget", new TGraphAsymmErrors());
  fHisto.GetTGraph("MassScan/ErrorYieldMassTarget")->SetNameTitle("MassScan/ErrorYieldMassTarget", "Target");

  BookHisto("MassScan/ErrorYieldMassTAX", new TGraphAsymmErrors());
  fHisto.GetTGraph("MassScan/ErrorYieldMassTAX")->SetNameTitle("MassScan/ErrorYieldMassTAX", "TAX");
    
  // Total scan

  BookHisto("TotalScan/Exclusion", new TGraphAsymmErrors());    
  fHisto.GetTGraph("TotalScan/Exclusion")->SetNameTitle("TotalScan/Exclusion", "Sensitivity as a function of coupling and N mass");

  // Toy-MC comparison

  BookHisto("ToyMC/DS/hpDS",   new TH1D("pDS",   "D_{S} momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 40, 0., 200.));
  BookHisto("ToyMC/DS/hpDS1",  new TH1D("pDS1",  "D_{S} momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 20, 40., 200.));
  BookHisto("ToyMC/DS/hptDS",  new TH1D("ptDS",  "D_{S} transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 10, 0., 5.));
  BookHisto("ToyMC/DS/hptDS1", new TH1D("ptDS1", "D_{S} transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 6, 0., 3.));

  BookHisto("ToyMC/DS/hpNS",   new TH1D("pNS",   "N momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 50, 0., 300.));
  BookHisto("ToyMC/DS/hpNS1",  new TH1D("pNS1",  "N momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 50, 0., 300.));
  BookHisto("ToyMC/DS/hptNS",  new TH1D("ptNS",  "N transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 20, 0., 2.));
  BookHisto("ToyMC/DS/hptNS1", new TH1D("ptNS1", "N transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 20, 0., 2.));

  BookHisto("ToyMC/DS/hpmuS",   new TH1D("pmuS",   "Muon daughter momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 50, 0., 300.));
  BookHisto("ToyMC/DS/hpmuS1",  new TH1D("pmuS1",  "Muon daughter momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 50, 0., 300.));
  BookHisto("ToyMC/DS/hptmuS",  new TH1D("ptmuS",  "Muon daughter transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, all HNLs", 20, 0., 2.));
  BookHisto("ToyMC/DS/hptmuS1", new TH1D("ptmuS1", "Muon daughter transverse momentum from D_{S} #rightarrow N#mu, N #rightarrow #pi#mu, good HNLs only", 20, 0., 2.));

  BookHisto("ToyMC/D0/hpD0",   new TH1D("pD0",   "D^{0} momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 40, 0., 200.));
  BookHisto("ToyMC/D0/hpD01",  new TH1D("pD01",  "D^{0} momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 40, 60., 300.));
  BookHisto("ToyMC/D0/hptD0",  new TH1D("ptD0",  "D^{0} transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 10, 0., 5.));
  BookHisto("ToyMC/D0/hptD01", new TH1D("ptD01", "D^{0} transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 6, 0., 3.));

  BookHisto("ToyMC/D0/hpN0",   new TH1D("pN0",   "N momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 50, 0., 300.));
  BookHisto("ToyMC/D0/hpN01",  new TH1D("pN01",  "N momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 50, 0., 300.));
  BookHisto("ToyMC/D0/hptN0",  new TH1D("ptN0",  "N transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 20, 0., 2.));
  BookHisto("ToyMC/D0/hptN01", new TH1D("ptN01", "N transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 20, 0., 2.));

  BookHisto("ToyMC/D0/hpmu0",   new TH1D("pmu0",   "Muon daughter momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 50, 0., 300.));
  BookHisto("ToyMC/D0/hpmu01",  new TH1D("pmu01",  "Muon daughter momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 50, 0., 300.));
  BookHisto("ToyMC/D0/hptmu0",  new TH1D("ptmu0",  "Muon daughter transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, all HNLs", 20, 0., 2.));
  BookHisto("ToyMC/D0/hptmu01", new TH1D("ptmu01", "Muon daughter transverse momentum from D^{0} #rightarrow K#muN, N #rightarrow #pi#mu, good HNLs only", 20, 0., 2.));

  // Needed for error graphs

  fHisto.GetTH1("SingleValue/hG")->Sumw2();
  fHisto.GetTH1("SingleValue/hA")->Sumw2();
  fHisto.GetTH1("SingleValue/hN")->Sumw2();
  fHisto.GetTH1("CouplingScan/hG")->Sumw2();
  fHisto.GetTH1("CouplingScan/hG1")->Sumw2();
  fHisto.GetTH1("CouplingScan/hG2")->Sumw2();
  fHisto.GetTH1("CouplingScan/hA")->Sumw2();
  fHisto.GetTH1("CouplingScan/hA1")->Sumw2();
  fHisto.GetTH1("CouplingScan/hA2")->Sumw2();
  fHisto.GetTH1("CouplingScan/hN")->Sumw2();
  fHisto.GetTH1("CouplingScan/hN1")->Sumw2();
  fHisto.GetTH1("CouplingScan/hN2")->Sumw2();
  fHisto.GetTH1("MassScan/hG")->Sumw2();
  fHisto.GetTH1("MassScan/hG1")->Sumw2();
  fHisto.GetTH1("MassScan/hG2")->Sumw2();
  fHisto.GetTH1("MassScan/hA")->Sumw2();
  fHisto.GetTH1("MassScan/hA1")->Sumw2();
  fHisto.GetTH1("MassScan/hA2")->Sumw2();
  fHisto.GetTH1("MassScan/hN")->Sumw2();
  fHisto.GetTH1("MassScan/hN1")->Sumw2();
  fHisto.GetTH1("MassScan/hN2")->Sumw2();
}

void HeavyNeutrinoScan::Process(Int_t) {
  
  Double_t fCoupling  = -999.;  
  Double_t MN         = 0.;
  Double_t HNLTau     = 0.;
  Double_t gammaTot   = 0.;
  Double_t momN       = 0.;
  Double_t momBin     = 0.;
  Double_t NDecayProb = 0.;
  Double_t NReachProb = 0.;
  Double_t Weight     = 0.;
  Double_t scale      = 1.E10;
  Double_t N = 0;
  Double_t D = 0;
  Bool_t isGood       = false;

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
	fGammaTot[round(MN)][fCoupling] = gammaTot;
	fTau[round(MN)][fCoupling] = HNLTau;
	fMasses[round(MN)] = round(MN);

	if (fNevents[round(MN)].count(fCoupling) == 0)
	  fNevents[round(MN)][fCoupling] = 0;
	fNevents[round(MN)][fCoupling]++;
	
	if (IsHNLGood == true && isGood == true)
	  fSumGood[round(MN)][fCoupling] += Weight;

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

	  FillHisto("CouplingScan/hReachCoupling",  fCoupling, NReachProb);
	  FillHisto("CouplingScan/hDecayCoupling",  fCoupling, NDecayProb);
	  FillHisto("CouplingScan/hProbCoupling",   fCoupling, NReachProb*NDecayProb);
	  FillHisto("CouplingScan/hWeightCoupling", fCoupling, /*Weight*/1);	
	  FillHisto("CouplingScan/hEnergyCoupling", fCoupling, momN/1000.*momN/1000. + MN/1000.*MN/1000.);
	}

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
	  
	  FillHisto("MassScan/hReachMass",  round(MN)/1000., NReachProb);
	  FillHisto("MassScan/hDecayMass",  round(MN)/1000., NDecayProb);
	  FillHisto("MassScan/hProbMass",   round(MN)/1000., NReachProb*NDecayProb);
	  FillHisto("MassScan/hWeightMass", round(MN)/1000., Weight);
	  FillHisto("MassScan/hEnergyMass", round(MN)/1000., momN/1000.*momN/1000. + MN/1000.*MN/1000.);
	}

	if (fCoupling == fCouplingForSingleValue && round(MN)/1000. == fMassForSingleValue) {
	  FillHisto("SingleValue/hCoupling",     fCoupling);
	  FillHisto("SingleValue/hMass",         round(MN)/1000.);
	  FillHisto("SingleValue/hHNLReachProb", NReachProb);
	  FillHisto("SingleValue/hHNLDecayProb", NDecayProb);
	  FillHisto("SingleValue/hWeight",       Weight);
	  D++;
	  if (IsHNLGood == true && isGood == true) 
	    N += Weight;
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

        if ((TMath::Sqrt(p->GetInitial4Momentum().E()*p->GetInitial4Momentum().E() - p->GetInitial4Momentum().P()*p->GetInitial4Momentum().P()))/1000. == fMassForSingleValue) {
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

  fHisto.GetTH2("SingleValue/hDThetaMom")->GetYaxis()->SetTitle("Polar angle [rad]");
  fHisto.GetTH2("SingleValue/hXYDecay")->GetYaxis()->SetTitle("Y [m]");
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
  
  // Acceptance and yield computation                                             

  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCoupling"))->Divide(fHisto.GetTH1("CouplingScan/hG"), fHisto.GetTH1("CouplingScan/hA"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingTarget"))->Divide(fHisto.GetTH1("CouplingScan/hG1"), fHisto.GetTH1("CouplingScan/hA1"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingTAX"))->Divide(fHisto.GetTH1("CouplingScan/hG2"), fHisto.GetTH1("CouplingScan/hA2"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCoupling"))->Divide(fHisto.GetTH1("CouplingScan/hG"), fHisto.GetTH1("CouplingScan/hN"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTarget"))->Divide(fHisto.GetTH1("CouplingScan/hG1"), fHisto.GetTH1("CouplingScan/hN1"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTAX"))->Divide(fHisto.GetTH1("CouplingScan/hG2"), fHisto.GetTH1("CouplingScan/hN2"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMass"))->Divide(fHisto.GetTH1("MassScan/hG"), fHisto.GetTH1("MassScan/hA"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassTarget"))->Divide(fHisto.GetTH1("MassScan/hG1"), fHisto.GetTH1("MassScan/hA1"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassTAX"))->Divide(fHisto.GetTH1("MassScan/hG2"), fHisto.GetTH1("MassScan/hA2"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMass"))->Divide(fHisto.GetTH1("MassScan/hG"), fHisto.GetTH1("MassScan/hN"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTarget"))->Divide(fHisto.GetTH1("MassScan/hG1"), fHisto.GetTH1("MassScan/hN1"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTAX"))->Divide(fHisto.GetTH1("MassScan/hG2"), fHisto.GetTH1("MassScan/hN2"), "cl=0.683 mode");    
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMom"))->Divide(fHisto.GetTH1("SingleValue/hG"), fHisto.GetTH1("SingleValue/hA"), "cl=0.683 mode");
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorYieldMom"))->Divide(fHisto.GetTH1("SingleValue/hG"), fHisto.GetTH1("SingleValue/hN"), "cl=0.683 mode");

  Double_t Coupling     = 0.;
  Double_t MN           = 0.;
  Int_t couplingCounter = 0;
  Int_t massCounter     = 0;

  for (auto it = fMasses.begin(); it != fMasses.end(); it++) {
    MN = it->first;
    for (auto it1 = fCouplings.begin(); it1 != fCouplings.end(); it1++) {
      Coupling = it1->first;

      if (MN/1000. == fMassForSingleValue) {
	fHisto.GetTGraph("CouplingScan/GammaTotCoupling")->SetPoint(couplingCounter, Coupling, fGammaTot[fMassForSingleValue*1000.][Coupling]);
	fHisto.GetTGraph("CouplingScan/TauCoupling")->SetPoint(couplingCounter, Coupling, fTau[fMassForSingleValue*1000.][Coupling]);
	couplingCounter++;
      }

      fYield[MN][Coupling] = fSumGood[MN][Coupling]/fNevents[MN][Coupling];
    }
    
    Coupling = fCouplingForSingleValue;
    fHisto.GetTGraph("MassScan/GammaTotMass")->SetPoint(massCounter, MN/1000., fGammaTot[MN][fCouplingForSingleValue]);
    fHisto.GetTGraph("MassScan/TauMass")->SetPoint(massCounter, MN/1000., fTau[MN][fCouplingForSingleValue]);
    massCounter++;
  }

  // Sensitivity

  Double_t POT = 1.E18;
  Double_t M = -999.;
  Double_t Umax = 999.;
  Double_t Umin = 999.;
  Double_t min = 1.E50;
  Double_t max = 1.E-50;
  Double_t mean = 0.;
  Double_t startM = 0.;
  Double_t stopM = 0.;
  Int_t meanN = 0;
  Int_t count = 0;
  Int_t counter = 1;

  for (auto it = fMasses.begin(); it != fMasses.end(); it++) {
    MN = it->first;
    for (auto it1 = fCouplings.begin(); it1 != fCouplings.end(); it1++) {
      Coupling = it1->first;
      if (fYield[MN][Coupling]*POT > 2.3) {
	M = MN;
	mean += -Coupling;
	meanN++;
      }
      cout<<M<<" "<<mean<<" "<<min<<" "<<Umin<<" "<<max<<" "<<Umax<<" "<<fYield[MN][Coupling]*POT<<endl;
      if (fYield[MN][Coupling]*POT > max && fYield[MN][Coupling]*POT > 2.3) {
	Umax = Coupling;
	max = fYield[MN][Coupling]*POT;
      }
      if (fYield[MN][Coupling]*POT < min && fYield[MN][Coupling]*POT > 2.3) {
	Umin = Coupling;
	min = fYield[MN][Coupling]*POT;
      }
    }

    if (meanN != 0.) {
      mean = -mean/meanN;
    }

    if (-mean > 0.) {
      count++;
      if (count == 1) {
	startM = M;
	static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("TotalScan/Exclusion"))->SetPoint(0, startM/1000., fCouplingStop);
      }
      static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("TotalScan/Exclusion"))->SetPoint(counter, M/1000., mean);
      static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("TotalScan/Exclusion"))->SetPointError(counter, 0., 0., TMath::Abs(mean-Umin), TMath::Abs(Umax-mean));
    //cout<<M<<" "<<mean<<" "<<Umin<<" "<<TMath::Abs(mean-Umin)<<" "<<Umax<<" "<<TMath::Abs(Umax-mean)<<endl;
      counter++;
      stopM = M;
      mean = 0.;
      meanN = 0;
      min = 1.E50;
      max = 1.E-50;
      M = -999.;
      Umax = 999.;
      Umin = 999.;
    }
    mean = 0.;
    meanN = 0;
    min = 1.E50;
    max = 1.E-50;
    M = -999.;
    Umax = 999.;
    Umin = 999.;
  }
  
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("TotalScan/Exclusion"))->SetPoint(counter, stopM/1000., fCouplingStop);
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("TotalScan/Exclusion"))->SetPoint(counter+1, startM/1000., fCouplingStop);
  static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("TotalScan/Exclusion"))->Print();

  // Mean probability vs mass and coupling

  for (Int_t i = 0; i < fHisto.GetTH2("CouplingScan/hProbCoupling")->GetNbinsX()-1; i++) {
    TH1D *h = fHisto.GetTH2("CouplingScan/hProbCoupling")->ProjectionY("", i, i+1, "");
    fHisto.GetTGraph("CouplingScan/MeanCoupling")->SetPoint(i, fHisto.GetTH2("CouplingScan/hProbCoupling")->GetXaxis()->GetBinCenter(i), h->GetMean());
  }
  
  for (Int_t i = 0; i < fHisto.GetTH2("MassScan/hProbMass")->GetNbinsX()-1; i++) {      
    TH1D *h = fHisto.GetTH2("MassScan/hProbMass")->ProjectionY("", i, i+1, "");           
    fHisto.GetTGraph("MassScan/MeanMass")->SetPoint(i, fHisto.GetTH2("MassScan/hProbMass")->GetXaxis()->GetBinCenter(i), h->GetMean());
  }  
  
  // Cosmetics
  
  CosmeticsGraph(fHisto.GetTGraph("CouplingScan/GammaTotCoupling"), "Log(U^{2})", "Total decay width [MeV]", 2);
  CosmeticsGraph(fHisto.GetTGraph("CouplingScan/TauCoupling"),  "Log(U^{2})", "Lifetime [ns]", 2);
  CosmeticsGraph(fHisto.GetTGraph("CouplingScan/MeanCoupling"),  "Log(U^{2})", "Mean probability", 2);
  CosmeticsGraph(fHisto.GetTGraph("MassScan/GammaTotMass"),  "N mass [GeV/c^{2}]", "Total decay width [MeV]", 2);
  CosmeticsGraph(fHisto.GetTGraph("MassScan/TauMass"),  "N mass [GeV/c^{2}]", "Lifetime [ns]", 2);
  CosmeticsGraph(fHisto.GetTGraph("MassScan/MeanMass"),  "N mass [GeV/c^{2}]", "Mean probability", 2);
  
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCoupling")), "Log(U^{2})", "Acceptance", 2);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingTarget")), "Log(U^{2})", "Acceptance", 8);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorAccCouplingTAX")), "Log(U^{2})", "Acceptance", 9);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCoupling")), "Log(U^{2})", "Yield per POT", 2);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTarget")), "Log(U^{2})", "Yield per POT", 8);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTAX")), "Log(U^{2})", "Yield per POT", 9);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMass")), "N mass [GeV/c^{2}]", "Acceptance", 2);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassTarget")), "N mass [GeV/c^{2}]", "Acceptance", 8);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassTAX")), "N mass [GeV/c^{2}]", "Acceptance", 9);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMass")), "N mass [GeV/c^{2}]", "Yield per POT", 2);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTarget")), "N mass [GeV/c^{2}]", "Yield per POT", 8);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTAX")), "N mass [GeV/c^{2}]", "Yield per POT", 9);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMom")), "P [GeV/c]", "Acceptance", 2);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorYieldMom")), "P [GeV/c]", "Yield per POT", 2);
  CosmeticsGraph( static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("TotalScan/Exclusion")), "N mass [GeV/c^{2}]", "Log(U^{2})", 2);

  SaveAllPlots();
  
  return;
}

// Cosmetics for plots

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
      
  if (!strcmp(g->GetName(), "TotalScan/Exclusion")) {
    g->SetLineColor(col);
    g->SetLineWidth(1);
    //g->SetFillColor(col); 
    //g->SetFillStyle(3001); 
    //g->Draw("ALF");
    g->Draw("AP");
  }
  else {
    g->SetLineColor(1);
    g->SetLineWidth(1);
    g->SetMarkerStyle(8);
    g->SetMarkerColor(col);
    g->Draw("AP");
  }
  
  return;
}
