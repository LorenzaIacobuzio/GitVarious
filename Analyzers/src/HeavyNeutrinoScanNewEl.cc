// ---------------------------------------------------------------        
//                                                                                  
// History:                                                                         
//                                                                             
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018    
//                                                                                
// ---------------------------------------------------------------              
/// \class HeavyNeutrinoScanNewEl
/// \Brief                                                                       
/// Produce expected sensitivity curves for HNL MC samples of different masses, according to the chosen theoretical model
/// \EndBrief                                                                  
//
/// \Detailed                                                                        
/// After retrieving all HNL weights using the ToolsLib called HNLWeight, and checking if an event passed the selection implemented in the analyzer HeavyNeutrino, the analyzer produces plots necessary to the sensitivity computation. Running the analyzer a second time in histo mode generates the final curve and all the plots listed below (written in the ROOT file in different directories): 
/// SingleValue: for all HNL masses together and a fixed coupling value (default is -6 = Log(10^-6)), quantities related to the D meson, the HNL and its product are produced; for a fixed HNL mass (default is 1 GeV) and coupling (default is -6), acceptances and yield per proton on traget (POT) are also computed as a function of the HNL momentum;
/// CouplingScan: for a fixed HNL mass, several plots, including acceptances and yield, are produced as a function of the coupling; the initial and final values of the coupling and the step of the coupling scan can be set as parameters of the analyzer (CouplingStart, CouplingStop, CouplingStep);
/// MassScan: for a fixed coupling, several plots, including acceptances and yield, are produced as a function of the HNL mass; the initial and final values of the mass and the step of the mass scan can be set as parameters of the analyzer, in GeV, and must correspond to the same parameters set at the MC event generation stage (MassStart, MassStop, MassStep);
/// TotalScan: the expected sensitivity curve is produced as a function of the HNL mass and coupling, according to the chosen theoretical model;
/// ToyMC: several plots are produced for toy-MC comparison, for two decay chains, DS->Nmu; N->pimu, and D0->KNmu; N->pimu.
///                  
/// This analyzer makes use of two ToolsLib, called HNLFunctions and HNLWeight.
/// Some parameters can be set: the values of the ratios between specific-flavour couplings (UeSquaredRatio, UmuSquaredRatio, UtauSquaredRatio); the values of the initial and final momentum and its step, used for the plots in SingleValue (MomStart, MomStop, MomStep); length of the beginning Z of the fiducial volume, length of the FV, decay mode of the HNLs, as set in the MC macro (InitialFV, LFV, Mode); Mode = 0 corresponds to pi-mu states, 1 to pi-e, 2 to rho-mu, 3 to rho-e.
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
#include <TRandom3.h>
#include "MCSimple.hh"
#include "MCInfo.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "BeamParameters.hh"
#include "GeometricAcceptance.hh"
#include "HNLFunctions.hh"
#include "HNLWeight.hh"
#include "HeavyNeutrinoScanNewEl.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

/// \class HeavyNeutrinoScanNewEl

HeavyNeutrinoScanNewEl::HeavyNeutrinoScanNewEl(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrinoScanNewEl") {

  fReadingData = GetIsTree();

  RequestAllMCTrees();

  AddParam("USquared", &fUSquared, 1.E-6); // change accordingly
  AddParam("InitialUeSquaredRatio", &fInitialUeSquaredRatio, 1.); // change accordingly
  AddParam("InitialUmuSquaredRatio", &fInitialUmuSquaredRatio, 1.); // change accordingly
  AddParam("InitialUtauSquaredRatio", &fInitialUtauSquaredRatio, 0.); // change accordingly
  AddParam("CouplingStart", &fCouplingStart, -10.); // -10
  AddParam("CouplingStop", &fCouplingStop, -1.); // -1 (do not put 0)
  AddParam("CouplingStep", &fCouplingStep, 0.1); // 0.1
  AddParam("MassStart", &fMassStart, 0.250); // keep this min
  AddParam("MassStop", &fMassStop, 1.960); // keep this max
  AddParam("MassStep", &fMassStep, 0.01); // keep this step
  AddParam("InitialFV", &fInitialFV, 102425.); // keep
  AddParam("LFV", &fLFV, 77575.); // keep
  AddParam("Mode", &fMode, 0);
  AddParam("MomStop", &fMomStop, 300.);
  AddParam("MomStart", &fMomStart, 0.);
  AddParam("MomStep", &fMomStep, 10.);
  AddParam("MassForSingleValue", &fMassForSingleValue, 1.);
  AddParam("CouplingForSingleValue", &fCouplingForSingleValue, -6.);

  fUeSquaredRatio = fInitialUeSquaredRatio;
  fUmuSquaredRatio = fInitialUmuSquaredRatio;
  fUtauSquaredRatio = fInitialUtauSquaredRatio;
  fNMom = round((std::abs(fMomStop - fMomStart))/fMomStep);
  fN = round((std::abs(fCouplingStop-fCouplingStart))/fCouplingStep);
  fNMass = round((std::abs(fMassStop-fMassStart))/fMassStep);
  fNContours = 20;
  r = new TRandom3(0);

  fCDAcomp = new TwoLinesCDA();
  fDistcomp = new PointLineDistance();
}

void HeavyNeutrinoScanNewEl::InitOutput() {

  OpenNewTree("Scan", "Scan");

  AddBranch("Scan", "Mass", &fTMass);
  AddBranch("Scan", "Coupling", &fTCoupling);
  AddBranch("Scan", "NEvents", &fTNEvents);
  AddBranch("Scan", "SumGood", &fTSumGood);
}

void HeavyNeutrinoScanNewEl::InitHist() {

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
    BookHisto("SingleValue/hHNLTheta",      new TH1D("HNLTheta", "N polar angle", 100., 0., 0.5));
    BookHisto("SingleValue/hHNLMom",        new TH1D("HNLMom", "N momentum", 100., -0.5, 200.));
    BookHisto("SingleValue/hTransverse",    new TH1D("Transverse", "N transverse position at decay point", 30, 0., 2.));
    BookHisto("SingleValue/hDProd",         new TH2D("DProd", "D meson production point in the TAXes", 50, -20., 20., 50, -35., -10.));
    BookHisto("SingleValue/hXYProd",        new TH2D("XYProd", "N production point", 40, -20., 20., 40, -35., 5.));
    BookHisto("SingleValue/hDThetaMom",     new TH2D("DThetaMom", "D meson polar angle vs momentum", 300, 0., 150., 50, 0., 0.5));
    BookHisto("SingleValue/hXYDecay",       new TH2D("XYDecay", "X, Y of N at decay point", 20, -1., 1., 20, -1., 1.));
    BookHisto("SingleValue/hXYDecay2Body",  new TH2D("XYDecay2Body", "X, Y of N at decay point", 20, -1., 1., 20, -1., 1.));
    BookHisto("SingleValue/hXYDecay3Body",  new TH2D("XYDecay3Body", "X, Y of N at decay point", 20, -1., 1., 20, -1., 1.));
    BookHisto("SingleValue/hXYDecayTom",    new TH2D("XYDecayTom", "X, Y of N at decay point", 20, -1., 1., 20, -1., 1.));
    BookHisto("SingleValue/hXYMuon",        new TH2D("XYMuon", "X, Y of muon daughter at CH1", 100, -1., 1., 100, -1., 1.));
    BookHisto("SingleValue/hXYPion",        new TH2D("XYPion", "X, Y of pion daughter at CH1", 100, -1., 1., 100, -1., 1.));
    BookHisto("SingleValue/hBeamvsTarTrue", new TH2D("BeamvsTarTrue", "N true trajectory (target)", 20, 0., 0.2, 100, 0., 1.)); 
    BookHisto("SingleValue/hBeamvsTAXTrue", new TH2D("BeamvsTAXTrue", "N true trajectory (TAX)", 20, 0., 0.2, 100, 0., 1.)); 

    BookHisto("SingleValue/hG", new TH1D("sG", "", fNMom, fMomStart, fMomStop));
    BookHisto("SingleValue/hR", new TH1D("sR", "", fNMom, fMomStart, fMomStop));
    BookHisto("SingleValue/hA", new TH1D("sA", "", fNMom, fMomStart, fMomStop));
    BookHisto("SingleValue/hN", new TH1D("sN", "", fNMom, fMomStart, fMomStop));

    // Coupling scan 

    BookHisto("CouplingScan/hReachCoupling",  new TH2D("ReachCoupling", "Probability of N reaching the FV vs coupling", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2., 100, -0.1, 1.1));
    BookHisto("CouplingScan/hDecayCoupling",  new TH2D("DecayCoupling", "Probability of N decaying in the FV vs coupling", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2., 100, -0.1, 1.1));
    BookHisto("CouplingScan/hProbCoupling",   new TH2D("ProbCoupling", "Probability of N reaching and decaying in the FV vs coupling", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2., 20, -20., 0.));
    BookHisto("CouplingScan/hWeightCoupling", new TH2D("WeightCoupling", "N weight vs coupling", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2., 100, -25., -7.));
    BookHisto("CouplingScan/hGoodWeightCoupling", new TH2D("GoodWeightCoupling", "Good N weight vs coupling", fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2., 1000, 1.E-12, 1.E-8));
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

    BookHisto("MassScan/hReachMass",     new TH2D("ReachMass", "Probability of N reaching the FV vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 100, -0.1, 1.1));
    BookHisto("MassScan/hDecayMass",     new TH2D("DecayMass", "Probability of N decaying in the FV vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 100, -0.1, 1.1));
    BookHisto("MassScan/hProbMass",      new TH2D("ProbMass", "Probability of N reaching and decaying in the FV vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 20, -8., 0.));
    BookHisto("MassScan/hWeightMass",    new TH2D("WeightMass", "N weight vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 20, -25., -10.));
    BookHisto("MassScan/hGoodWeightMass",    new TH2D("GoodWeightMass", "Good N weight vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 10000, 1.E-13, 1.E-9));
    BookHisto("MassScan/hEnergyMass",    new TH2D("EnergyMass", "N energy vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 100, 0., 100.));
    BookHisto("MassScan/hEnergyMassTom", new TH2D("EnergyMassTom", "N energy vs N mass", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., 100, 0., 100.));

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
    ImportAllInputHistogram("HeavyNeutrinoScanNewEl/SingleValue", false, "SingleValue");
    ImportAllInputHistogram("HeavyNeutrinoScanNewEl/CouplingScan", false, "CouplingScan");
    ImportAllInputHistogram("HeavyNeutrinoScanNewEl/MassScan", false, "MassScan");
    ImportAllInputHistogram("HeavyNeutrinoScanNewEl/ToyMC/DS", false, "ToyMC/DS");
    ImportAllInputHistogram("HeavyNeutrinoScanNewEl/ToyMC/D0", false, "ToyMC/D0");
    
    // One value
    /*
    BookHisto("SingleValue/ErrorAccMomSel", new TGraphAsymmErrors());
    fHisto.GetTGraph("SingleValue/ErrorAccMomSel")->SetNameTitle("ErrorAccMomSel", "Selection acceptance vs N momentum");
    BookHisto("SingleValue/ErrorAccMomReg", new TGraphAsymmErrors());
    fHisto.GetTGraph("SingleValue/ErrorAccMomReg")->SetNameTitle("ErrorAccMomReg", "Regeneration acceptance vs N momentum");
    BookHisto("SingleValue/ErrorAccMomFV", new TGraphAsymmErrors());
    fHisto.GetTGraph("SingleValue/ErrorAccMomFV")->SetNameTitle("ErrorAccMomFV", "FV acceptance vs N momentum");

    BookHisto("SingleValue/ErrorYieldMom", new TGraphAsymmErrors());
    fHisto.GetTGraph("SingleValue/ErrorYieldMom")->SetNameTitle("ErrorYieldMom", "Yield per POT vs N momentum");
    */
    // Coupling

    BookHisto("CouplingScan/MeanCoupling", new TGraph());
    fHisto.GetTGraph("CouplingScan/MeanCoupling")->SetNameTitle("MeanCoupling", "Mean probability of N reaching and decaying in the FV vs coupling");
        
    BookHisto("CouplingScan/GammaTotCoupling", new TGraph());
    fHisto.GetTGraph("CouplingScan/GammaTotCoupling")->SetNameTitle("GammaTotCoupling", "N total decay width vs coupling");
    BookHisto("CouplingScan/TauCoupling", new TGraph());
    fHisto.GetTGraph("CouplingScan/TauCoupling")->SetNameTitle("TauCoupling", "N mean lifetime vs coupling");

    BookHisto("CouplingScan/ErrorAccCouplingSelTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSelTarget")->SetNameTitle("ErrorAccCouplingSelTarget", "Selection Target");
    BookHisto("CouplingScan/ErrorAccCouplingRegTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingRegTarget")->SetNameTitle("ErrorAccCouplingRegTarget", "Regeneration Target");
    BookHisto("CouplingScan/ErrorAccCouplingFVTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFVTarget")->SetNameTitle("ErrorAccCouplingFVTarget", "FV Target");
    BookHisto("CouplingScan/ErrorYieldCouplingTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTarget")->SetNameTitle("ErrorYieldCouplingTarget", "Target");

    BookHisto("CouplingScan/ErrorAccCouplingSelTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSelTAX")->SetNameTitle("ErrorAccCouplingSelTAX", "Selection TAX");
    BookHisto("CouplingScan/ErrorAccCouplingRegTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingRegTAX")->SetNameTitle("ErrorAccCouplingRegTAX", "Regeneration TAX");
    BookHisto("CouplingScan/ErrorAccCouplingFVTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFVTAX")->SetNameTitle("ErrorAccCouplingFVTAX", "FV TAX");
    BookHisto("CouplingScan/ErrorYieldCouplingTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTAX")->SetNameTitle("ErrorYieldCouplingTAX", "TAX");

    BookHisto("CouplingScan/ErrorAccCouplingSel", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingSel")->SetNameTitle("ErrorAccCouplingSel", "Selection");
    BookHisto("CouplingScan/ErrorAccCouplingReg", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingReg")->SetNameTitle("ErrorAccCouplingReg", "Regeneration");
    BookHisto("CouplingScan/ErrorAccCouplingFV", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorAccCouplingFV")->SetNameTitle("ErrorAccCouplingFV", "FV");
    BookHisto("CouplingScan/ErrorYieldCoupling", new TGraphAsymmErrors());
    fHisto.GetTGraph("CouplingScan/ErrorYieldCoupling")->SetNameTitle("ErrorYieldCoupling", "All");

    // Mass

    BookHisto("MassScan/MeanMass", new TGraph());
    fHisto.GetTGraph("MassScan/MeanMass")->SetNameTitle("MeanMass", "Mean probability of N reaching and decaying in the FV vs N mass");
    BookHisto("MassScan/GammaTotMass", new TGraph());
    fHisto.GetTGraph("MassScan/GammaTotMass")->SetNameTitle("GammaTotMass", "N total decay width vs N mass");
    BookHisto("MassScan/TauMass", new TGraph());
    fHisto.GetTGraph("MassScan/TauMass")->SetNameTitle("TauMass", "N mean lifetime vs N mass");

    BookHisto("MassScan/ErrorAccMassSelTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassSelTarget")->SetNameTitle("ErrorAccMassSelTarget", "Selection Target");
    BookHisto("MassScan/ErrorAccMassRegTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassRegTarget")->SetNameTitle("ErrorAccMassRegTarget", "Regeneration Target");
    BookHisto("MassScan/ErrorAccMassFVTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassFVTarget")->SetNameTitle("ErrorAccMassFVTarget", "FV Target");
    BookHisto("MassScan/ErrorYieldMassTarget", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorYieldMassTarget")->SetNameTitle("ErrorYieldMassTarget", "Target");

    BookHisto("MassScan/ErrorAccMassSelTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassSelTAX")->SetNameTitle("ErrorAccMassSelTAX", "Selection TAX");
    BookHisto("MassScan/ErrorAccMassRegTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassRegTAX")->SetNameTitle("ErrorAccMassRegTAX", "Regeneration TAX");
    BookHisto("MassScan/ErrorAccMassFVTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassFVTAX")->SetNameTitle("ErrorAccMassFVTAX", "FV TAX");
    BookHisto("MassScan/ErrorYieldMassTAX", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorYieldMassTAX")->SetNameTitle("ErrorYieldMassTAX", "TAX");

    BookHisto("MassScan/ErrorAccMassSel", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassSel")->SetNameTitle("ErrorAccMassSel", "Selection");
    BookHisto("MassScan/ErrorAccMassReg", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassReg")->SetNameTitle("ErrorAccMassReg", "Regeneration");
    BookHisto("MassScan/ErrorAccMassFV", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorAccMassFV")->SetNameTitle("ErrorAccMassFV", "FV");
    BookHisto("MassScan/ErrorYieldMass", new TGraphAsymmErrors());
    fHisto.GetTGraph("MassScan/ErrorYieldMass")->SetNameTitle("ErrorYieldMass", "All");

    // Total

    BookHisto("TotalScan/hExclusion", new TH2D("Exclusion", "Sensitivity as a function of N mass and coupling", fNMass+1, fMassStart-fMassStep/2., fMassStop+fMassStep/2., fN+1, fCouplingStart-fCouplingStep/2., fCouplingStop+fCouplingStep/2.));

    BookHisto("TotalScan/Contours", new TGraph());
    fHisto.GetTGraph("TotalScan/Contours")->SetNameTitle("Contours", "Contours for sensitivity");
  }
}

void HeavyNeutrinoScanNewEl::Process(Int_t) {

  if (fReadingData) {

    Double_t fCoupling = -999.;  
    Double_t MN = 0.;
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
	  momN = Weights[i]["Momentum"];
	  NReachProb = Weights[i]["ReachProb"];
	  NDecayProb = Weights[i]["DecayProb"];
	  Weight = Weights[i]["Weight"];   
	  isGood = Weights[i]["IsGood"];
	  fMasses[round(MN)] = round(MN);

	  // Maps for exclusion plot

	  if (fNEvents[round(MN)].count(fCoupling) == 0)
	    fNEvents[round(MN)][fCoupling] = 0;
	  fNEvents[round(MN)][fCoupling]++;

	  if (IsHNLGood == true && isGood == true) {
	    fSumGood[round(MN)][fCoupling] += Weight;
	  }

	  // Histos for acceptance and yield: coupling

	  if (round(MN)/1000. == fMassForSingleValue) {

	    // Scan on the N momentum

	    if (fCoupling == fCouplingForSingleValue) {
	      if (momN/1000 >= fMomStart && momN/1000 <= fMomStop) {
		Double_t momN1 = momN/1000.;
		momBin = fMomStep*trunc(momN1/fMomStep);
		
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
	    FillHisto("CouplingScan/hProbCoupling", fCoupling, TMath::Log10(NReachProb*NDecayProb));
	    FillHisto("CouplingScan/hWeightCoupling", fCoupling, TMath::Log10(Weight));	
	    FillHisto("CouplingScan/hEnergyCoupling", fCoupling, TMath::Sqrt(momN/1000.*momN/1000. + MN/1000.*MN/1000.));

	    if (isGood)
	      FillHisto("CouplingScan/hGoodWeightCoupling", fCoupling, Weight);	
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
	    FillHisto("MassScan/hProbMass", round(MN)/1000., TMath::Log10(NReachProb*NDecayProb));
	    FillHisto("MassScan/hWeightMass", round(MN)/1000., TMath::Log10(Weight));
	    FillHisto("MassScan/hEnergyMass", round(MN)/1000., TMath::Sqrt(momN/1000.*momN/1000. + round(MN)/1000.*round(MN)/1000.));

	    if (isGood)
	      FillHisto("MassScan/hGoodWeightMass", round(MN)/1000., Weight);

	    for (Int_t j = 0; j < evt->GetNKineParts(); j++) {
	      KinePart *p = evt->GetKinePart(j);
	      if (p->GetParentID() == -1 && p->GetPDGcode() == 999) {
		if (p->GetParticleName() == "D+->e+N" || p->GetParticleName() == "D-->e-N" || p->GetParticleName() == "D+->mu+N" || p->GetParticleName() == "D-->mu-N" || p->GetParticleName() == "DS+->e+N" || p->GetParticleName() == "DS-->e-N"|| p->GetParticleName() == "DS+->mu+N" || p->GetParticleName() == "DS-->mu-N") {
		  if (p->GetProdPos().Z() <= 400. && p->GetProdPos().Z() >= -400.) {
		    FillHisto("MassScan/hEnergyMassTom", round(MN)/1000., TMath::Sqrt(momN/1000.*momN/1000. + round(MN)/1000.*round(MN)/1000.));	
		  }
		}
	      }
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
	if (WeightsGaia[i]["IsGood"] == 1) {
	  goodWeightGaia = 1.;
	}
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
	  WeightGaia = 1.;

	  if (round(Weights[kineCounter]["Mass"])/1000. == fMassForSingleValue) {
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
	    FillHisto("SingleValue/hDProd", p->GetMomAtCheckPoint(8).X(), p->GetMomAtCheckPoint(8).Y());
	    FillHisto("SingleValue/hXYProd", p->GetProdPos().X(), p->GetProdPos().Y());
	    
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
	  }
	    
	  // Tommaso comparison

	  if (round(TMath::Sqrt(p->GetInitial4Momentum().E()*p->GetInitial4Momentum().E() - p->GetInitial4Momentum().P()*p->GetInitial4Momentum().P()))/1000. == 0.8) {
	    if (p->GetProdPos().Z() <= 400. && p->GetProdPos().Z() >= -400.) {
	      FillHisto("SingleValue/hXYDecay", p->GetEndPos().X()/1000., p->GetEndPos().Y()/1000., WeightTom);
	      if (p->GetParticleName().Contains("0") || ((p->GetParticleName().Contains("tau") && (p->GetParticleName().Contains("e") || p->GetParticleName().Contains("mu"))))) {
		FillHisto("SingleValue/hXYDecay3Body", p->GetEndPos().X()/1000., p->GetEndPos().Y()/1000., WeightTom);
	      }
	      else {
		FillHisto("SingleValue/hXYDecay2Body", p->GetEndPos().X()/1000., p->GetEndPos().Y()/1000., WeightTom);
	      }
	      if (p->GetParticleName() == "D+->e+N" || p->GetParticleName() == "D-->e-N" || p->GetParticleName() == "D+->mu+N" || p->GetParticleName() == "D-->mu-N" || p->GetParticleName() == "DS+->e+N" || p->GetParticleName() == "DS-->e-N"|| p->GetParticleName() == "DS+->mu+N" || p->GetParticleName() == "DS-->mu-N") {
		FillHisto("SingleValue/hDThetaMom", p->GetMomAtCheckPoint(1).Y()/1000., p->GetPosAtCheckPoint(1).x(), WeightTom);
		FillHisto("SingleValue/hXYDecayTom", p->GetEndPos().X()/1000., p->GetEndPos().Y()/1000., WeightTom);
		FillHisto("SingleValue/hTransverse", TMath::Sqrt(p->GetEndPos().X()/1000.*p->GetEndPos().X()/1000. + p->GetEndPos().Y()/1000.*p->GetEndPos().Y()/1000.));
	      }
	    }
	  }
		  
	  // Toy-MC comparison: all

	  if (round(TMath::Sqrt(p->GetInitial4Momentum().E()*p->GetInitial4Momentum().E() - p->GetInitial4Momentum().P()*p->GetInitial4Momentum().P()))/1000. == fMassForSingleValue && p->GetProdPos().Z() <= 400. && p->GetProdPos().Z() >= -400.) {	  
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

void HeavyNeutrinoScanNewEl::EndOfJobUser() {
 
  if (fReadingData) {

    // X axis
 
    fHisto.GetTH1("SingleValue/hZMotherProd")->GetXaxis()->SetTitle("Z coordinate [m]");
    fHisto.GetTH1("SingleValue/hZDProd")->GetXaxis()->SetTitle("Z coordinate [m]");
    fHisto.GetTH1("SingleValue/hZTauProd")->GetXaxis()->SetTitle("Z coordinate [m]");
    fHisto.GetTH1("SingleValue/hZDDecay")->GetXaxis()->SetTitle("Z coordinate [m]");
    fHisto.GetTH1("SingleValue/hDTheta")->GetXaxis()->SetTitle("Polar angle [rad]");
    fHisto.GetTH1("SingleValue/hDLambda")->GetXaxis()->SetTitle("Decay length [mm]");
    fHisto.GetTH1("SingleValue/hDPath")->GetXaxis()->SetTitle("Z [mm]");
    fHisto.GetTH1("SingleValue/hDMom")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("SingleValue/hZHNLDecay")->GetXaxis()->SetTitle("Z coordinate [m]");
    fHisto.GetTH1("SingleValue/hHNLGamma")->GetXaxis()->SetTitle("Lorentz factor");
    fHisto.GetTH1("SingleValue/hHNLTheta")->GetXaxis()->SetTitle("Polar angle [rad]");
    fHisto.GetTH1("SingleValue/hHNLMom")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH1("SingleValue/hTransverse")->GetXaxis()->SetTitle("Distance [m]");
    fHisto.GetTH2("SingleValue/hDThetaMom")->GetXaxis()->SetTitle("P [GeV/c]");
    fHisto.GetTH2("SingleValue/hXYDecay")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("SingleValue/hXYProd")->GetXaxis()->SetTitle("X [mm]");
    fHisto.GetTH2("SingleValue/hXYDecay2Body")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("SingleValue/hXYDecay3Body")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("SingleValue/hXYDecayTom")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("SingleValue/hXYMuon")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("SingleValue/hXYPion")->GetXaxis()->SetTitle("X [m]");
    fHisto.GetTH2("SingleValue/hDProd")->GetXaxis()->SetTitle("X [mm]");
    fHisto.GetTH2("SingleValue/hBeamvsTarTrue")->GetXaxis()->SetTitle("Impact parameter [m]");
    fHisto.GetTH2("SingleValue/hBeamvsTAXTrue")->GetXaxis()->SetTitle("Impact parameter [m]");
    fHisto.GetTH2("CouplingScan/hReachCoupling")->GetXaxis()->SetTitle("Log(U^{2})");
    fHisto.GetTH2("CouplingScan/hDecayCoupling")->GetXaxis()->SetTitle("Log(U^{2})");
    fHisto.GetTH2("CouplingScan/hProbCoupling")->GetXaxis()->SetTitle("Log(U^{2})");
    fHisto.GetTH2("CouplingScan/hWeightCoupling")->GetXaxis()->SetTitle("Log(U^{2})");
    fHisto.GetTH2("CouplingScan/hGoodWeightCoupling")->GetXaxis()->SetTitle("Log(U^{2})");
    fHisto.GetTH2("CouplingScan/hEnergyCoupling")->GetXaxis()->SetTitle("Log(U^{2})");
    fHisto.GetTH2("MassScan/hReachMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("MassScan/hDecayMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("MassScan/hProbMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("MassScan/hWeightMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("MassScan/hGoodWeightMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("MassScan/hEnergyMass")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("MassScan/hEnergyMassTom")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
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
    fHisto.GetTH2("SingleValue/hXYDecay2Body")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("SingleValue/hXYDecay3Body")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("SingleValue/hXYDecayTom")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("SingleValue/hXYMuon")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("SingleValue/hXYPion")->GetYaxis()->SetTitle("Y [m]");
    fHisto.GetTH2("SingleValue/hDProd")->GetYaxis()->SetTitle("Y [mm]");
    fHisto.GetTH2("SingleValue/hXYProd")->GetYaxis()->SetTitle("Y [mm]");
    fHisto.GetTH2("SingleValue/hBeamvsTarTrue")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
    fHisto.GetTH2("SingleValue/hBeamvsTAXTrue")->GetYaxis()->SetTitle("Vertex-beam axis distance [m]");
    fHisto.GetTH2("CouplingScan/hReachCoupling")->GetYaxis()->SetTitle("Reach probability");
    fHisto.GetTH2("CouplingScan/hDecayCoupling")->GetYaxis()->SetTitle("Decay probability");
    fHisto.GetTH2("CouplingScan/hProbCoupling")->GetYaxis()->SetTitle("Log(reach and decay probability)");
    fHisto.GetTH2("CouplingScan/hWeightCoupling")->GetYaxis()->SetTitle("Weight");
    fHisto.GetTH2("CouplingScan/hGoodWeightCoupling")->GetYaxis()->SetTitle("Weight");
    fHisto.GetTH2("CouplingScan/hEnergyCoupling")->GetYaxis()->SetTitle("N energy [GeV]");
    fHisto.GetTH2("MassScan/hReachMass")->GetYaxis()->SetTitle("Reach probability");
    fHisto.GetTH2("MassScan/hDecayMass")->GetYaxis()->SetTitle("Decay probability");
    fHisto.GetTH2("MassScan/hProbMass")->GetYaxis()->SetTitle("Log(each and decay probability)");
    fHisto.GetTH2("MassScan/hWeightMass")->GetYaxis()->SetTitle("Weight");
    fHisto.GetTH2("MassScan/hGoodWeightMass")->GetYaxis()->SetTitle("Weight");
    fHisto.GetTH2("MassScan/hEnergyMass")->GetYaxis()->SetTitle("N energy [GeV]");
    fHisto.GetTH2("MassScan/hEnergyMassTom")->GetYaxis()->SetTitle("N energy [GeV]");
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

    // Z axis

    fHisto.GetTH2("SingleValue/hDThetaMom")->GetZaxis()->SetTitle("Normalized to POT, 0.5 GeV/c and 10 mrad");
    fHisto.GetTH2("SingleValue/hXYDecay")->GetZaxis()->SetTitle("Normalized to POT and 4 cm^{2}");
    fHisto.GetTH2("SingleValue/hBeamvsTarTrue")->GetZaxis()->SetTitle("Normalized to POT and 1 cm^{2}");
    fHisto.GetTH2("SingleValue/hBeamvsTAXTrue")->GetZaxis()->SetTitle("Normalized to POT and 1 cm^{2}");

    Double_t Coupling = 0.;
    Double_t MN = 0.;

    for (auto it = fMasses.begin(); it != fMasses.end(); it++) {
      MN = it->first;
      for (auto it1 = fCouplings.begin(); it1 != fCouplings.end(); it1++) {
	Coupling = it1->first;
	fTMass = MN;
	fTCoupling = Coupling;
	fTNEvents = fNEvents[MN][Coupling];
	fTSumGood = fSumGood[MN][Coupling];
	FillTree("Scan");
      }
    }
  }
  else {

    // Acceptance and yield computation: momentum
    /*    
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomSel"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sG", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sR", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomReg"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sR", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sA", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomFV"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sA", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sN", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorYieldMom"))->Divide((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sG", true), (TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sN", true), "cl=0.683 mode");
    */

    Normalize(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sG", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sR", true)), "Selection acceptance vs N momentum", "ErrorAccMomSel");    
    Normalize(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sR", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sA", true)), "Regeneration acceptance vs N momentum", "ErrorAccMomReg");
    Normalize(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sA", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sN", true)), "FV acceptance vs N momentum", "ErrorAccMomFV");  
    Normalize(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sG", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "SingleValue/sN", true)), "Yield per POT vs N momentum", "ErrorYieldMom");

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
    SumGraphs(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCoupling")), static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTarget")), static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("CouplingScan/ErrorYieldCouplingTAX")));

    // Acceptance and yield computation: mass

    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassSelTarget"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mG1", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR1", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassRegTarget"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR1", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA1", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassFVTarget"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA1", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mN1", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTarget"))->Divide((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mG1", true), (TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mN1", true), "cl=0.683 mode");

    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassSelTAX"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mG2", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR2", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassRegTAX"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR2", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA2", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassFVTAX"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA2", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mN2", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTAX"))->Divide((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mG2", true), (TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mN2", true), "cl=0.683 mode");

    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassSel"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mG", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassReg"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mR", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA", true)), "cl=0.683 mode");
    static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorAccMassFV"))->Divide(static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mA", true)), static_cast<TH1D*>((TH1D*)RequestHistogram(fAnalyzerName, "MassScan/mN", true)), "cl=0.683 mode");
    SumGraphs(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMass")), static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTarget")), static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("MassScan/ErrorYieldMassTAX")));
    
    // Retrieve file info for exclusion plot
    
    Double_t Coupling = 0.;
    Double_t MN = 0.;
    Double_t POT = 1.E18;
    Int_t couplingCounter = 0;
    Int_t massCounter = 0;
    string line;

    TTree *tree = (TTree*)GetCurrentFile()->Get("Scan");

    tree->SetBranchAddress("Mass", &fTMass);
    tree->SetBranchAddress("Coupling", &fTCoupling);
    tree->SetBranchAddress("NEvents", &fTNEvents);
    tree->SetBranchAddress("SumGood", &fTSumGood);

    for (Int_t i = 0; i < tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      fNEvents[fTMass][fTCoupling] += fTNEvents;
      fSumGood[fTMass][fTCoupling] += fTSumGood;
    }

    // Lifetime and exclusion plots

    for (Int_t i = fMassStart*1000.; i <= fMassStop*1000.; i += fMassStep*1000) {
      MN = i;
      for (Int_t j  = fCouplingStart*10; j <= fCouplingStop*10; j += fCouplingStep*10) {
	Coupling = j/10.;

	if (MN/1000. == fMassForSingleValue) {
	  fUSquared = TMath::Power(10., -10.+couplingCounter/10.);
	  fHisto.GetTGraph("CouplingScan/GammaTotCoupling")->SetPoint(couplingCounter, Coupling, GammaTot(MN, false));
	  fHisto.GetTGraph("CouplingScan/TauCoupling")->SetPoint(couplingCounter, Coupling, tauN(MN, false));
	  couplingCounter++;
	}

	if (fSumGood[MN][Coupling]/fNEvents[MN][Coupling] != 0.0/0.0 && fSumGood[MN][Coupling]/fNEvents[MN][Coupling] != 1.0/0.0)
	  fYield[MN][Coupling] = fSumGood[MN][Coupling]/fNEvents[MN][Coupling];
	FillHisto("TotalScan/hExclusion", MN/1000., Coupling, fYield[MN][Coupling]*POT);
      }
    
      Coupling = fCouplingForSingleValue;
      fUSquared = TMath::Power(10., fCouplingForSingleValue);
      fHisto.GetTGraph("MassScan/GammaTotMass")->SetPoint(massCounter, MN/1000., GammaTot(MN, false));
      fHisto.GetTGraph("MassScan/TauMass")->SetPoint(massCounter, MN/1000., tauN(MN, false));
      massCounter++;
    }

    EvaluateUL(fHisto.GetTH2("TotalScan/hExclusion"), fHisto.GetTGraph("TotalScan/Contours"));

    Int_t N = fHisto.GetTGraph("TotalScan/Contours")->GetN();
    Double_t X[N], Y[N];

    if (fInitialUeSquaredRatio != 0. && fInitialUmuSquaredRatio != 0. && fInitialUtauSquaredRatio != 0.) {
      for (Int_t i = 0; i < N; i++)
	fHisto.GetTGraph("TotalScan/Contours")->GetPoint(i, X[i], Y[i]);
      
      fHisto.GetTGraph("TotalScan/Contours")->SetPoint(N,   1.96, Y[N-1]-0.2);
      fHisto.GetTGraph("TotalScan/Contours")->SetPoint(N+1, 1.96, Y[N-1]-0.4);
      fHisto.GetTGraph("TotalScan/Contours")->SetPoint(N+2, 1.96, Y[N-1]-0.8);
      fHisto.GetTGraph("TotalScan/Contours")->SetPoint(N+3, 1.96, Y[N-1]-1.1);
      fHisto.GetTGraph("TotalScan/Contours")->SetPoint(N+4, 1.96, Y[N-1]-1.4);
      fHisto.GetTGraph("TotalScan/Contours")->SetPoint(N+5, 1.96, Y[N-1]-1.9);
      fHisto.GetTGraph("TotalScan/Contours")->SetPoint(N+6, 1.96, Y[N-1]-2.0);
      fHisto.GetTGraph("TotalScan/Contours")->SetPoint(N+7, 1.96, Y[0]+0.1);
      fHisto.GetTGraph("TotalScan/Contours")->SetPoint(N+8, 1.96, Y[0]);
    }

    // Mean probability vs mass and coupling

    for (Int_t i = 1; i < static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "CouplingScan/ProbCoupling", true))->GetNbinsX()+1; i++) {
      TH1D *h = static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "CouplingScan/ProbCoupling", true))->ProjectionY("", i, i+1, "");
      fHisto.GetTGraph("CouplingScan/MeanCoupling")->SetPoint(i-1, static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "CouplingScan/ProbCoupling", true))->GetXaxis()->GetBinCenter(i), h->GetMean());
    }

    for (Int_t i = 1; i < static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "MassScan/ProbMass", true))->GetNbinsX()+1; i++) {      
      TH1D *h = static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "MassScan/ProbMass", true))->ProjectionY("", i, i+1, "");           
      fHisto.GetTGraph("MassScan/MeanMass")->SetPoint(i-1, static_cast<TH2D*>((TH2D*)RequestHistogram(fAnalyzerName, "MassScan/ProbMass", true))->GetXaxis()->GetBinCenter(i), h->GetMean());
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
    /*
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomSel")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomReg")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorAccMomFV")), "Log(U^{2})", "Acceptance", 2);
    CosmeticsGraph(static_cast<TGraphAsymmErrors*>(fHisto.GetTGraph("SingleValue/ErrorYieldMom")), "Log(U^{2})", "Yield per POT", 2);
    */

    CosmeticsHisto(static_cast<TH1D*>(fHisto.GetTH1("SingleValue/ErrorAccMomSel")), "N momentum [GeV/c]", "Acceptance", 2);                                               
    CosmeticsHisto(static_cast<TH1D*>(fHisto.GetTH1("SingleValue/ErrorAccMomReg")), "N momentum [GeV/c]", "Acceptance", 2);                                                                   
    CosmeticsHisto(static_cast<TH1D*>(fHisto.GetTH1("SingleValue/ErrorAccMomFV")), "N momentum [GeV/c]", "Acceptance", 2);   
    CosmeticsHisto(static_cast<TH1D*>(fHisto.GetTH1("SingleValue/ErrorYieldMom")), "N momentum [GeV/c]", "Yield per POT", 2);

    fHisto.GetTH2("TotalScan/hExclusion")->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    fHisto.GetTH2("TotalScan/hExclusion")->GetYaxis()->SetTitle("Log(U^{2})");

    TTree* scanTree = static_cast<TTree*>(GetCurrentFile()->Get("Scan"))->CloneTree();
    TDirectory *d = gDirectory;
    gFile->cd();
    scanTree->Write();
    gDirectory = d;
  }

  SaveAllPlots();
  
  return;
}

void HeavyNeutrinoScanNewEl::CosmeticsGraph(TGraph* g, const char* x, const char* y, Int_t col) {

  g->GetXaxis()->SetTitle(x);
  g->GetYaxis()->SetTitle(y);
  g->SetLineColor(col);
  g->SetLineWidth(3);
  g->Draw("AC");
  
  return;
}

void HeavyNeutrinoScanNewEl::CosmeticsGraph(TGraphAsymmErrors* g, const char* x, const char* y, Int_t col) {

  g->GetXaxis()->SetTitle(x);
  g->GetYaxis()->SetTitle(y);
  g->SetLineColor(1);
  g->SetLineWidth(1);
  g->SetMarkerStyle(8);
  g->SetMarkerColor(col);
  g->Draw("AP");
  
  return;
}

void HeavyNeutrinoScanNewEl::CosmeticsHisto(TH1* h, const char* x, const char* y, Int_t col) {

  h->GetXaxis()->SetTitle(x);
  h->GetYaxis()->SetTitle(y);
  h->SetLineColor(1);
  h->SetLineWidth(1);
  h->SetMarkerStyle(8);
  h->SetMarkerColor(col);
  h->Draw("E");

  return;
}

void HeavyNeutrinoScanNewEl::EvaluateUL(TH2* h, TGraph* gr) {

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
    Double_t sigLevel = -TMath::Log(1-CLs[i]); // maximum number of events, given the CL, for Nobs = 0 and Bkg = 0                     

    contours[i] = sigLevel;
  }

  h->SetContour(nContours, contours);

  std::vector<TGraph*> contoursHisto = ExtractContours(h);

  if (contoursHisto.size() == 0) {
    cout<<"Graf for UL is empty"<<endl;
    return;
  }

  Int_t level = 1;
  //Int_t level = 0; // contour index to work on
  Int_t nn = contoursHisto[level]->GetN();

  TCanvas *c = new TCanvas();
  contoursHisto[level]->Draw();
  c->Print("cont.pdf");
  
  Double_t* xx = new Double_t[nn];
  Double_t* yy = new Double_t[nn];

  xx = contoursHisto[level]->GetX();
  yy = contoursHisto[level]->GetY();

  for (Int_t i = 0; i < nn; i++) {
    gr->SetPoint(i, xx[i], yy[i]);
  }
  
  gr->SetLineColor(kRed);
  gr->SetLineWidth(3);
  gr->Draw("AP*");
  
  return;
}
  
std::vector<TGraph*> HeavyNeutrinoScanNewEl::ExtractContours(TH2* h) {
  
  // Draw contours as filled regions and save points 
  
  TCanvas *c = new TCanvas();
  h->Draw("cont z list");
  c->Print("contz.pdf");

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

  return graf;
}

void HeavyNeutrinoScanNewEl::SumGraphs(TGraphAsymmErrors* res, TGraphAsymmErrors* g1, TGraphAsymmErrors* g2) {

  for (Int_t i = 0; i < g1->GetN(); i++) {
    Double_t x1, y1, x2, y2;
    g1->GetPoint(i, x1, y1);
    g2->GetPoint(i, x2, y2);
    res->SetPoint(i, x1, y1+y2);
    res->SetPointError(i, g1->GetErrorXlow(i), g1->GetErrorXhigh(i), TMath::Sqrt(g1->GetErrorYlow(i)*g1->GetErrorYlow(i)+g2->GetErrorYlow(i)*g2->GetErrorYlow(i)), TMath::Sqrt(g1->GetErrorYhigh(i)*g1->GetErrorYhigh(i)+g2->GetErrorYhigh(i)*g2->GetErrorYhigh(i)));
  }

  return;
}

// Function to compute yield vs mom and renormalize it to the total yield for a certain coupling value 

void HeavyNeutrinoScanNewEl::Normalize(TH1D* h1, TH1D* h2, TString title, TString name) {

  Double_t I = h2->Integral();
  TH1D* res = (TH1D*)h1->Clone(name);
  res->Scale(1./I);
  res->SetTitle(title);
  BookHisto("SingleValue/" + name, res);

  return;
}
