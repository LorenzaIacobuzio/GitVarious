// ---------------------------------------------------------------                                      
//                                                                                                      
// History:                                                                                             
//                                                                                                      
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018                               
//                                                                                                      
// ---------------------------------------------------------------                                      
/// \class HeavyNeutrinoScan
/// \Brief                                                                                              
/// Produce expected yield per POT for HNL MC samples, for single values of HNL mass and coupling, for 
/// a scan on the coupling or on the mass, or for both
/// \EndBrief                                                                                           
/// \Detailed                                                                                           
/// After retrieving all HNL weights using the tool called HNLWeight, 
/// and checking if an event passed the selection implemented in the analyzer HeavyNeutrino,  
/// acceptance and yield per POT are computed. 
/// Different sets of plots are produced in different cases:
/// if the sample contains only one mass HNLs and the coupling has been set to a constant, acceptance
/// and yield are computed as single values; the yield is also computed as a function of the HNL momentum;
/// if a scan on the coupling is enabled, plots are produced as a function of the coupling;
/// if the MC sample contains several masses, plots are produced as a function of the HNL mass;
/// finally, if a scan on the coupling is enabled on a MC sample with several masses, an expected
/// exclusion plot is produced as a function of HNL mass and coupling.                      
/// Thus, the analyzer is able to run either on MC samples of just one HNL mass or on samples
/// containing HNLs of different masses.                                                                
/// This analyzer makes use of two ToolsLib, called HNLFunctions and HNLWeight.  
/// The values of the ratios between specific-flavour couplings can be either set as external           
/// parameters from command line or taken as default values.                                            
/// For example, if the user sets USquared = 1.E-10, UeSquaredRatio = 5., UmuSquaredRatio = 1., 
/// UtauSquaredRatio = 3.5, the specific-flavour coupling values will be: UeSquared = 5.25E-11, 
/// UmuSquared = 1.05E-11, UtauSquared = 3.68E-11.   
/// Default values are: USquared = 1.E-6, UeSquaredRatio = 1., UmuSquaredRatio = 16.,         
/// UtauSquaredRatio = 3.8.
/// The values of the initial and final coupling and the scan step can be either set as external      
/// parameters from command line or taken as default values.                                           
/// For example, if the user assigns -3. to the starting coupling, -2. to the final one                
/// and 0.5 to the step, the scan will be performed for Log(U2) = [-3., -2.5, -2.], this meaning       
/// U2 = [1.E-3., 1.E-2.5, 1.E-2.].                                                                    
/// Default values are: CouplingStart = -10., CouplingStop = 0., CouplingStep = 0.1.
/// A boolean to enable/disable the coupling scan can be also set as an external parameter.
/// The coupling scan is enabled by default.
/// This analyzer produces .txt files needed for computing error bars for plots generated during
/// the first step, and it can be run in --histo mode to produce these plots:
/// ./myExecutable -i outFile.root -o outFile1.root --histo
/// The files are deleted at the end of the second step.
/// The values of the initial and final momentum and its step can be either set as external         
/// parameters from command line or taken as default values. These are used to produce plots 
/// of the acceptance and yield per POT as a function of the HNL momentum.
/// Default values are: MomStart = 0., MomStop = 200., MomStep = 10.
/// Other parameters to be set from command line are:
/// mass value to produce plots at a fixed mass (coupling scan and single value plots); default value
/// is MassForSingleValue = 1000.;
/// number of events in each chunk the sample is split into, to compute error bars in --histo mode;
/// default value is SplitStep = 10.
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
  AddParam("CouplingStep", &fCouplingStep, 5.);
  AddParam("EnableCouplingScan", &fEnableCouplingScan, true);
  AddParam("InitialFV", &fInitialFV, 102500.);
  AddParam("LFV", &fLFV, 77500.);
  AddParam("MomStop", &fMomStop, 200.);
  AddParam("MomStart", &fMomStart, 0.);
  AddParam("MomStep", &fMomStep, 10.);
  AddParam("MassForSingleValue", &fMassForSingleValue, 1000.);
  AddParam("SplitStep", &fSplitStep, 10);

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
    
    // Scan histos
    
    fhReachCoupling       = nullptr;
    fhDecayCoupling       = nullptr;
    fhProbCoupling        = nullptr;
    fhWeightCoupling      = nullptr;
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
    
    BookHisto("SingleValue/hZMotherProd",  new TH1D("ZMotherProd", "Z of N mother production point", 20000., -250., 33000.));
    BookHisto("SingleValue/hZDProd",       new TH1D("ZDProd",      "Z of D meson production point",  20000., -250., 33000.));
    BookHisto("SingleValue/hZTauProd",     new TH1D("ZTauProd",    "Z of tau production point",      20000., -250., 33000.));
    BookHisto("SingleValue/hZDDecay",      new TH1D("ZDDecay",     "Z of N production point",        20000., -250., 33000.));
    BookHisto("SingleValue/hDTheta",       new TH1D("DTheta",  "N mother theta",        100,  0., 0.3));
    BookHisto("SingleValue/hDLambda",      new TH1D("DLambda", "N mother decay length", 100, -1., 40.));
    BookHisto("SingleValue/hDPath",        new TH1D("DPath",   "N mother path in Z",    100, -1., 50.));
    BookHisto("SingleValue/hDMom",         new TH1D("DMom",    "N mother momentum",     100, -1., 170.));
    BookHisto("SingleValue/hZHNLDecay",    new TH1D("ZHNLDecay", "Z of N decay point", 100., 90., 190.));
    BookHisto("SingleValue/hHNLGamma",     new TH1D("HNLGamma", "N Lorentz gamma", 50., 0., 170.));
    BookHisto("SingleValue/hHNLDecayProb", new TH1D("HNLDecayProb", "N decay probability", 100., 0., 0.0065));
    BookHisto("SingleValue/hHNLReachProb", new TH1D("HNLReachProb", "N probability of reaching FV", 100., 0.99, 1.001));
    BookHisto("SingleValue/hHNLTheta",     new TH1D("HNLTheta", "N theta", 100., 0., 0.5));
    BookHisto("SingleValue/hHNLMom",       new TH1D("HNLMom", "N momentum", 100., -0.5, 200.));
    BookHisto("SingleValue/hWeight",       new TH1D("Weight", "N Weight", 1000, 1.E-40, 1.E-5));
    BookHisto("SingleValue/hCoupling",     new TH1D("Coupling", "Coupling", 100, -10., 0.));
    BookHisto("SingleValue/hMass",         new TH1D("Mass", "Mass", 10, 0.3, 1.2));
    BookHisto("SingleValue/hAcc",          new TH1D("Acc",       "Acceptance for one value of N mass and coupling",                             1000, 1.E-30, 1.));
    BookHisto("SingleValue/hAccTarget",    new TH1D("AccTarget", "Acceptance for one value of N mass and coupling, for target-produced events", 1000, 1.E-30, 1.));
    BookHisto("SingleValue/hAccTAX",       new TH1D("AccTAX",    "Acceptance for one value of N mass and coupling, for TAX-produced events",    1000, 1.E-30, 1.));
    BookHisto("SingleValue/hYield",        new TH1D("Yield",       "Yield per POT for one value of N mass and coupling",                                1000, 1.E-50, 1.E-10));
    BookHisto("SingleValue/hYieldTarget",  new TH1D("YieldTarget", "Yield per POT for one value of N mass and coupling, for target-produced events",    1000, 1.E-50, 1.E-10));
    BookHisto("SingleValue/hYieldTAX",     new TH1D("YieldTAX",    "Yield per POT for one value of N mass and coupling, for TAX-produced events",       1000, 1.E-50, 1.E-10));

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
    BookHisto("CouplingScan/hWeightCoupling", new TH2D("WeightCoupling", "N weight vs coupling",                                         fN, fCouplingStart, fCouplingStop, 1000, 1.E-25, 1.E-5));
    
    fgGammaTotCoupling = new TGraphErrors();
    fgGammaTotCoupling->SetNameTitle("CouplingScan/GammaTotCoupling", "N total decay width vs coupling");
    BookHisto(fgGammaTotCoupling);
    
    fgTauCoupling = new TGraphErrors();
    fgTauCoupling->SetNameTitle("CouplingScan/TauCoupling", "N lifetime vs coupling");
    BookHisto(fgTauCoupling);
    
    fgAccCoupling = new TGraphErrors();
    fgAccCoupling->SetNameTitle("CouplingScan/AccCoupling", "Acceptance vs coupling");
    BookHisto(fgAccCoupling);
    
    fgAccCouplingTarget = new TGraphErrors();
    fgAccCouplingTarget->SetNameTitle("CouplingScan/AccCouplingTarget", "Acceptance per POT vs coupling, for target-produced events");
    BookHisto(fgAccCouplingTarget);
    
    fgAccCouplingTAX = new TGraphErrors();
    fgAccCouplingTAX->SetNameTitle("CouplingScan/AccCouplingTAX", "Acceptance per POT vs coupling, for TAX-produced events");
    BookHisto(fgAccCouplingTAX);
    
    fgYieldCoupling = new TGraphErrors();
    fgYieldCoupling->SetNameTitle("CouplingScan/YieldCoupling", "Yield per POT vs coupling");
    BookHisto(fgYieldCoupling);
    
    fgYieldCouplingTarget = new TGraphErrors();
    fgYieldCouplingTarget->SetNameTitle("CouplingScan/YieldCouplingTarget", "Yield per POT vs coupling, for target-produced events");
    BookHisto(fgYieldCouplingTarget);
    
    fgYieldCouplingTAX = new TGraphErrors();
    fgYieldCouplingTAX->SetNameTitle("CouplingScan/YieldCouplingTAX", "Yield per POT vs coupling, for TAX-produced events");
    BookHisto(fgYieldCouplingTAX);
    
    // Mass scan
    
    BookHisto("MassScan/hReachMass",  new TH2D("ReachMass",  "Probability of N reaching the FV vs N mass",                 10, 0.25, 1.25, 1000, -0.1, 1.1));
    BookHisto("MassScan/hDecayMass",  new TH2D("DecayMass",  "Probability of N decaying in the FV vs N mass",              10, 0.25, 1.25, 1000, -0.1, 1.1));
    BookHisto("MassScan/hProbMass",   new TH2D("ProbMass",   "Probability of N reaching and decaying in the FV vs N mass", 10, 0.25, 1.25, 1000, -0.1, 1.1));
    BookHisto("MassScan/hWeightMass", new TH2D("WeightMass", "N weight vs N mass",                                         10, 0.25, 1.25, 1000, 1.E-25, 1.E-5));
    
    fgGammaTotMass = new TGraphErrors();
    fgGammaTotMass->SetNameTitle("MassScan/GammaTotMass", "N total decay width vs N mass");
    BookHisto(fgGammaTotMass);
    
    fgTauMass = new TGraphErrors();
    fgTauMass->SetNameTitle("MassScan/TauMass", "N lifetime vs N mass");
    BookHisto(fgTauMass);
    
    fgAccMass = new TGraphErrors();
    fgAccMass->SetNameTitle("MassScan/AccMass", "Acceptance vs N mass");
    BookHisto(fgAccMass);
    
    fgAccMassTarget = new TGraphErrors();
    fgAccMassTarget->SetNameTitle("MassScan/AccMassTarget", "Acceptance per POT vs N mass, for target-produced events");
    BookHisto(fgAccMassTarget);
    
    fgAccMassTAX = new TGraphErrors();
    fgAccMassTAX->SetNameTitle("MassScan/AccMassTAX", "Acceptance per POT vs N mass, for TAX-produced events");
    BookHisto(fgAccMassTAX);
    
    fgYieldMass = new TGraphErrors();
    fgYieldMass->SetNameTitle("MassScan/YieldMass", "Yield per POT vs N mass");
    BookHisto(fgYieldMass);
    
    fgYieldMassTarget = new TGraphErrors();
    fgYieldMassTarget->SetNameTitle("MassScan/YieldMassTarget", "Yield per POT vs N mass, for target-produced events");
    BookHisto(fgYieldMassTarget);
    
    fgYieldMassTAX = new TGraphErrors();
    fgYieldMassTAX->SetNameTitle("MassScan/YieldMassTAX", "Yield per POT vs N mass, for TAX-produced events");
    BookHisto(fgYieldMassTAX);
    
    // Total scan
    
    fgExclusion = new TGraphErrors();
    fgExclusion->SetNameTitle("TotalScan/Exclusion", "Sensitivity as a function of coupling and N mass");
    BookHisto(fgExclusion);
  }
  else {

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

  if (fEnableCouplingScan == false) {
    fCouplingStart = TMath::Log10(fUSquared);
    fCouplingStop = fCouplingStart;
  }

  //Bool_t IsHNLGood = *(Bool_t*)GetOutput("HeavyNeutrino.Output");
  //REMOVE
  Bool_t IsHNLGood = true;

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

      std::vector<std::map<std::string, Double_t>> Weights = ComputeWeight(evt, fUSquared, fUeSquaredRatio, fUmuSquaredRatio, fUtauSquaredRatio, fLInitialFV, fLFV);

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
	  if (Weights[i]["ProdProb"] == fDBeProdProb) 
	    fSumGoodTarget[round(MN)][fCoupling] += Weight;
	  else if (Weights[i]["ProdProb"] == fDCuProdProb) 
	    fSumGoodTAX[round(MN)][fCoupling] += Weight;
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
	    fErrorFile << round(MN) << "\t" << fCoupling << "\t" << "0." << "\t" << "0." << endl;
	  if (fSumAllTAX[round(MN)][fCoupling] != 0. && fNeventsTAX[round(MN)][fCoupling] != 0)
	    fErrorFileTAX << round(MN) << "\t" << fCoupling << "\t" << fSumGoodTAX[round(MN)][fCoupling]/fSumAllTAX[round(MN)][fCoupling] << "\t" << fSumGoodTAX[round(MN)][fCoupling]/fNeventsTAX[round(MN)][fCoupling] << endl;
	  else
	    fErrorFile << round(MN) << "\t" << fCoupling << "\t" << "0." << "\t" << "0." << endl;
	}

	FillHisto("CouplingScan/hReachCoupling",  fCoupling,       NReachProb);
	FillHisto("CouplingScan/hDecayCoupling",  fCoupling,       NDecayProb);
        FillHisto("CouplingScan/hProbCoupling",   fCoupling,       NReachProb*NDecayProb);
	FillHisto("CouplingScan/hWeightCoupling", fCoupling,       Weight);	
	FillHisto("MassScan/hReachMass",          round(MN)/1000., NReachProb);
	FillHisto("MassScan/hDecayMass",          round(MN)/1000., NDecayProb);
        FillHisto("MassScan/hProbMass",           round(MN)/1000., NReachProb*NDecayProb);
	FillHisto("MassScan/hWeightMass",         round(MN)/1000., Weight);

	if (fEnableCouplingScan == false || (fEnableCouplingScan == true && fCoupling == (fCouplingStart - fCouplingStop)/2. && round(MN) == fMassForSingleValue)) {
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

    fUSquared = TMath::Power(10., (fCouplingStart - fCouplingStop)/2.);
    std::vector<std::map<std::string, Double_t>> Weights = ComputeWeight(evt, fUSquared, fUeSquaredRatio, fUmuSquaredRatio, fUtauSquaredRatio, fLInitialFV, fLFV);
    
    for (UInt_t i = 0; i < Weights.size(); i++) {
      MN =  round(Weights[i]["Mass"]);
      isGood = Weights[i]["IsGood"];
      if (MN == fMassForSingleValue) {
	momN = Weights[i]["Momentum"]/1000.;
	momBin = fMomStep*trunc(momN/fMomStep);
	fMomenta[momBin] = momBin;

	if (fNeventsMom.count(momBin) == 0)
	  fNeventsMom[momBin] = 0;
	fNeventsMom[momBin]++;
        fSumAllMom[momBin] += Weight;

	if (IsHNLGood == true && isGood == true)
	  fSumGoodMom[momBin] += Weight;

        if (i == 0) {
          if (fEvtCounterMom.count(momBin) == 0)
            fEvtCounterMom[momBin] = 0;
          fEvtCounterMom[momBin]++;
        }

        if (fEvtCounterMom[momBin] % fSplitStep == 0 && i == Weights.size() - 1) {
	  if (fSumAllMom[momBin] != 0. && fNeventsMom[momBin] != 0)
	    fErrorFileMom << momBin << "\t" << fSumGoodMom[momBin]/fSumAllMom[momBin] << "\t" << fSumGoodMom[momBin]/fNeventsMom[momBin] << endl;
	  else
	    fErrorFileMom << momBin << "\t" << "0." << "\t" << "0." << endl;
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
	FillHisto(  "SingleValue/hZMotherProd", p->GetPosAtCheckPoint(0).z());

	if (p->GetPosAtCheckPoint(1).z() != 0.)
	  FillHisto("SingleValue/hZDProd",      p->GetPosAtCheckPoint(1).z());
	if (p->GetPosAtCheckPoint(2).z() != 0.)
	  FillHisto("SingleValue/hZTauProd",    p->GetPosAtCheckPoint(2).z());

	FillHisto(  "SingleValue/hZDDecay",     p->GetProdPos().Z());
	FillHisto(  "SingleValue/hDTheta",      p->GetPosAtCheckPoint(1).x());
	FillHisto(  "SingleValue/hDLambda",     p->GetPosAtCheckPoint(1).y());
	FillHisto(  "SingleValue/hDPath",       p->GetMomAtCheckPoint(1).X());
	FillHisto(  "SingleValue/hDMom",        p->GetMomAtCheckPoint(1).Y()/1000.);
	FillHisto(  "SingleValue/hZHNLDecay",   p->GetEndPos().Z()/1000.);
	FillHisto(  "SingleValue/hHNLGamma",    p->GetInitial4Momentum().Gamma());
	FillHisto(  "SingleValue/hHNLTheta",    p->GetMomAtCheckPoint(0).Z());
	FillHisto(  "SingleValue/hHNLMom",      p->GetMomAtCheckPoint(0).T()/1000.);
      }
    }
  }
}

void HeavyNeutrinoScan::EndOfJobUser() {

  if (fReadingData) {
  
    // Retrieve histos

    fhZMotherProd      = (TH1D*)fHisto.GetTH1("SingleValue/hZMotherProd");
    fhZDProd           = (TH1D*)fHisto.GetTH1("SingleValue/hZDProd");
    fhZTauProd         = (TH1D*)fHisto.GetTH1("SingleValue/hZTauProd");
    fhZDDecay          = (TH1D*)fHisto.GetTH1("SingleValue/hZDDecay");
    fhDTheta           = (TH1D*)fHisto.GetTH1("SingleValue/hDTheta");
    fhDLambda          = (TH1D*)fHisto.GetTH1("SingleValue/hDLambda");
    fhDPath            = (TH1D*)fHisto.GetTH1("SingleValue/hDPath");
    fhDMom             = (TH1D*)fHisto.GetTH1("SingleValue/hDMom");
    fhZHNLDecay        = (TH1D*)fHisto.GetTH1("SingleValue/hZHNLDecay");
    fhHNLGamma         = (TH1D*)fHisto.GetTH1("SingleValue/hHNLGamma");
    fhHNLDecayProb     = (TH1D*)fHisto.GetTH1("SingleValue/hHNLDecayProb");
    fhHNLReachProb     = (TH1D*)fHisto.GetTH1("SingleValue/hHNLReachProb");
    fhHNLTheta         = (TH1D*)fHisto.GetTH1("SingleValue/hHNLTheta");
    fhHNLMom           = (TH1D*)fHisto.GetTH1("SingleValue/hHNLMom");
    fhWeight           = (TH1D*)fHisto.GetTH1("SingleValue/hWeight");
    fhCoupling         = (TH1D*)fHisto.GetTH1("SingleValue/hCoupling");
    fhMass             = (TH1D*)fHisto.GetTH1("SingleValue/hMass");
    fhAcc              = (TH1D*)fHisto.GetTH1("SingleValue/hAcc");
    fhAccTarget        = (TH1D*)fHisto.GetTH1("SingleValue/hAccTarget");
    fhAccTAX           = (TH1D*)fHisto.GetTH1("SingleValue/hAccTAX");
    fhYield            = (TH1D*)fHisto.GetTH1("SingleValue/hYield");
    fhYieldTarget      = (TH1D*)fHisto.GetTH1("SingleValue/hYieldTarget");
    fhYieldTAX         = (TH1D*)fHisto.GetTH1("SingleValue/hYieldTAX");
    fhReachCoupling    = (TH2D*)fHisto.GetTH2("CouplingScan/hReachCoupling");
    fhDecayCoupling    = (TH2D*)fHisto.GetTH2("CouplingScan/hDecayCoupling");
    fhProbCoupling     = (TH2D*)fHisto.GetTH2("CouplingScan/hProbCoupling");
    fhWeightCoupling   = (TH2D*)fHisto.GetTH2("CouplingScan/hWeightCoupling");
    fhReachMass        = (TH2D*)fHisto.GetTH2("MassScan/hReachMass");
    fhDecayMass        = (TH2D*)fHisto.GetTH2("MassScan/hDecayMass");
    fhProbMass         = (TH2D*)fHisto.GetTH2("MassScan/hProbMass");
    fhWeightMass       = (TH2D*)fHisto.GetTH2("MassScan/hWeightMass");

    // X axis title

    fhZMotherProd     ->GetXaxis()->SetTitle("Z [mm]");
    fhZDProd          ->GetXaxis()->SetTitle("Z [mm]");
    fhZTauProd        ->GetXaxis()->SetTitle("Z [mm]");
    fhZDDecay         ->GetXaxis()->SetTitle("Z [mm]");
    fhDTheta          ->GetXaxis()->SetTitle("Theta [rad]");
    fhDLambda         ->GetXaxis()->SetTitle("Decay length [mm]");
    fhDPath           ->GetXaxis()->SetTitle("Z [mm]");
    fhDMom            ->GetXaxis()->SetTitle("P [GeV]");
    fhZHNLDecay       ->GetXaxis()->SetTitle("Z [m]");
    fhHNLGamma        ->GetXaxis()->SetTitle("Lorentz gamma");
    fhHNLDecayProb    ->GetXaxis()->SetTitle("Decay probability");
    fhHNLReachProb    ->GetXaxis()->SetTitle("Reach probability");
    fhHNLTheta        ->GetXaxis()->SetTitle("Theta [rad]");
    fhHNLMom          ->GetXaxis()->SetTitle("P [GeV]");
    fhWeight          ->GetXaxis()->SetTitle("Weight");
    fhCoupling        ->GetXaxis()->SetTitle("Coupling");
    fhMass            ->GetXaxis()->SetTitle("Mass [GeV]");
    fhAcc             ->GetXaxis()->SetTitle("Acceptance");
    fhAccTarget       ->GetXaxis()->SetTitle("Acceptance");
    fhAccTAX          ->GetXaxis()->SetTitle("Acceptance");
    fhYield           ->GetXaxis()->SetTitle("Yield per POT");
    fhYieldTarget     ->GetXaxis()->SetTitle("Yield per POT");
    fhYieldTAX        ->GetXaxis()->SetTitle("Yield per POT");
    fhReachCoupling   ->GetXaxis()->SetTitle("Log of coupling");
    fhDecayCoupling   ->GetXaxis()->SetTitle("Log of coupling");
    fhProbCoupling    ->GetXaxis()->SetTitle("Log of coupling");
    fhWeightCoupling  ->GetXaxis()->SetTitle("Log of coupling");
    fhReachMass       ->GetXaxis()->SetTitle("N mass [GeV]");
    fhDecayMass       ->GetXaxis()->SetTitle("N mass [GeV]");
    fhProbMass        ->GetXaxis()->SetTitle("N mass [GeV]");
    fhWeightMass      ->GetXaxis()->SetTitle("N mass [GeV]");

    // Y axis title

    fhReachCoupling   ->GetYaxis()->SetTitle("Reach probability");
    fhDecayCoupling   ->GetYaxis()->SetTitle("Decay probability");
    fhProbCoupling    ->GetYaxis()->SetTitle("Reach and decay probability");
    fhWeightCoupling  ->GetYaxis()->SetTitle("Weight");
    fhReachMass       ->GetYaxis()->SetTitle("Reach probability");
    fhDecayMass       ->GetYaxis()->SetTitle("Decay probability");
    fhProbMass        ->GetYaxis()->SetTitle("Reach and decay probability");
    fhWeightMass      ->GetYaxis()->SetTitle("Weight");

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

	if (fEnableCouplingScan == false || (fEnableCouplingScan == true && Coupling == (fCouplingStart - fCouplingStop)/2. && MN == fMassForSingleValue)) {      

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
    
      Coupling = (fCouplingStart - fCouplingStop)/2.;
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

    CosmeticsGraph(fgGammaTotCoupling,    "Log of coupling",  "Total decay width [MeV]", 2);
    CosmeticsGraph(fgTauCoupling,         "Log of coupling",  "Lifetime [ns]",           2);
    CosmeticsGraph(fgAccCoupling,         "Log of coupling",  "Acceptance",              2);
    CosmeticsGraph(fgAccCouplingTarget,   "Log of coupling",  "Acceptance",              8);
    CosmeticsGraph(fgAccCouplingTAX,      "Log of coupling",  "Acceptance",              9);
    CosmeticsGraph(fgYieldCoupling,       "Log of coupling",  "Yield per POT",           2);
    CosmeticsGraph(fgYieldCouplingTarget, "Log of coupling",  "Yield per POT",           8);
    CosmeticsGraph(fgYieldCouplingTAX,    "Log of coupling",  "Yield per POT",           9);
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
    CosmeticsGraph(fgExclusion,           "N mass [GeV]",     "Log of coupling",         2);
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
  
  for (auto it = errorListAcc.begin(); it != errorListAcc.end(); it++)
    errorResAcc[it->first] = ComputeRMS(it->second)[1];

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

  for (auto it = v.begin(); it != v.end(); it++)
    Sum += *it;

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
    fhReachCoupling       = nullptr;
    fhDecayCoupling       = nullptr;
    fhProbCoupling        = nullptr;
    fhWeightCoupling      = nullptr;
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