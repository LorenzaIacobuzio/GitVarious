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
/// and yield are computed as single values;
/// if a scan on the coupling is enabled on a one mass MC sample, plots are produced as a function
/// of the coupling;
/// if the coupling is fixed on a MC sample with several masses, plots are produced as a function
/// of the HNL mass;
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
  
  if (!GetIsTree()) return;

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

  fN = round((std::abs(fCouplingStop-fCouplingStart))/fCouplingStep);
  fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
  fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
  fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;
  fMassForSingleValue = 1000.;
  errorCounter       = 0;
  errorCounterTarget = 0;
  errorCounterTAX    = 0;
  errorStep          = 10;
  errorFile.open("ErrorBars.txt");
  errorFileTarget.open("ErrorBarsTarget.txt");
  errorFileTAX.open("ErrorBarsTAX.txt");

  // One value histos

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
  fhWeightMass          = nullptr;
  fgAccMass             = nullptr;
  fgAccMassTarget       = nullptr;
  fgAccMassTAX          = nullptr;
  fgYieldMass           = nullptr;
  fgYieldMassTarget     = nullptr;
  fgYieldMassTAX        = nullptr;
  fgGammaTotMass        = nullptr;
  fgTauMass             = nullptr;
  fgExclusion           = nullptr;
}

void HeavyNeutrinoScan::InitHist() {

  // One value of mass and coupling

  BookHisto("SingleValue/hZDProd",       new TH1D("ZDProd", "Z of N mother production point", 20000., -250., 33000.));
  BookHisto("SingleValue/hZDDecay",      new TH1D("ZDDecay", "Z of N mother decay point",     20000., -250., 33000.));
  BookHisto("SingleValue/hDTheta",       new TH1D("DTheta",  "N mother theta",        100,  0., 0.3));
  BookHisto("SingleValue/hDLambda",      new TH1D("DLambda", "N mother decay length", 100, -1., 40.));
  BookHisto("SingleValue/hDPath",        new TH1D("DPath",   "N mother path in Z",    100, -1., 50.));
  BookHisto("SingleValue/hDMom",         new TH1D("DMom",    "N mother momentum",     100, -1., 170.));
  BookHisto("SingleValue/hZHNLDecay",    new TH1D("ZHNLDecay", "Z of HNL decay point", 100., 90., 190.));
  BookHisto("SingleValue/hHNLGamma",     new TH1D("HNLGamma", "Lorentz gamma of HNL", 50., 0., 170.));
  BookHisto("SingleValue/hHNLDecayProb", new TH1D("HNLDecayProb", "HNL decay probability", 100., 0., 0.0065));
  BookHisto("SingleValue/hHNLReachProb", new TH1D("HNLReachProb", "HNL probability of reaching FV", 100., 0.99, 1.001));
  BookHisto("SingleValue/hHNLTheta",     new TH1D("HNLTheta", "HNL theta", 100., 0., 0.5));
  BookHisto("SingleValue/hHNLMom",       new TH1D("HNLMom", "HNL momentum", 100., -0.5, 200.));
  BookHisto("SingleValue/hWeight",       new TH1D("Weight", "Weight", 1000, 1.E-40, 1.E-5));
  BookHisto("SingleValue/hCoupling",     new TH1D("Coupling", "Coupling", 100, -10., 0.));
  BookHisto("SingleValue/hMass",         new TH1D("Mass", "Mass", 10, 0.3, 1.2));
  BookHisto("SingleValue/hAcc",          new TH1D("Acc", "Acceptance for one value of N mass and coupling",                                           1000, 1.E-30, 1.));
  BookHisto("SingleValue/hAccTarget",    new TH1D("AccTarget", "Acceptance per POT for one value of N mass and coupling, for target-produced events", 1000, 1.E-30, 1.));
  BookHisto("SingleValue/hAccTAX",       new TH1D("AccTAX", "Acceptance per POT for one value of N mass and coupling, for TAX-produced events",       1000, 1.E-30, 1.));
  BookHisto("SingleValue/hYield",        new TH1D("Yield", "Yield per POT for one value of N mass and coupling",                                      1000, 1.E-50, 1.E-10));
  BookHisto("SingleValue/hYieldTarget",  new TH1D("YieldTarget", "Yield per POT for one value of N mass and coupling, for target-produced events",    1000, 1.E-50, 1.E-10));
  BookHisto("SingleValue/hYieldTAX",     new TH1D("YieldTAX", "Yield per POT for one value of N mass and coupling, for TAX-produced events",          1000, 1.E-50, 1.E-10));

  // Coupling scan 

  BookHisto("CouplingScan/hReachCoupling",  new TH2D("ReachCoupling", "Probability of N reaching the FV vs coupling",    fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("CouplingScan/hDecayCoupling",  new TH2D("DecayCoupling", "Probability of N decaying in the FV vs coupling", fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("CouplingScan/hWeightCoupling", new TH2D("WeightCoupling", "N weight vs coupling",                           fN, fCouplingStart, fCouplingStop, 1000, 1.E-25, 1.E-5));

  fgGammaTotCoupling = new TGraph();
  fgGammaTotCoupling->SetNameTitle("CouplingScan/GammaTotCoupling", "N total decay width vs coupling");
  BookHisto(fgGammaTotCoupling);

  fgTauCoupling = new TGraph();
  fgTauCoupling->SetNameTitle("CouplingScan/TauCoupling", "N lifetime vs coupling");
  BookHisto(fgTauCoupling);

  fgAccCoupling = new TGraph();
  fgAccCoupling->SetNameTitle("CouplingScan/AccCoupling", "Acceptance vs coupling");
  BookHisto(fgAccCoupling);

  fgAccCouplingTarget = new TGraph();
  fgAccCouplingTarget->SetNameTitle("CouplingScan/AccCouplingTarget", "Acceptance per POT vs coupling, for target-produced events");
  BookHisto(fgAccCouplingTarget);

  fgAccCouplingTAX = new TGraph();
  fgAccCouplingTAX->SetNameTitle("CouplingScan/AccCouplingTAX", "Acceptance per POT vs coupling, for TAX-produced events");
  BookHisto(fgAccCouplingTAX);

  fgYieldCoupling = new TGraph();
  fgYieldCoupling->SetNameTitle("CouplingScan/YieldCoupling", "Yield per POT vs coupling");
  BookHisto(fgYieldCoupling);

  fgYieldCouplingTarget = new TGraph();
  fgYieldCouplingTarget->SetNameTitle("CouplingScan/YieldCouplingTarget", "Yield per POT vs coupling, for target-produced events");
  BookHisto(fgYieldCouplingTarget);

  fgYieldCouplingTAX = new TGraph();
  fgYieldCouplingTAX->SetNameTitle("CouplingScan/YieldCouplingTAX", "Yield per POT vs coupling, for TAX-produced events");
  BookHisto(fgYieldCouplingTAX);

  // Mass scan

  BookHisto("MassScan/hReachMass",  new TH2D("ReachMass", "Probability of N reaching the FV vs N mass",    10, 0.3, 1.2, 1000, -0.1, 1.1));
  BookHisto("MassScan/hDecayMass",  new TH2D("DecayMass", "Probability of N decaying in the FV vs N mass", 10, 0.3, 1.2, 1000, -0.1, 1.1));
  BookHisto("MassScan/hWeightMass", new TH2D("WeightMass", "N weight vs N mass",                           10, 0.3, 1.2, 1000, 1.E-25, 1.E-5));

  fgGammaTotMass = new TGraph();
  fgGammaTotMass->SetNameTitle("MassScan/GammaTotMass", "N total decay width vs N mass");
  BookHisto(fgGammaTotMass);

  fgTauMass = new TGraph();
  fgTauMass->SetNameTitle("MassScan/TauMass", "N lifetime vs N mass");
  BookHisto(fgTauMass);

  fgAccMass = new TGraph();
  fgAccMass->SetNameTitle("MassScan/AccMass", "Acceptance vs N mass");
  BookHisto(fgAccMass);

  fgAccMassTarget = new TGraph();
  fgAccMassTarget->SetNameTitle("MassScan/AccMassTarget", "Acceptance per POT vs N mass, for target-produced events");
  BookHisto(fgAccMassTarget);

  fgAccMassTAX = new TGraph();
  fgAccMassTAX->SetNameTitle("MassScan/AccMassTAX", "Acceptance per POT vs N mass, for TAX-produced events");
  BookHisto(fgAccMassTAX);

  fgYieldMass = new TGraph();
  fgYieldMass->SetNameTitle("MassScan/YieldMass", "Yield per POT vs N mass");
  BookHisto(fgYieldMass);

  fgYieldMassTarget = new TGraph();
  fgYieldMassTarget->SetNameTitle("MassScan/YieldMassTarget", "Yield per POT vs N mass, for target-produced events");
  BookHisto(fgYieldMassTarget);

  fgYieldMassTAX = new TGraph();
  fgYieldMassTAX->SetNameTitle("MassScan/YieldMassTAX", "Yield per POT vs N mass, for TAX-produced events");
  BookHisto(fgYieldMassTAX);

  // Total scan

  fgExclusion = new TGraph();
  fgExclusion->SetNameTitle("TotalScan/Exclusion", "Sensitivity vs N mass and N mass");
  BookHisto(fgExclusion);
}

void HeavyNeutrinoScan::Process(Int_t) {

  Double_t fCoupling          = -999.;  
  Double_t MN                 = 0.;
  Double_t HNLTau             = 0.;
  Double_t gammaTot           = 0.;
  Double_t NDecayProb         = 0.;
  Double_t NReachProb         = 0.;
  Double_t Weight             = 0.;
  Bool_t isGood               = false;

  if (fEnableCouplingScan == false) {
    fCouplingStart = TMath::Log10(fUSquared);
    fCouplingStop = fCouplingStart;
  }

  // Scan on the coupling                                                                       
  
  for(Int_t couplingIndex  = fCouplingStart*10; couplingIndex <= fCouplingStop*10; couplingIndex += fCouplingStep*10) {
    fCoupling = couplingIndex/10.;

    fUSquared = TMath::Power(10, fCoupling);
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
	gammaTot = GammaTot(MN);
	NReachProb = Weights[i]["ReachProb"];
	NDecayProb = Weights[i]["DecayProb"];
	Weight = Weights[i]["Weight"];    
	fMasses[round(MN)] = round(MN);
	fGammaTot[round(MN)][fCoupling] = gammaTot;
	fTau[round(MN)][fCoupling] = HNLTau;
	fSumAll[round(MN)][fCoupling] += Weight;

	if (fNevents[round(MN)].count(fCoupling) == 0)
	  fNevents[round(MN)][fCoupling] = 0;
	fNevents[round(MN)][fCoupling]++;
	errorCounter++;

	if (Weights[i]["ProdProb"] == fDBeProdProb) {
	  if (fNeventsTarget[round(MN)].count(fCoupling) == 0)
	    fNeventsTarget[round(MN)][fCoupling] = 0;
	  fNeventsTarget[round(MN)][fCoupling]++;
	  fSumAllTarget[round(MN)][fCoupling] += Weight;
	  errorCounterTarget++;
	}
	else if (Weights[i]["ProdProb"] == fDCuProdProb) {
	  if (fNeventsTAX[round(MN)].count(fCoupling) == 0)
	    fNeventsTAX[round(MN)][fCoupling] = 0;
	  fNeventsTAX[round(MN)][fCoupling]++;
	  fSumAllTAX[round(MN)][fCoupling] += Weight;
	  errorCounterTAX++;
	}

	FillHisto("CouplingScan/hReachCoupling",  fCoupling, NReachProb);
	FillHisto("CouplingScan/hDecayCoupling",  fCoupling, NDecayProb);
	FillHisto("CouplingScan/hWeightCoupling", fCoupling, Weight);	
	FillHisto("MassScan/hReachMass",          MN/1000.,  NReachProb);
	FillHisto("MassScan/hDecayMass",          MN/1000.,  NDecayProb);
	FillHisto("MassScan/hWeightMass",         MN/1000.,  Weight);

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

  if (GetWithMC()) {
    Event *evt = GetMCEvent();
    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);
      if (p->GetParentID() == -1 && p->GetPDGcode() == 999) {
	FillHisto("SingleValue/hZDProd",    p->GetPosAtCheckPoint(0).z());
	FillHisto("SingleValue/hZDDecay",   p->GetProdPos().Z());
	FillHisto("SingleValue/hDTheta",    p->GetPosAtCheckPoint(0).x());
	FillHisto("SingleValue/hDLambda",   p->GetPosAtCheckPoint(0).y());
	FillHisto("SingleValue/hDPath",     p->GetMomAtCheckPoint(0).X());
	FillHisto("SingleValue/hDMom",      p->GetMomAtCheckPoint(0).Y()/1000.);
	FillHisto("SingleValue/hZHNLDecay", p->GetEndPos().Z()/1000.);
	FillHisto("SingleValue/hHNLGamma",  p->GetInitial4Momentum().Gamma());
	FillHisto("SingleValue/hHNLTheta",  p->GetMomAtCheckPoint(0).Z());
	FillHisto("SingleValue/hHNLMom",    p->GetMomAtCheckPoint(0).T()/1000.);
      }
    }
  }

  //Bool_t IsHNLGood = *(Bool_t*)GetOutput("HeavyNeutrino.Output");
  //REMOVE
  Bool_t IsHNLGood = true;
  
  if (IsHNLGood == true) {
    for(Int_t couplingIndex  = fCouplingStart*10; couplingIndex <= fCouplingStop*10; couplingIndex += fCouplingStep*10) {
      fCoupling = couplingIndex/10.;
      fUSquared = TMath::Power(10, fCoupling);
      fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
      fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
      fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;

      if (GetWithMC()) {
	Event *evt = GetMCEvent();
	std::vector<std::map<std::string, Double_t>> Weights = ComputeWeight(evt, fUSquared, fUeSquaredRatio, fUmuSquaredRatio, fUtauSquaredRatio, fLInitialFV, fLFV);

	for (UInt_t i = 0; i < Weights.size(); i++) {
	  isGood = Weights[i]["IsGood"];
	  if (isGood == true) {
	    Weight = Weights[i]["Weight"];
	    fSumGood[round(MN)][fCoupling] += Weight;
	    if (Weights[i]["ProdProb"] == fDBeProdProb) 
	      fSumGoodTarget[round(MN)][fCoupling] += Weight;
	    else if (Weights[i]["ProdProb"] == fDCuProdProb) 
	      fSumGoodTAX[round(MN)][fCoupling] += Weight;
	  }
	}
      }
      
      if (errorCounter%errorStep == 0)
	errorFile << round(MN) << "\t" << fCoupling << "\t" << fSumAll[round(MN)][fCoupling] << "\t" << fSumGood[round(MN)][fCoupling] << "\t" << fNevents[round(MN)][fCoupling] << endl;
      if (errorCounterTarget%errorStep == 0)
	errorFileTarget << round(MN) << "\t" << fCoupling << "\t" << fSumAllTarget[round(MN)][fCoupling] << "\t" << fSumGoodTarget[round(MN)][fCoupling] << "\t" << fNeventsTarget[round(MN)][fCoupling] << endl;
      if (errorCounterTAX%errorStep == 0)
	errorFileTAX << round(MN) << "\t" << fCoupling << "\t" << fSumAllTAX[round(MN)][fCoupling] << "\t" << fSumGoodTAX[round(MN)][fCoupling] << "\t" << fNeventsTAX[round(MN)][fCoupling] << endl;
    }
  }
}

void HeavyNeutrinoScan::EndOfJobUser() {
  
  // Retrieve histos

  fhZDProd           = (TH1D*)fHisto.GetTH1("SingleValue/hZDProd");
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
  fhWeightCoupling   = (TH2D*)fHisto.GetTH2("CouplingScan/hWeightCoupling");
  fhReachMass        = (TH2D*)fHisto.GetTH2("MassScan/hReachMass");
  fhDecayMass        = (TH2D*)fHisto.GetTH2("MassScan/hDecayMass");
  fhWeightMass       = (TH2D*)fHisto.GetTH2("MassScan/hWeightMass");

  // X axis title

  fhZDProd          ->GetXaxis()->SetTitle("Z [mm]");
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
  fhWeightCoupling  ->GetXaxis()->SetTitle("Log of coupling");
  fhReachMass       ->GetXaxis()->SetTitle("N mass [GeV]");
  fhDecayMass       ->GetXaxis()->SetTitle("N mass [GeV]");
  fhWeightMass      ->GetXaxis()->SetTitle("N mass [GeV]");

  // Y axis title

  fhReachCoupling   ->GetYaxis()->SetTitle("Reach probability");
  fhDecayCoupling   ->GetYaxis()->SetTitle("Decay probability");
  fhWeightCoupling  ->GetYaxis()->SetTitle("Weight");
  fhReachMass       ->GetYaxis()->SetTitle("Reach probability");
  fhDecayMass       ->GetYaxis()->SetTitle("Decay probability");
  fhWeightMass      ->GetYaxis()->SetTitle("Weight");

  // Acceptance computation                                                                       

  Double_t Coupling          = 0.;
  Double_t MN                = 0.;
  Int_t counter              = 0;
  Int_t couplingCounter      = 0;
  Int_t massCounter          = 0;

  for (auto it = fMasses.begin(); it != fMasses.end(); it++) {
    MN = it->first;
    for (auto it1 = fCouplings.begin(); it1 != fCouplings.end(); it1++) {
      Coupling = it1->first;
      fAcc[MN][Coupling]         =       fSumGood[MN][Coupling]/       fSumAll[MN][Coupling];
      fAccTarget[MN][Coupling]   = fSumGoodTarget[MN][Coupling]/ fSumAllTarget[MN][Coupling];
      fAccTAX[MN][Coupling]      =    fSumGoodTAX[MN][Coupling]/    fSumAllTAX[MN][Coupling];
      fYield[MN][Coupling]       =       fSumGood[MN][Coupling]/      fNevents[MN][Coupling];
      fYieldTarget[MN][Coupling] = fSumGoodTarget[MN][Coupling]/fNeventsTarget[MN][Coupling];
      fYieldTAX[MN][Coupling]    =    fSumGoodTAX[MN][Coupling]/   fNeventsTAX[MN][Coupling];
      fProb[MN][Coupling]        =        fSumAll[MN][Coupling]/      fNevents[MN][Coupling];
      
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
      
      if (fYield[MN][Coupling]*1.E18 > 2.3)
	fgExclusion->SetPoint(counter, MN/1000., Coupling);
      else
	fgExclusion->SetPoint(counter, counter, counter);

      counter++;
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

  // Titles

  fgGammaTotCoupling   ->GetXaxis()->SetTitle("Log of coupling");
  fgTauCoupling        ->GetXaxis()->SetTitle("Log of coupling");
  fgAccCoupling        ->GetXaxis()->SetTitle("Log of coupling");
  fgAccCouplingTarget  ->GetXaxis()->SetTitle("Log of coupling");
  fgAccCouplingTAX     ->GetXaxis()->SetTitle("Log of coupling");
  fgYieldCoupling      ->GetXaxis()->SetTitle("Log of coupling");
  fgYieldCouplingTarget->GetXaxis()->SetTitle("Log of coupling");
  fgYieldCouplingTAX   ->GetXaxis()->SetTitle("Log of coupling");
  fgGammaTotCoupling   ->GetYaxis()->SetTitle("Total decay width [MeV]");
  fgTauCoupling        ->GetYaxis()->SetTitle("Lifetime [ns]");
  fgAccCoupling        ->GetYaxis()->SetTitle("Acceptance");
  fgAccCouplingTarget  ->GetYaxis()->SetTitle("Acceptance");
  fgAccCouplingTAX     ->GetYaxis()->SetTitle("Acceptance");
  fgYieldCoupling      ->GetYaxis()->SetTitle("Yield");
  fgYieldCouplingTarget->GetYaxis()->SetTitle("Yield");
  fgYieldCouplingTAX   ->GetYaxis()->SetTitle("Yield");
  fgGammaTotCoupling   ->SetLineColor(2);
  fgTauCoupling        ->SetLineColor(2);
  fgAccCoupling        ->SetLineColor(2);
  fgAccCouplingTarget  ->SetLineColor(8);
  fgAccCouplingTAX     ->SetLineColor(9);
  fgYieldCoupling      ->SetLineColor(2);
  fgYieldCouplingTarget->SetLineColor(8);
  fgYieldCouplingTAX   ->SetLineColor(9);
  fgGammaTotCoupling   ->SetLineWidth(3);
  fgTauCoupling        ->SetLineWidth(3);
  fgAccCoupling        ->SetLineWidth(3);
  fgAccCouplingTarget  ->SetLineWidth(3);
  fgAccCouplingTAX     ->SetLineWidth(3);
  fgYieldCoupling      ->SetLineWidth(3);
  fgYieldCouplingTarget->SetLineWidth(3);
  fgYieldCouplingTAX   ->SetLineWidth(3);
  fgGammaTotCoupling   ->Draw("AC");
  fgTauCoupling        ->Draw("AC");
  fgAccCoupling        ->Draw("AC");
  fgAccCouplingTarget  ->Draw("AC");
  fgAccCouplingTAX     ->Draw("AC");
  fgYieldCoupling      ->Draw("AC");
  fgYieldCouplingTarget->Draw("AC");
  fgYieldCouplingTAX   ->Draw("AC");
  fgGammaTotCoupling   ->Write();
  fgTauCoupling        ->Write();
  fgAccCoupling        ->Write();
  fgAccCouplingTarget  ->Write();
  fgAccCouplingTAX     ->Write();
  fgYieldCoupling      ->Write();
  fgYieldCouplingTarget->Write();
  fgYieldCouplingTAX   ->Write();

  fgGammaTotMass   ->GetXaxis()->SetTitle("N mass [GeV]");
  fgTauMass        ->GetXaxis()->SetTitle("N mass [GeV]");
  fgAccMass        ->GetXaxis()->SetTitle("N mass [GeV]");
  fgAccMassTarget  ->GetXaxis()->SetTitle("N mass [GeV]");
  fgAccMassTAX     ->GetXaxis()->SetTitle("N mass [GeV]");
  fgYieldMass      ->GetXaxis()->SetTitle("N mass [GeV]");
  fgYieldMassTarget->GetXaxis()->SetTitle("N mass [GeV]");
  fgYieldMassTAX   ->GetXaxis()->SetTitle("N mass [GeV]");
  fgGammaTotMass   ->GetYaxis()->SetTitle("Total decay width [MeV]");
  fgTauMass        ->GetYaxis()->SetTitle("Lifetime [ns]");
  fgAccMass        ->GetYaxis()->SetTitle("Acceptance");
  fgAccMassTarget  ->GetYaxis()->SetTitle("Acceptance");
  fgAccMassTAX     ->GetYaxis()->SetTitle("Acceptance");
  fgYieldMass      ->GetYaxis()->SetTitle("Yield");
  fgYieldMassTarget->GetYaxis()->SetTitle("Yield");
  fgYieldMassTAX   ->GetYaxis()->SetTitle("Yield");
  fgGammaTotMass   ->SetLineColor(2);
  fgTauMass        ->SetLineColor(2);
  fgAccMass        ->SetLineColor(2);
  fgAccMassTarget  ->SetLineColor(8);
  fgAccMassTAX     ->SetLineColor(9);
  fgYieldMass      ->SetLineColor(2);
  fgYieldMassTarget->SetLineColor(8);
  fgYieldMassTAX   ->SetLineColor(9);
  fgGammaTotMass   ->SetLineWidth(3);
  fgTauMass        ->SetLineWidth(3);
  fgAccMass        ->SetLineWidth(3);
  fgAccMassTarget  ->SetLineWidth(3);
  fgAccMassTAX     ->SetLineWidth(3);
  fgYieldMass      ->SetLineWidth(3);
  fgYieldMassTarget->SetLineWidth(3);
  fgYieldMassTAX   ->SetLineWidth(3);
  fgGammaTotMass   ->Draw("AC");
  fgTauMass        ->Draw("AC");
  fgAccMass        ->Draw("AC");
  fgAccMassTarget  ->Draw("AC");
  fgAccMassTAX     ->Draw("AC");
  fgYieldMass      ->Draw("AC");
  fgYieldMassTarget->Draw("AC");
  fgYieldMassTAX   ->Draw("AC");
  fgGammaTotMass   ->Write();
  fgTauMass        ->Write();
  fgAccMass        ->Write();
  fgAccMassTarget  ->Write();
  fgAccMassTAX     ->Write();
  fgYieldMass      ->Write();
  fgYieldMassTarget->Write();
  fgYieldMassTAX   ->Write();

  fgExclusion->GetXaxis()->SetTitle("N mass [GeV]");
  fgExclusion->GetYaxis()->SetTitle("Log of coupling");
  fgExclusion->SetLineColor(2);
  fgExclusion->SetLineWidth(3);
  fgExclusion->Draw("AC");
  fgExclusion->Write();

  SaveAllPlots();

  return;
}

HeavyNeutrinoScan::~HeavyNeutrinoScan() {

  errorFile.close();
  errorFileTarget.close();
  errorFileTAX.close();

  fhZDProd              = nullptr;
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
  fhWeightMass          = nullptr;
  fgAccMass             = nullptr;
  fgAccMassTarget       = nullptr;
  fgAccMassTAX          = nullptr;
  fgYieldMass           = nullptr;
  fgYieldMassTarget     = nullptr;
  fgYieldMassTAX        = nullptr;
  fgGammaTotMass        = nullptr;
  fgTauMass             = nullptr;
  fgExclusion           = nullptr;
}
