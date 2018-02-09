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
#include "HeavyNeutrinoAllIncluded.hh"
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

/// \class HeavyNeutrinoAllIncluded

HeavyNeutrinoAllIncluded::HeavyNeutrinoAllIncluded(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrinoAllIncluded") {

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

void HeavyNeutrinoAllIncluded::InitHist() {

  BookHisto("CouplingScan/hReachCoupling",    new TH2D("ReachCoupling", "Probability of N reaching the FV vs coupling",    fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("CouplingScan/hDecayCoupling",    new TH2D("DecayCoupling", "Probability of N decaying in the FV vs coupling", fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("CouplingScan/hWeightCoupling",   new TH2D("WeightCoupling", "N weight vs coupling",                           fN, fCouplingStart, fCouplingStop, 1000, 1.E-12, 1.E-9));

  fgGammaTotCoupling = new TGraph();
  fgGammaTotCoupling->SetNameTitle("CouplingScan/GammaTotCoupling", "N total decay width vs coupling");
  BookHisto(fgGammaTotCoupling);

  fgTauCoupling = new TGraph();
  fgTauCoupling->SetNameTitle("CouplingScan/TauCoupling", "N lifetime vs coupling");
  BookHisto(fgTauCoupling);

  fgAccCoupling = new TGraph();
  fgAccCoupling->SetNameTitle("CouplingScan/AccCoupling", "Acceptance vs coupling");
  BookHisto(fgAccCoupling);

  fgYieldCoupling = new TGraph();
  fgYieldCoupling->SetNameTitle("CouplingScan/YieldCoupling", "Yield per POT vs coupling");
  BookHisto(fgYieldCoupling);

  BookHisto("MassScan/hReachMass",    new TH2D("ReachMass", "Probability of N reaching the FV vs coupling",    fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("MassScan/hDecayMass",    new TH2D("DecayMass", "Probability of N decaying in the FV vs coupling", fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("MassScan/hWeightMass",   new TH2D("WeightMass", "N weight vs coupling",                           fN, fCouplingStart, fCouplingStop, 1000, 1.E-12, 1.E-9));

  fgGammaTotMass = new TGraph();
  fgGammaTotMass->SetNameTitle("MassScan/GammaTotMass", "N total decay width vs coupling");
  BookHisto(fgGammaTotMass);

  fgTauMass = new TGraph();
  fgTauMass->SetNameTitle("MassScan/TauMass", "N lifetime vs coupling");
  BookHisto(fgTauMass);

  fgAccMass = new TGraph();
  fgAccMass->SetNameTitle("MassScan/AccMass", "Acceptance vs coupling");
  BookHisto(fgAccMass);

  fgYieldMass = new TGraph();
  fgYieldMass->SetNameTitle("MassScan/YieldMass", "Yield per POT vs coupling");
  BookHisto(fgYieldMass);

  fgExclusion = new TGraph();
  fgExclusion->SetNameTitle("TotalScan/Exclusion", "Sensitivity vs N mass and coupling");
  BookHisto(fgExclusion);
}

void HeavyNeutrinoAllIncluded::Process(Int_t) {
  
  Double_t MN             = 0.;
  Double_t HNLTau         = 0.;
  Double_t gammaTot       = 0.;
  Double_t NDecayProb     = 0.;
  Double_t NReachProb     = 0.;
  Double_t Weight         = 0.;
  Bool_t isGood           = false;
  
  // Scan on the coupling                                                                       
  
  for(Double_t fCoupling = fCouplingStart; fCoupling < fCouplingStop-fCouplingStep; fCoupling += fCouplingStep) {
    
    fUSquared = TMath::Power(10, fCoupling);
    fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
    fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
    fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;
    fCouplings[fCoupling] = fCoupling;
    
    SetCouplingForWeight(fUSquared, fUeSquared, fUmuSquared, fUtauSquared);
    
    if (GetWithMC()) {
      std::vector<std::map<std::string, Double_t>> Weights = *(std::vector<std::map<std::string, Double_t>>*)GetOutput("HeavyNeutralLeptonWeight.Output");
      for (UInt_t i = 0; i < Weights.size(); i++) {
	MN =  Weights[i]["Mass"];
	HNLTau = Weights[i]["Lifetime"];
	gammaTot = GammaTot(MN);
	NReachProb = Weights[i]["ReachProb"];
	NDecayProb = Weights[i]["DecayProb"];
	Weight = Weights[i]["Weight"];    
	fMasses[round(MN)] = round(MN);
	if (fNevents[round(MN)].count(fCoupling) == 0)
	  fNevents[round(MN)][fCoupling] = 0;
	fNevents[round(MN)][fCoupling]++;
	fGammaTot[round(MN)][fCoupling] = gammaTot;
	fTau[round(MN)][fCoupling] = HNLTau;
	fSumAll[round(MN)][fCoupling] += Weight;
	
	FillHisto("CouplingScan/hReachCoupling",  fCoupling, NReachProb);
	FillHisto("CouplingScan/hDecayCoupling",  fCoupling, NDecayProb);
	FillHisto("CouplingScan/hWeightCoupling", fCoupling, Weight);
	
	FillHisto("MassScan/hReachMass",  fCoupling, NReachProb);
	FillHisto("MassScan/hDecayMass",  fCoupling, NDecayProb);
	FillHisto("MassScan/hWeightMass", fCoupling, Weight);
      }
    }
  }
  
  Bool_t IsHNLGood = *(Bool_t)GetOutput("HeavyNeutrino.Output");
  
  if (IsHNLGood == true) {
    for(Double_t fCoupling = fCouplingStart; fCoupling < fCouplingStop-fCouplingStep; fCoupling += fCouplingStep) {
      if (GetWithMC()) {
	std::vector<std::map<std::string, Double_t>> Weights = *(std::vector<std::map<std::string, Double_t>>*)GetOutput("HeavyNeutralLeptonWeight.Output");
	for (UInt_t i = 0; i < Weights.size(); i++) {
	  Weight = Weights[i]["Weight"];
	  fSumAll[round(MN)][fCoupling] += Weight;
	}
      }
    }
  }
}

void HeavyNeutrinoAllIncluded::EndOfJobUser() {
  
  // Retrieve scan histos

  fhReachCoupling    = (TH2D*)fHisto.GetTH2("hReachCoupling");
  fhDecayCoupling    = (TH2D*)fHisto.GetTH2("hDecayCoupling");
  fhWeightCoupling   = (TH2D*)fHisto.GetTH2("hWeightCoupling");

  fhReachMass    = (TH2D*)fHisto.GetTH2("hReachMass");
  fhDecayMass    = (TH2D*)fHisto.GetTH2("hDecayMass");
  fhWeightMass   = (TH2D*)fHisto.GetTH2("hWeightMass");

  // X axis title

  fhReachCoupling   ->GetXaxis()->SetTitle("Log of coupling");
  fhDecayCoupling   ->GetXaxis()->SetTitle("Log of coupling");
  fhWeightCoupling  ->GetXaxis()->SetTitle("Log of coupling");

  fhReachMass   ->GetXaxis()->SetTitle("N mass [GeV]");
  fhDecayMass   ->GetXaxis()->SetTitle("N mass [GeV]");
  fhWeightMass  ->GetXaxis()->SetTitle("N mass [GeV]");

  // Y axis title

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

  SaveAllPlots();

  return;
}

HeavyNeutrinoAllIncluded::~HeavyNeutrinoAllIncluded() {

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
