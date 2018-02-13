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
  AddParam("CouplingStep", &fCouplingStep, 0.1);
  AddParam("EnableCouplingScan", &fEnableCouplingScan, true);

  fN = round((std::abs(fCouplingStop-fCouplingStart))/fCouplingStep);
  fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
  fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
  fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;

  // Scan histos

  fhAcc                 = nullptr;
  fhYield               = nullptr;
  fhYieldTarget         = nullptr;
  fhYieldTAX            = nullptr;
  fhReachCoupling       = nullptr;
  fhDecayCoupling       = nullptr;
  fhWeightCoupling      = nullptr;
  fgAccCoupling         = nullptr;
  fgYieldCoupling       = nullptr;
  fgYieldCouplingTarget = nullptr;
  fgYieldCouplingTAX    = nullptr;
  fgGammaTotCoupling    = nullptr;
  fgTauCoupling         = nullptr;
  fhReachMass           = nullptr;
  fhDecayMass           = nullptr;
  fhWeightMass          = nullptr;
  fgAccMass             = nullptr;
  fgYieldMass           = nullptr;
  fgYieldMassTarget     = nullptr;
  fgYieldMassTAX        = nullptr;
  fgGammaTotMass        = nullptr;
  fgTauMass             = nullptr;
  fgExclusion           = nullptr;
}

void HeavyNeutrinoScan::InitHist() {

  // One value of mass and coupling

  BookHisto("SingleValue/hAcc",         new TH1D("Acc", "Acceptance for one value of N mass and coupling",                                        1000, 1.E-30, 1.));
  BookHisto("SingleValue/hYield",       new TH1D("Yield", "Yield per POT for one value of N mass and coupling",                                   1000, 1.E-50, 1.E-10));

  BookHisto("SingleValue/hYieldTarget", new TH1D("YieldTarget", "Yield per POT for one value of N mass and coupling, for target-produced events", 1000, 1.E-50, 1.E-10));
  BookHisto("SingleValue/hYieldTAX",    new TH1D("YieldTAX", "Yield per POT for one value of N mass and coupling, for TAX-produced events",       1000, 1.E-50, 1.E-10));

  // Coupling scan 

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

  fgYieldCouplingTarget = new TGraph();
  fgYieldCouplingTarget->SetNameTitle("CouplingScan/YieldCouplingTarget", "Yield per POT vs coupling, for target-produced events");
  BookHisto(fgYieldCouplingTarget);

  fgYieldCouplingTAX = new TGraph();
  fgYieldCouplingTAX->SetNameTitle("CouplingScan/YieldCouplingTAX", "Yield per POT vs coupling, for TAX-produced events");
  BookHisto(fgYieldCouplingTAX);

  // Mass scan

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

  fgYieldMassTarget = new TGraph();
  fgYieldMassTarget->SetNameTitle("MassScan/YieldMassTarget", "Yield per POT vs coupling, for target-produced events");
  BookHisto(fgYieldMassTarget);

  fgYieldMassTAX = new TGraph();
  fgYieldMassTAX->SetNameTitle("MassScan/YieldMassTAX", "Yield per POT vs coupling, for TAX-produced events");
  BookHisto(fgYieldMassTAX);

  fgExclusion = new TGraph();
  fgExclusion->SetNameTitle("TotalScan/Exclusion", "Sensitivity vs N mass and coupling");
  BookHisto(fgExclusion);
}

void HeavyNeutrinoScan::Process(Int_t) {
  
  Double_t MN             = 0.;
  Double_t HNLTau         = 0.;
  Double_t gammaTot       = 0.;
  Double_t NDecayProb     = 0.;
  Double_t NReachProb     = 0.;
  Double_t Weight         = 0.;
  Bool_t isGood           = false;

  if (fEnableCouplingScan == false) {

    fCouplingStart = fUSquared;
    fCouplingStep = 0.1;
    fCouplingStop = fUSquared + 2*fCouplingStep;
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

	if (Weights[i]["DProdProb"] == fDBeProdProb) {
	  if (fNeventsTarget[round(MN)].count(fCoupling) == 0)
	    fNeventsTarget[round(MN)][fCoupling] = 0;
	  fNeventsTarget[round(MN)][fCoupling]++;
	  fSumAllTarget[round(MN)][fCoupling] += Weight;
	}
	else if (Weights[i]["DProdProb"] == fDCuProdProb) {
	  if (fNeventsTAX[round(MN)].count(fCoupling) == 0)
	    fNeventsTAX[round(MN)][fCoupling] = 0;
	  fNeventsTAX[round(MN)][fCoupling]++;
	  fSumAllTAX[round(MN)][fCoupling] += Weight;
	}
	
	FillHisto("CouplingScan/hReachCoupling",  fCoupling, NReachProb);
	FillHisto("CouplingScan/hDecayCoupling",  fCoupling, NDecayProb);
	FillHisto("CouplingScan/hWeightCoupling", fCoupling, Weight);	
	FillHisto("MassScan/hReachMass",          fCoupling, NReachProb);
	FillHisto("MassScan/hDecayMass",          fCoupling, NDecayProb);
	FillHisto("MassScan/hWeightMass",         fCoupling, Weight);
      }
    }
  }
  
  Bool_t IsHNLGood = *(Bool_t*)GetOutput("HeavyNeutrino.Output");
  
  if (IsHNLGood == true) {
    for(Double_t fCoupling = fCouplingStart; fCoupling < fCouplingStop-fCouplingStep; fCoupling += fCouplingStep) {
      if (GetWithMC()) {
	Event *evt = GetMCEvent();
	std::vector<std::map<std::string, Double_t>> Weights = ComputeWeight(evt, fUSquared, fUeSquaredRatio, fUmuSquaredRatio, fUtauSquaredRatio, fLInitialFV, fLFV);

	for (UInt_t i = 0; i < Weights.size(); i++) {
	  isGood = Weights[i]["IsHNLGood"];
	  if (isGood >= 1.) {
	    Weight = Weights[i]["Weight"];
	    fSumGood[round(MN)][fCoupling] += Weight;
	    if (Weights[i]["DProdProb"] == fDBeProdProb) 
	      fSumGoodTarget[round(MN)][fCoupling] += Weight;
	    else if (Weights[i]["DProdProb"] == fDCuProdProb) 
	      fSumGoodTAX[round(MN)][fCoupling] += Weight;
	  }
	}
      }
    }
  }
}

void HeavyNeutrinoScan::EndOfJobUser() {
  
  // Retrieve scan histos

  fhAcc              = (TH1D*)fHisto.GetTH1("SingleValue/hAcc");
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

  fhAcc             ->GetXaxis()->SetTitle("Acceptance");
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

  Double_t Coupling = 0.;
  Double_t MN       = 0.;
  Int_t counter     = 0;

  for (auto it = fMasses.begin(); it != fMasses.end(); it++) {
    MN = it->first;
    for (auto it1 = fCouplings.begin(); it1 != fCouplings.end(); it1++) {
      Coupling = it1->first;
      fAcc[MN][Coupling]         =       fSumGood[MN][Coupling]/       fSumAll[MN][Coupling];
      fProb[MN][Coupling]        =        fSumAll[MN][Coupling]/      fNevents[MN][Coupling];
      fYield[MN][Coupling]       =       fSumGood[MN][Coupling]/      fNevents[MN][Coupling];
      fYieldTarget[MN][Coupling] = fSumGoodTarget[MN][Coupling]/fNeventsTarget[MN][Coupling];
      fYieldTAX[MN][Coupling]    =    fSumGoodTAX[MN][Coupling]/   fNeventsTAX[MN][Coupling];

      FillHisto("SingleValue/hAcc",                 fAcc[MN][Coupling]);
      FillHisto("SingleValue/hYield",             fYield[MN][Coupling]);
      FillHisto("SingleValue/hYieldTarget", fYieldTarget[MN][Coupling]);
      FillHisto("SingleValue/hYieldTAX",       fYieldTAX[MN][Coupling]);

      fgGammaTotCoupling   ->SetPoint(counter, Coupling, fGammaTot   [MN][Coupling]);
      fgTauCoupling        ->SetPoint(counter, Coupling, fTau        [MN][Coupling]);
      fgAccCoupling        ->SetPoint(counter, Coupling, fAcc        [MN][Coupling]);
      fgYieldCoupling      ->SetPoint(counter, Coupling, fYield      [MN][Coupling]);
      fgYieldCouplingTarget->SetPoint(counter, Coupling, fYieldTarget[MN][Coupling]);
      fgYieldCouplingTAX   ->SetPoint(counter, Coupling, fYieldTAX   [MN][Coupling]);
      fgGammaTotMass       ->SetPoint(counter, MN/1000., fGammaTot   [MN][Coupling]);
      fgTauMass            ->SetPoint(counter, MN/1000., fTau        [MN][Coupling]);
      fgAccMass            ->SetPoint(counter, MN/1000., fAcc        [MN][Coupling]);
      fgYieldMass          ->SetPoint(counter, MN/1000., fYield      [MN][Coupling]);
      fgYieldMassTarget    ->SetPoint(counter, MN/1000., fYieldTarget[MN][Coupling]);
      fgYieldMassTAX       ->SetPoint(counter, MN/1000., fYieldTAX   [MN][Coupling]);

      if (fYield[MN][Coupling]*1.E18 > 2.3)
	fgExclusion->SetPoint(counter, MN/1000., Coupling);
      counter++;
    }
  }

  // Titles

  fgGammaTotCoupling   ->GetXaxis()->SetTitle("Log of coupling");
  fgTauCoupling        ->GetXaxis()->SetTitle("Log of coupling");
  fgAccCoupling        ->GetXaxis()->SetTitle("Log of coupling");
  fgYieldCoupling      ->GetXaxis()->SetTitle("Log of coupling");
  fgYieldCouplingTarget->GetXaxis()->SetTitle("Log of coupling");
  fgYieldCouplingTAX   ->GetXaxis()->SetTitle("Log of coupling");
  fgGammaTotCoupling   ->GetYaxis()->SetTitle("Total decay width [MeV]");
  fgTauCoupling        ->GetYaxis()->SetTitle("Lifetime [ns]");
  fgAccCoupling        ->GetYaxis()->SetTitle("Acceptance");
  fgYieldCoupling      ->GetYaxis()->SetTitle("Yield");
  fgYieldCouplingTarget->GetYaxis()->SetTitle("Yield");
  fgYieldCouplingTAX   ->GetYaxis()->SetTitle("Yield");
  fgGammaTotCoupling   ->SetLineColor(2);
  fgTauCoupling        ->SetLineColor(2);
  fgAccCoupling        ->SetLineColor(2);
  fgYieldCoupling      ->SetLineColor(2);
  fgYieldCouplingTarget->SetLineColor(8);
  fgYieldCouplingTAX   ->SetLineColor(9);
  fgGammaTotCoupling   ->SetLineWidth(3);
  fgTauCoupling        ->SetLineWidth(3);
  fgAccCoupling        ->SetLineWidth(3);
  fgYieldCoupling      ->SetLineWidth(3);
  fgYieldCouplingTarget->SetLineWidth(3);
  fgYieldCouplingTAX   ->SetLineWidth(3);
  fgGammaTotCoupling   ->Draw("AC");
  fgTauCoupling        ->Draw("AC");
  fgAccCoupling        ->Draw("AC");
  fgYieldCoupling      ->Draw("AC");
  fgYieldCouplingTarget->Draw("AC");
  fgYieldCouplingTAX   ->Draw("AC");
  fgGammaTotCoupling   ->Write();
  fgTauCoupling        ->Write();
  fgAccCoupling        ->Write();
  fgYieldCoupling      ->Write();
  fgYieldCouplingTarget->Write();
  fgYieldCouplingTAX   ->Write();

  fgGammaTotMass   ->GetXaxis()->SetTitle("N mass [GeV]");
  fgTauMass        ->GetXaxis()->SetTitle("N mass [GeV]");
  fgAccMass        ->GetXaxis()->SetTitle("N mass [GeV]");
  fgYieldMass      ->GetXaxis()->SetTitle("N mass [GeV]");
  fgYieldMassTarget->GetXaxis()->SetTitle("N mass [GeV]");
  fgYieldMassTAX   ->GetXaxis()->SetTitle("N mass [GeV]");
  fgGammaTotMass   ->GetYaxis()->SetTitle("Total decay width [MeV]");
  fgTauMass        ->GetYaxis()->SetTitle("Lifetime [ns]");
  fgAccMass        ->GetYaxis()->SetTitle("Acceptance");
  fgYieldMass      ->GetYaxis()->SetTitle("Yield");
  fgYieldMassTarget->GetYaxis()->SetTitle("Yield");
  fgYieldMassTAX   ->GetYaxis()->SetTitle("Yield");
  fgGammaTotMass   ->SetLineColor(2);
  fgTauMass        ->SetLineColor(2);
  fgAccMass        ->SetLineColor(2);
  fgYieldMass      ->SetLineColor(2);
  fgYieldMassTarget->SetLineColor(8);
  fgYieldMassTAX   ->SetLineColor(9);
  fgGammaTotMass   ->SetLineWidth(3);
  fgTauMass        ->SetLineWidth(3);
  fgAccMass        ->SetLineWidth(3);
  fgYieldMass      ->SetLineWidth(3);
  fgYieldMassTarget->SetLineWidth(3);
  fgYieldMassTAX   ->SetLineWidth(3);
  fgGammaTotMass   ->Draw("AC");
  fgTauMass        ->Draw("AC");
  fgAccMass        ->Draw("AC");
  fgYieldMass      ->Draw("AC");
  fgYieldMassTarget->Draw("AC");
  fgYieldMassTAX   ->Draw("AC");
  fgGammaTotMass   ->Write();
  fgTauMass        ->Write();
  fgAccMass        ->Write();
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

  fhAcc                 = nullptr;
  fhYield               = nullptr;
  fhYieldTarget         = nullptr;
  fhYieldTAX            = nullptr;
  fhReachCoupling       = nullptr;
  fhDecayCoupling       = nullptr;
  fhWeightCoupling      = nullptr;
  fgAccCoupling         = nullptr;
  fgYieldCoupling       = nullptr;
  fgYieldCouplingTarget = nullptr;
  fgYieldCouplingTAX    = nullptr;
  fgGammaTotCoupling    = nullptr;
  fgTauCoupling         = nullptr;
  fhReachMass           = nullptr;
  fhDecayMass           = nullptr;
  fhWeightMass          = nullptr;
  fgAccMass             = nullptr;
  fgYieldMass           = nullptr;
  fgYieldMassTarget     = nullptr;
  fgYieldMassTAX        = nullptr;
  fgGammaTotMass        = nullptr;
  fgTauMass             = nullptr;
  fgExclusion           = nullptr;
}
