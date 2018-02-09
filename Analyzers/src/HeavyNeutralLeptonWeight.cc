// ---------------------------------------------------------------                                
//                                                                                                      
// History:                                                                                             
//                                                                                                      
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018                               
//                                                                                                      
// ---------------------------------------------------------------                                      

/// \class HeavyNeutralLeptonWeight
/// \Brief                                                  
/// Compute the weight associated to each heavy neutral lepton in a MC sample, as an output for the user
/// \EndBrief                                                                                           
/// \Detailed
/// For each HNL in the MC sample, a weight is computed for the user.
/// This analyzer makes use of an external library implemented in PhysicsObjects, called HNLFunctions.
/// The value of the squared HNL coupling and the values of the ratios between specific-flavour        
/// couplings can be either set as external parameters from command line or taken as default values.   
/// For example, if the user sets U2 = 1.E-10, and UeRatio = 5., UmuRatio = 1., UtauRatio = 3.5,       
/// the specific-flavour coupling values will be: Ue2 = 5.25E-11, Umu2 = 1.05E-11, Utau2 = 3.68E-11.
/// The values of the beginning of the fiducial volume and its length (identical to the ones set 
/// in the MC macro for production) can be either set as external parameters from command line 
/// or taken as default values.
/// For example, if the user assigns 100000. to the beginning of the FV and 80000. to its length, the FV
/// will begin at 100 m from the target centre and end at 180 m.
/// A vector containing all the weight associated to the HNLs in a certain event is produced
/// as output of the analyzer. Each vector element is a triplet containing HNL mass, 
/// coupling with the lepton the HNL is produced with, and weight. 
/// The output can be accessed from another analyzer in the following way:
/// \code
/// std::vector<std::map<std::string, Double_t>> Weights = 
///   *(std::vector<std::map<std::string, Double_t>>*)
///   GetOutput("HeavyNeutralLeptonWeight.Output");
/// \endcode
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
#include "HeavyNeutralLeptonWeight.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "BeamParameters.hh"
#include "HNLFunctions.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

/// \class HeavyNeutralLeptonWeight

HeavyNeutralLeptonWeight::HeavyNeutralLeptonWeight(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutralLeptonWeight") {

  if (!GetIsTree()) return;

  RequestAllMCTrees();
  RequestAllRecoTrees();

  AddParam("USquared", &fUSquared, 1.E-6);
  AddParam("UeSquaredRatio", &fUeSquaredRatio, 1.);
  AddParam("UmuSquaredRatio", &fUmuSquaredRatio, 16.);
  AddParam("UtauSquaredRatio", &fUtauSquaredRatio, 3.8);
  AddParam("LInitialFV", &fLInitialFV, 102500.);
  AddParam("LFV", &fLFV, 77500.);

  fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
  fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
  fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;
}

void HeavyNeutralLeptonWeight::InitOutput() {

  RegisterOutput("Output", &fWeightContainer);
}

void HeavyNeutralLeptonWeight::Process(Int_t) {

  Double_t MN             = 0.;
  Double_t HNLTau         = 0.;
  Double_t NDecayProb     = 0.;
  Double_t NReachProb     = 0.;
  Double_t LReach         = 0.;
  Double_t LeptonUSquared = 0.;
  Double_t DecayFactor    = 0.;
  Double_t ProdFactor     = 0.;
  Double_t Weight         = 0.;
  Double_t DProdProb      = 0.;
  Bool_t isGood           = false;
  TVector3 point1;
  TVector3 point2;
  TVector3 momentum1;

  fWeightContainer.clear();

  if (GetWithMC()) {
    Event *evt = GetMCEvent();
    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);      
      if (p->GetParentID() == -1 && p->GetPDGcode() == 999) {
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
	else if (p->GetParticleName() == "DS->Ntau" || p->GetParticleName().Contains("nu_tau") || (p->GetParticleName().Contains("tau") && (p->GetParticleName().Contains("rho") || p->GetParticleName().Contains("pi")|| p->GetParticleName().Contains("K"))))
	  LeptonUSquared = fUtauSquared;

	if (p->GetParticleName().Contains("tau->")) {
	  if (p->GetParticleName().Contains("DS"))
	    ProdFactor = fDStoTauBR;
	  else
	    ProdFactor = fDtoTauBR;
	}
	else
	  ProdFactor = 1.;

	if (p->GetProdPos().Z() < fTAXDistance)
	  DProdProb = fDBeProdProb;
	else if (p->GetProdPos().Z() >= fTAXDistance)
	  DProdProb = fDCuProdProb;

	if (p->GetEndProcessName() == "good")
	  isGood = true;

	Weight = DProdProb*fDDecayProb*NReachProb*NDecayProb*DecayFactor*ProdFactor*LeptonUSquared;

	std::map<std::string, Double_t> SingleHNL;

	SingleHNL["Mass"] = MN;
	SingleHNL["Lifetime"] = HNLTau;
	SingleHNL["DecayFactor"] = DecayFactor;
	SingleHNL["ReachProb"] = NReachProb;
	SingleHNL["DecayProb"] = NDecayProb;
	SingleHNL["LeptonUSquared"] = LeptonUSquared;
	SingleHNL["ProdFactor"] = ProdFactor;
	SingleHNL["ProdProb"] = DProdProb;
	SingleHNL["IsGood"] = isGood;
	SingleHNL["Weight"] = Weight;

	fWeightContainer.push_back(SingleHNL);
      }
    }
  }
}

