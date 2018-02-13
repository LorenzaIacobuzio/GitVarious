// ---------------------------------------------------------------                                
//                                                                                                      
// History:                                                                                             
//                                                                                                      
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018                               
//                                                                                                      
// ---------------------------------------------------------------                                      

/// \class HNLWeight
/// \Brief                                                  
/// Compute and return the weight associated to each heavy neutral lepton in a MC sample
/// \EndBrief                                                                                           
/// \Detailed
/// For each HNL in the MC sample, a weight is computed for the user.
/// This analyzer makes use of a ToolsLib, called HNLFunctions.
/// A pointer to the current processed event is needed as first argument of the function ComputeWeight.
/// The value of the squared HNL coupling and the values of the ratios between specific-flavour        
/// couplings are also passed as arguments of the tool.   
/// For example, if the user sets U2 = 1.E-10, and UeRatio = 5., UmuRatio = 1., UtauRatio = 3.5,       
/// the specific-flavour coupling values will be: Ue2 = 5.25E-11, Umu2 = 1.05E-11, Utau2 = 3.68E-11.
/// The values of the beginning of the fiducial volume and its length (must be identical to the ones set
/// in the MC macro for production) are also passed as arguments.
/// For example, if the user assigns 100000. to the beginning of the FV and 80000. to its length, the FV
/// will begin at 100 m from the target centre and end at 180 m.
/// A vector of maps is produced and returned. 
/// Each map contains quantities associated to each HNL (between brackets is the name by which the user
/// can retrieve them from the map): mass (Mass), lifetime (Lifetime), 
/// factor associated to its decay (DecayFactor), probability of reaching the FV (ReachProb) and
/// decaying into it (DecayProb), coupling with the lepton produced in pair with (LeptonUSquared), 
/// factor associated to its production (ProdFactor), probability of proton producing D meson in 
/// target/TAX (DProdProb), weight for each HNL (Weight), and a boolean to check if the HNL 
/// is associated to the two daughters of the event (IsGood).
/// The returned vector can be retrieved in the following way:
/// \code
/// std::vector<std::map<std::string, Double_t>> = 
/// ComputeWeight(evt, 1.E-10, 5., 1., 3.5, 100000., 80000.);
/// \endcode
///
/// \author Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch)                                               
/// \EndDetailed

#include "Event.hh"
#include "HNLFunctions.hh"
#include "HNLWeight.hh"

using namespace std;

/// \class HNLWeight

std::vector<std::map<std::string, Double_t>> fWeightContainer;

std::vector<std::map<std::string, Double_t>> ComputeWeight(Event* evt, Double_t USquared, Double_t UeSquaredRatio, Double_t UmuSquaredRatio, Double_t UtauSquaredRatio, Double_t LInitialFV, Double_t LFV) {

  fUSquared = USquared;
  fUeSquaredRatio = UeSquaredRatio;
  fUmuSquaredRatio = UmuSquaredRatio;
  fUtauSquaredRatio = UtauSquaredRatio;
  fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
  fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
  fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;
  fLInitialFV = LInitialFV;
  fLFV = LFV;

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
  Bool_t IsGood           = false;
  TVector3 point1;
  TVector3 point2;
  TVector3 momentum1;

  fWeightContainer.clear();
  
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
	IsGood = true;
      
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
      SingleHNL["IsGood"] = IsGood;
      SingleHNL["Weight"] = Weight;
      
      fWeightContainer.push_back(SingleHNL);
    }
  }

  return fWeightContainer;
}
