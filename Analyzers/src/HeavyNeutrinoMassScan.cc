#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include "HeavyNeutrinoMassScan.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "GeometricAcceptance.hh"
#include "TriggerConditions.hh"
#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/GaussLegendreIntegrator.h"
#include "TF1.h"
#include "TF2.h"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

#define TRIGGER_L0_PHYSICS_TYPE 1
#define MCTriggerMask 0xFF
#define LabelSize 0.05

/// \class HeavyNeutrinoMassScan

HeavyNeutrinoMassScan::HeavyNeutrinoMassScan(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrinoMassScan") {
  
  RequestAllMCTrees();
  RequestAllRecoTrees();
  RequestL0Data();
  RequestL1Data();

  AddParam("USquared", &fUSquared, 1.E-6);
  AddParam("UeSquaredRatio", &fUeSquaredRatio, 1.);
  AddParam("UmuSquaredRatio", &fUmuSquaredRatio, 16.);
  AddParam("UtauSquaredRatio", &fUtauSquaredRatio, 3.8);

  fCDAcomp     = new TwoLinesCDA();
  fDistcomp    = new PointLineDistance();
  fLAVMatching = new LAVMatching();
  fSAVMatching = new SAVMatching();

  // Scan variables

  fMN = 0.;

  for (Int_t i = 0; i < fN; i++) {
    fSumGood[i]   = 0.;
    fSumAll[i]    = 0.;
    fNevents[i]   = 0;
    fMasses[i]    = 0.;
    fAcc[i]       = 0.;
    fGammaTot[i]  = 0.;
    fTau[i]       = 0.;
    fProb[i]      = 0.;
    fYield[i]     = 0.;
  }

  // Masses                                                                                           

  fMe        = 0.511;
  fMmu       = 105.66;
  fMtau      = 1776.82;
  fMpi       = 139.57;
  fMpi0      = 134.98;
  fMrho      = 775.45;
  fMrho0     = 775.49;
  fMeta      = 547.86;
  fMetaprime = 957.78;
  fMD        = 1869.62;
  fMDS       = 1968.47;
  fMD0       = 1864.84;
  fMK        = 493.68;
  fMK0       = 497.65;
  fMp        = 938.27;
  fMKStar    = 891.76;
  fMK0Star   = 895.55;

  // Lifetimes                                                                                         

  fDlife   = 1.04E-3;
  fDSlife  = 5.E-4;
  fD0life  = 4.1E-4;
  ftaulife = 2.91E-4;

  // Constants                                                                                        

  fhc       = 197.327E-12; // MeV mm                                                    
  fcLight   = 299.792; //mm/ns                                                                         
  fGF       = 1.166E-11; // MeV^-2                                                           
  fPi       = 130.41;  // MeV                                                          
  fRho      = 1.04E5; // MeV^2                                                                      
  fD        = 222.6;
  fDS       = 280.1;
  fK        = 159.8;
  fEta      = 1.2*fPi;
  fEtaprime = -0.45*fPi;
  fsigmacc  = 2.3*75.; //mubarn at sqrt(s) = 82 GeV (400 GeV proton on Be(9) (mBe = 9*1 GeV), taken from Gaia's note                                                                                 

  // CKM                                                                                            
 
  fVcs = 0.9734;
  fVcd = 0.2252;
  fVud = 0.9743;
  fVus = 0.2253;

  // Form factors, pseudoscalar and vector mesons                                            

  fDK0   = 0.745; // f+                                                                      
  fDpi0  = 0.648;
  fD0K   = 0.736;
  fD0pi  = 0.637;
  fgDK0  = -0.495; // f-                                                                   
  fgDpi0 = -0.435;
  fgD0K  = fgDK0;
  fgD0pi = fgDpi0;

  fA0D  = 0.398;
  fA1D  = 0.47;
  fA2D  = -1.24;
  fVD   = 0.66;
  fA0D0 = 0.4;
  fA1D0 = 0.47;
  fA2D0 = -1.24;
  fVD0  = 0.66;

  // Fragmentation fractions                                                                  

  ffD  = 0.246;
  ffD0 = 0.565;
  ffDS = 0.08;

  // NA62 parameters                                                                    

  fpMom                   = 400000.; // MeV                                                
  fBeA                    = 4;
  fBeDensity              = 1.85; // g/cm3                                                       
  fpBeLambda              = 421.; // mm                                                            
  ftargetLength           = 400.; // mm                                                              
  fCuA                    = 29;
  fCuDensity              = 8.96; // g/cm3                                                              
  fpCuLambda              = 153.; // mm                                                   
  fTAXLength              = 1615.; // mm                                                               
  fTAXDistance            = 24685.;
  fbeamLength             = 102500.0; // mm           
  fzCHOD                  = 239009.0;
  fzMUV3                  = 246800.0;
  fLFV                    = 77500.;
  fLInitialFV             = 102500.;
  fzStraw[0]              = 183508.0;
  fzStraw[1]              = 194066.0;
  fzStraw[2]              = 204459.0;
  fzStraw[3]              = 218885.0;
  fxStrawChamberCentre[0] = 101.2;
  fxStrawChamberCentre[1] = 114.4;
  fxStrawChamberCentre[2] = 92.4;
  fxStrawChamberCentre[3] = 52.8;
  frMinStraw              = 60.0; 
  frMaxStraw              = 1010.0;
  fzCHODPlane             = 239009.0;
  frMinCHOD               = 120.0;
  frMaxCHOD               = 1110.0;

  // Other parameters                                                                                   

  fDBeProdProb = 0.00069;
  fDCuProdProb = fDBeProdProb*TMath::Power((29./4.),1./3.); // ACu/ABe
  fDDecayProb  = 1.;
  fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
  fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
  fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;

  // Histos

  fhReach     = nullptr;
  fhDecay     = nullptr;
  fhWeight    = nullptr;
  fgAcc       = nullptr;
  fgYield     = nullptr;
  fgGammaTot  = nullptr;
  fgTau       = nullptr;
}

void HeavyNeutrinoMassScan::InitHist() {

  BookHisto("hReach",    new TH2D("Reach", "Probability of N reaching the FV vs N mass",    fN, fInitialMass, fFinalMass, 1000, -0.1, 1.1 ));
  BookHisto("hDecay",    new TH2D("Decay", "Probability of N decaying in the FV vs N mass", fN, fInitialMass, fFinalMass, 1000, -0.1, 1.1 ));
  BookHisto("hWeight",   new TH2D("Weight", "N weight vs N mass",                           fN, fInitialMass, fFinalMass, 1000, 0., 1.E-8 ));

  fgGammaTot = new TGraph();
  fgGammaTot->SetNameTitle("GammaTot", "N total decay width vs N mass");
  BookHisto(fgGammaTot);

  fgTau = new TGraph();
  fgTau->SetNameTitle("Tau", "N lifetime vs N mass");
  BookHisto(fgTau);

  fgAcc = new TGraph();
  fgAcc->SetNameTitle("Acc", "Acceptance vs N mass");
  BookHisto(fgAcc);

  fgYield = new TGraph();
  fgYield->SetNameTitle("Yield", "Yield per POT vs N mass");
  BookHisto(fgYield);
}

void HeavyNeutrinoMassScan::Process(Int_t) {

  TRecoLAVEvent*          LAVEvent  = (TRecoLAVEvent*)          GetEvent("LAV");
  TRecoIRCEvent*          IRCEvent  = (TRecoIRCEvent*)          GetEvent("IRC");
  TRecoSACEvent*          SACEvent  = (TRecoSACEvent*)          GetEvent("SAC");

  Double_t HNLTau         = 0.;
  Double_t gammaTot       = 0.;
  Double_t NDecayProb     = 0.;
  Double_t NReachProb     = 0.;
  Double_t LReach         = 0.;
  Double_t LeptonUSquared = 0.;
  Double_t ProdFactor     = 0.;
  Double_t DecayFactor    = 0.;
  Double_t Weight         = 0.;
  Double_t DProdProb      = 0.;
  TLorentzVector mom1;
  TLorentzVector mom2;
  TVector3 point1;
  TVector3 point2;
  TVector3 momentum1;

  fIndex = 0;
  fTemp = 0.;
  fCounter = 0;

  if (GetWithMC()) {
    Event *evt = GetMCEvent();
    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);
      if (p->GetParentID() == 0) {
	if (p->GetCharge() == 1.) {
	  mom1 = p->GetInitial4Momentum();
	}
	else if (p->GetCharge() == -1.) {
	  mom2 = p->GetInitial4Momentum();
	}
      }
    }
    
    // Computation of coupling-related quantities of all HNLs (good and bad) + scan on the mass
    
    for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
      KinePart *p = evt->GetKinePart(i);      
      if (p->GetParentID() == -1 && p->GetPDGcode() == 999) {
	point1.SetXYZ(p->GetProdPos().X(), p->GetProdPos().Y(), p->GetProdPos().Z());
	point2.SetXYZ(0., 0., fLInitialFV);
	momentum1.SetXYZ(p->GetInitial4Momentum().Px(), p->GetInitial4Momentum().Py(), p->GetInitial4Momentum().Pz());
	fMN = ComputeHNLMass(p);

	if (fCounter == 0) {
	  fTemp = fMN;
	  fCounter++;
	  fMasses[fIndex] = fMN;
	}

	if (fTemp != fMN)
	  fMasses[fIndex] = fMN;
	else
	  fIndex = 0;

	fNevents[fIndex]++;
	gammaTot = GammaTot(fMN);
	HNLTau = tauN(gammaTot);
	fGammaTot[fIndex] = gammaTot;
	fTau[fIndex] = HNLTau;
	LReach = ComputeL(point1, point2, momentum1);
	ProdFactor = ComputeProd(p, fMN);
	DecayFactor = ComputeDecay(fMN);
	NReachProb = ComputeNReachProb(p, HNLTau, LReach);
	NDecayProb = ComputeNDecayProb(p, HNLTau, fLFV);
	  
	if (p->GetParticleName().Contains("e") && !p->GetParticleName().Contains("nu_tau"))
	  LeptonUSquared = fUeSquared;
	else if (p->GetParticleName().Contains("mu") && !p->GetParticleName().Contains("nu_tau"))
	  LeptonUSquared = fUmuSquared;
	else if (p->GetParticleName() == "DS->Ntau" || p->GetParticleName().Contains("nu_tau") || (p->GetParticleName().Contains("tau") && (p->GetParticleName().Contains("rho") || p->GetParticleName().Contains("pi"))))
	  LeptonUSquared = fUtauSquared;
	  
	if (p->GetProdPos().Z() < fTAXDistance)
	  DProdProb = fDBeProdProb;
	else if (p->GetProdPos().Z() >= fTAXDistance)
	  DProdProb = fDCuProdProb;
	  
	// Weight to be associated to each HNL
	  
	Weight = DProdProb*fDDecayProb*NReachProb*NDecayProb*DecayFactor*ProdFactor*LeptonUSquared;
	fSumAll[fIndex] += Weight;
	fIndex++;

	FillHisto("hReach",  fMN, NReachProb);
	FillHisto("hDecay",  fMN, NDecayProb);
	FillHisto("hWeight", fMN, Weight);
      }
    }
  }

  // Select two-track events
    
  std::vector<DownstreamTrack> Tracks = *(std::vector<DownstreamTrack>*) GetOutput("DownstreamTrackBuilder.Output");

  if (Tracks.size() == 2) {
      
    // Track features
      
    Int_t Charge1                                 = Tracks[0].GetCharge();
    Int_t Charge2                                 = Tracks[1].GetCharge();
    Double_t ChiSquare1                           = Tracks[0].GetChi2();
    Double_t ChiSquare2                           = Tracks[1].GetChi2();
    TRecoSpectrometerCandidate* SpectrometerCand1 = Tracks[0].GetSpectrometerCandidate();
    TRecoSpectrometerCandidate* SpectrometerCand2 = Tracks[1].GetSpectrometerCandidate();
    TVector3 Mom1                                 = SpectrometerCand1->GetThreeMomentumBeforeMagnet();
    TVector3 Mom2                                 = SpectrometerCand2->GetThreeMomentumBeforeMagnet();
    TVector3 TotMom                               = Mom1 + Mom2;
      
    // Compute CDA of track1,2 wrt kaon axis, and between each other
      
    fCDAcomp->SetLine1Point1(0.0, 0.0, 102000.0);     // kaon axis
    fCDAcomp->SetDir1(1.2E-3, 0., 1.);
  
    fCDAcomp->SetLine2Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);     // track1
    fCDAcomp->SetDir2(Mom1);
    fCDAcomp->ComputeVertexCDA();
  
    Double_t CDA1     = fCDAcomp->GetCDA();     // kaon axis-track1
  
    fCDAcomp->SetLine2Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand2->yAtBeforeMagnet(10.0), 10.0);     // track2
    fCDAcomp->SetDir2(Mom2);
    fCDAcomp->ComputeVertexCDA();
  
    Double_t CDA2     = fCDAcomp->GetCDA();     // kaon axis-track2
  
    fCDAcomp->SetLine1Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);     // track1
    fCDAcomp->SetDir1(Mom1);
    fCDAcomp->SetLine2Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand2->yAtBeforeMagnet(10.0), 10.0);     // track2
    fCDAcomp->SetDir2(Mom2);
    fCDAcomp->ComputeVertexCDA();
  
    Double_t CDA     = fCDAcomp->GetCDA();
    TVector3 Vertex  = fCDAcomp->GetVertex();
    Double_t Zvertex = fCDAcomp->GetVertex().z();     // track1-track2
  
    // Compute distance of two-track vertex wrt beam axis
  
    fDistcomp->SetLinePoint1(0., 0., 102000.);
    fDistcomp->SetLineDir(0., 0., 1.);
    fDistcomp->SetPoint(Vertex);
  
    fDistcomp->ComputeDistance();
  
    Double_t BeamlineDist = fDistcomp->GetDistance();

    // Track selection, CUT 1: Two tracks in Spectrometer acceptance
      
    for (Int_t i = 0; i < 4; i++) {
      Double_t x1        = SpectrometerCand1->xAt(fzStraw[i]);
      Double_t y1        = SpectrometerCand1->yAt(fzStraw[i]);
      Double_t r1        = sqrt(x1*x1 + y1*y1);
      Double_t rShifted1 = sqrt(pow(x1-fxStrawChamberCentre[i],2) + y1*y1);
      Double_t x2        = SpectrometerCand2->xAt(fzStraw[i]);
      Double_t y2        = SpectrometerCand2->yAt(fzStraw[i]);
      Double_t r2        = sqrt(x2*x2 + y2*y2);
      Double_t rShifted2 = sqrt(pow(x2-fxStrawChamberCentre[i],2) + y2*y2);
      Bool_t inAcc       = false;
	
      if ((rShifted1 > frMinStraw && r1 < frMaxStraw) && (rShifted2 > frMinStraw && r2 < frMaxStraw))
	inAcc = true;
      if (inAcc)  {
	  
	// Track selection, CUT 2: Chi2 and momentum cuts
	  
	if (ChiSquare1 <= 20. && ChiSquare2 <= 20.) {
	  if (SpectrometerCand1->GetNChambers() >= 3 && SpectrometerCand2->GetNChambers() >= 3) {
	      
	    // Track selection, CUT 3: Opposite-charged tracks
  
	    if (Charge1 + Charge2 == 0) {
		
	      // Downstream track selection, CUT 4: Extrapolation and association to CHOD
		
	      Bool_t CHODAssoc = (Tracks[0].CHODAssociationExists() && Tracks[1].CHODAssociationExists());
	      x1      = SpectrometerCand1->xAtAfterMagnet(fzCHODPlane);
	      y1      = SpectrometerCand1->yAtAfterMagnet(fzCHODPlane);
	      r1      = sqrt(x1*x1+y1*y1);
	      x2      = SpectrometerCand2->xAtAfterMagnet(fzCHODPlane);
	      y2      = SpectrometerCand2->yAtAfterMagnet(fzCHODPlane);
	      r2      = sqrt(x2*x2+y2*y2);
	      inAcc     = false;
		
	      if ((r1 > frMinCHOD && r1 < frMaxCHOD) && (r2 > frMinCHOD && r2 < frMaxCHOD))
		inAcc = true;
	      if (inAcc) {
		if (CHODAssoc) {
		    
		  // Downstream track selection, CUT 5: Extrapolation and association to LKr
		    
		  Bool_t LKrAssoc = (Tracks[0].LKrAssociationExists() && Tracks[1].LKrAssociationExists());

		  if (GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kLKr) && GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kLKr)) {
		    if (LKrAssoc) {
			
		      // Downstream track selection, CUT 6: Extrapolation and association to MUV3
    
		      if (GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand1, kMUV3) && GeometricAcceptance::GetInstance()->InAcceptance(SpectrometerCand2, kMUV3)) {
			  
			Bool_t Assoc1 = Tracks[0].MUV3AssociationExists();
			Bool_t Assoc2 = Tracks[1].MUV3AssociationExists();
			Int_t Assoc   = 0;
			  
			if (Assoc1 && !Assoc2)
			  Assoc = 1;
			else if (!Assoc1 && Assoc2)
			  Assoc = 2;
			else
			  Assoc = 0;

			if (Assoc == 1 || Assoc == 2) {
			    
			  // Compute time of MUV3 and CHOD candidates for better resolution wrt Spectrometer tracks
  
			  Double_t MUV3Time;
  
			  if (Assoc == 1)
			    MUV3Time = Tracks[0].GetMUV3Time(0);
			  else if (Assoc == 2)
			    MUV3Time = Tracks[1].GetMUV3Time(0);
			    
			  Double_t CHODTime1 = Tracks[0].GetCHODTime();
			  Double_t CHODTime2 = Tracks[1].GetCHODTime();
			    
			  // Energy cuts, CUT 7: Cut on E/p in LKr
			    
			  Double_t EoP1  = Tracks[0].GetLKrEoP();
			  Double_t EoP2  = Tracks[1].GetLKrEoP();
			  Double_t MuEoP = 0.;
			  Double_t PiEoP = 0.;
			    
			  if (Assoc == 1) {
			    MuEoP = EoP1;
			    PiEoP = EoP2;
			  }
			  else if (Assoc == 2) {
			    MuEoP = EoP2;
			    PiEoP = EoP1;
			  }
  
			  if (MuEoP < 0.2) {
			    if (PiEoP > 0.2 && PiEoP < 0.8) {
				
			      // Veto cuts, CUT 8: LAV veto
				
			      fLAVMatching->SetReferenceTime((CHODTime1 + CHODTime2) / 2);
				
			      if (GetWithMC())
				fLAVMatching->SetTimeCuts(99999, 99999);     // LAV is not time-aligned in MC
			      else
				fLAVMatching->SetTimeCuts(-10., 10.);
				
			      if (!fLAVMatching->LAVHasTimeMatching(LAVEvent)) {
				  
				// Veto cuts, CUT 9: SAV veto
				  
				fSAVMatching->SetReferenceTime((CHODTime1 + CHODTime2) / 2);
				  
				if (GetWithMC()) {
				  fSAVMatching->SetIRCTimeCuts(99999, 99999);     // SAV is not time-aligned in MC
				  fSAVMatching->SetSACTimeCuts(99999, 99999);
				} else {
				  fSAVMatching->SetIRCTimeCuts(10.0, 10.0);
				  fSAVMatching->SetSACTimeCuts(10.0, 10.0);
				}
  
				if (!fSAVMatching->SAVHasTimeMatching(IRCEvent, SACEvent)) {
				    
				  // Geometrical cuts, CUT 11: Cut on CDA of two tracks
				    
				  if (CDA < 10.) {
				      
				    // Geometrical cuts, CUT 12: Cut on two-track vertex wrt beamline

				    if (BeamlineDist > 100.) {
					
				      // Geometrical cuts, CUT 13: Cut on Z of two-track vertex

				      if (Zvertex > 102500. && Zvertex < 180000.) {
					  
					// Geometrical cuts, CUT 14: Cut on CDA of each track wrt beam axis
  
					if (CDA1 > 50. && CDA2 > 50.) {

					  // Computation of invariant mass
					    
					  // Computation of coupling-related quantities of the only good HNL in each event
					    
					  if (GetWithMC()) {
					    Event *evt = GetMCEvent();
					    for (Int_t j = 0; j < evt->GetNKineParts(); j++) {
					      KinePart *p = evt->GetKinePart(j);
					      if (p->GetParentID() == -1 && p->GetPDGcode() == 999 && p->GetEndProcessName() == "good") {
						point1.SetXYZ(p->GetProdPos().X(), p->GetProdPos().Y(), p->GetProdPos().Z());
						point2.SetXYZ(0., 0., fLInitialFV);
						momentum1.SetXYZ(p->GetInitial4Momentum().Px(), p->GetInitial4Momentum().Py(), p->GetInitial4Momentum().Pz());
						fMN = ComputeHNLMass(p);
						gammaTot = GammaTot(fMN);
						HNLTau = tauN(gammaTot);
						LReach = ComputeL(point1, point2, momentum1);
						ProdFactor = ComputeProd(p, fMN);
						DecayFactor = ComputeDecay(fMN);
						NReachProb = ComputeNReachProb(p, HNLTau, LReach);
						NDecayProb = ComputeNDecayProb(p, HNLTau, fLFV);
						  
						if (p->GetParticleName().Contains("e") && !p->GetParticleName().Contains("nu_tau"))
						  LeptonUSquared = fUeSquared;
						else if (p->GetParticleName().Contains("mu") && !p->GetParticleName().Contains("nu_tau"))
						  LeptonUSquared = fUmuSquared;
						else if (p->GetParticleName() == "DS->Ntau" || p->GetParticleName().Contains("nu_tau") || (p->GetParticleName().Contains("tau") && (p->GetParticleName().Contains("rho") || p->GetParticleName().Contains("pi"))))
						  LeptonUSquared = fUtauSquared;
						  
						if (p->GetProdPos().Z() < fTAXDistance)
						  DProdProb = fDBeProdProb;
						else if(p->GetProdPos().Z() >= fTAXDistance)
						  DProdProb = fDCuProdProb;
						  
						Weight = DProdProb*fDDecayProb*NReachProb*NDecayProb*DecayFactor*ProdFactor*LeptonUSquared;
						fSumGood[fIndex-1] += Weight;
					      }
					    }
					  }
					}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void HeavyNeutrinoMassScan::EndOfJobUser() {
  
  // Retrieve histos

  fhReach    = (TH2D*)fHisto.GetTH2("hReach");
  fhDecay    = (TH2D*)fHisto.GetTH2("hDecay");
  fhWeight   = (TH2D*)fHisto.GetTH2("hWeight");

  fhReach   ->GetXaxis()->SetTitle("N mass [GeV]");
  fhDecay   ->GetXaxis()->SetTitle("N mass [GeV]");
  fhWeight  ->GetXaxis()->SetTitle("N mass [GeV]");

  fhReach   ->GetYaxis()->SetTitle("Reach probability");
  fhDecay   ->GetYaxis()->SetTitle("Decay probability");
  fhWeight  ->GetYaxis()->SetTitle("Weight");

  // Acceptance computation

  for (Int_t i = 0; i < fN; i++) {
    fAcc[i] = fSumGood[i]/fSumAll[i];
    fProb[i] = fSumAll[i]/fNevents[i];
    fYield[i] = fAcc[i]*fProb[i];

    fgGammaTot->SetPoint(i, fMasses[i]/1000., fGammaTot[i]);
    fgTau     ->SetPoint(i, fMasses[i]/1000., fTau[i]);
    fgAcc     ->SetPoint(i, fMasses[i]/1000., fAcc[i]);
    fgYield   ->SetPoint(i, fMasses[i]/1000., fYield[i]);
  }

  // Set titles, etc.

  fgGammaTot->GetXaxis()->SetTitle("N mass [GeV]");
  fgTau     ->GetXaxis()->SetTitle("N mass [GeV]");
  fgAcc     ->GetXaxis()->SetTitle("N mass [GeV]");
  fgYield   ->GetXaxis()->SetTitle("N mass [GeV]");
  fgGammaTot->GetYaxis()->SetTitle("Total decay width [MeV]");
  fgTau     ->GetYaxis()->SetTitle("Lifetime [ns]");
  fgAcc     ->GetYaxis()->SetTitle("Acceptance");
  fgYield   ->GetYaxis()->SetTitle("Yield");
  fgGammaTot->SetLineColor(2);
  fgTau     ->SetLineColor(2);
  fgAcc     ->SetLineColor(2);
  fgYield   ->SetLineColor(2);
  fgGammaTot->SetLineWidth(3);
  fgTau     ->SetLineWidth(3);
  fgAcc     ->SetLineWidth(3);
  fgYield   ->SetLineWidth(3);
  fgGammaTot->Draw("AC");
  fgTau     ->Draw("AC");
  fgAcc     ->Draw("AC");
  fgYield   ->Draw("AC");
  fgGammaTot->Write();
  fgTau     ->Write();
  fgAcc     ->Write();
  fgYield   ->Write();

  SaveAllPlots();

  return;
}

// HNL mass

Double_t HeavyNeutrinoMassScan::ComputeHNLMass(KinePart* p) {

  Double_t MN = TMath::Sqrt(p->GetInitial4Momentum().E()*p->GetInitial4Momentum().E() - p->GetInitial4Momentum().P()*p->GetInitial4Momentum().P());

  return MN;
}

// Distance between two points

Double_t HeavyNeutrinoMassScan::ComputeL(TVector3 p1, TVector3 p2, TVector3 mom1) {
 
  TVector3 r(p1.x() + mom1.Px()/mom1.Pz()*(p2.z()-p1.z()), p1.y() + mom1.Py()/mom1.Pz()*(p2.z()-p1.z()), p2.z());
  Double_t x = r.x()-p1.x();
  Double_t y = r.y()-p1.y();
  Double_t z = r.z()-p1.z();
  Double_t L = TMath::Sqrt(x*x + y*y + z*z);

  return L;
}

// Probability of HNL decaying in FV

Double_t HeavyNeutrinoMassScan::ComputeNDecayProb(KinePart* p, Double_t tau, Double_t l) {

  Double_t cLight = 299.9; // mm/ns
  Double_t NProb  = 1. - TMath::Exp(-l/(p->GetInitial4Momentum().Beta()*p->GetInitial4Momentum().Gamma()*cLight*tau));

  return NProb;
}

// Probability of HNL reaching FV

Double_t HeavyNeutrinoMassScan::ComputeNReachProb(KinePart* p, Double_t tau, Double_t l) {

  Double_t cLight = 299.9; // mm/ns
  Double_t NProb = TMath::Exp(-l/(p->GetInitial4Momentum().Beta()*p->GetInitial4Momentum().Gamma()*cLight*tau));

  return NProb;
}

// Phasespace for 2-body HNL production mode

Double_t HeavyNeutrinoMassScan::PhaseSpace(Double_t Mass1, Double_t Mass2, Double_t Mass3) {

  Double_t phaseSpace = TMath::Power(Mass1*Mass1 - Mass2*Mass2 - Mass3*Mass3, 2) - 4.*Mass2*Mass2*Mass3*Mass3;

  return phaseSpace;
}

// Phasespace factor for 2-body HNL production mode

Double_t HeavyNeutrinoMassScan::PhaseSpaceFactor(Double_t Mass1, Double_t Mass2, Double_t Mass3) {

  Double_t factor = 0.;
  Double_t phaseSpace = PhaseSpace(Mass1, Mass2, Mass3);

  if(phaseSpace > 0.) {
    factor = (Mass1*Mass1*(Mass2*Mass2 + Mass3*Mass3) - TMath::Power(Mass2*Mass2 - Mass3*Mass3,2))*TMath::Power(phaseSpace, 0.5)/(Mass3*Mass3*TMath::Power(Mass1*Mass1 - Mass3*Mass3, 2));
  }

  return factor;
}

// Total BR for 2-body HNL production mode

Double_t HeavyNeutrinoMassScan::TwoBodyBR(Double_t Mass1, Double_t Mass2, Double_t Mass3, Int_t Dorigin, Bool_t noU) {

  Double_t brt = 0.;
  Double_t life = 0.;
  Double_t V = 0.;
  Double_t f = 0.;
  Double_t a = 0.;
  Double_t b = 0.;
  Double_t c = 0.;
  Double_t d = 0.;
  Double_t U2 = 0.;

  if (Mass1 == fMD) {
    Dorigin = 0;
  }
  else if (Mass1 == fMDS) {
    Dorigin = 1;
  }

  if (Mass1 >= (Mass2 + Mass3) && PhaseSpaceFactor(Mass1, Mass2, Mass3) > 0.) {
    if (Mass1 == fMD) {
      life = fDlife;
      V = fVcd;
      f = fD;
    }
    else if (Mass1 == fMDS) {
      life = fDSlife;
      V = fVcs;
      f = fDS;
    }
    else if (Mass1 == fMtau) {
      life = ftaulife;
      f = fPi;
      V = fVud;
    }
    else {
      cout<<"[TwoBodyBR] Unknown mother hadron"<<endl;
      _exit(1);
    }

    if (Mass3 != fMe && Mass3 != fMmu && Mass3 != fMtau && Mass3 != fMpi && Mass3 != fMrho) {
      cout<<"[TwoBodyBR] Unknown 2-body decay"<<endl;
      _exit(1);
    }

    if (noU == true)
      U2 = 1.;
    else {
      if (Mass3 == fMe)
	U2 = fUeSquared;
      else if (Mass3 == fMmu)
	U2 = fUmuSquared;
      else if (Mass3 == fMtau)
	U2 = fUtauSquared;
      else if (Mass3 == fMpi || Mass3 == fMrho)
	U2 = fUtauSquared;
    }

    if (Mass1 != fMtau) { // D,DS->Nl
      a = U2*life*fGF*fGF*f*f*V*V*Mass1*Mass2*Mass2/(8.*TMath::Pi());
      b = 1. - Mass2*Mass2/(Mass1*Mass1) + 2.*Mass3*Mass3/(Mass1*Mass1);
      c = (1. - Mass3*Mass3/(Mass1*Mass1))*Mass3*Mass3/(Mass2*Mass2);
      d = TMath::Power(1. + Mass2*Mass2/(Mass1*Mass1) - Mass3*Mass3/(Mass1*Mass1), 2.) - 4.*Mass2*Mass2/(Mass1*Mass1);
      brt = a*(b+c)*TMath::Sqrt(d);
    }
    else if (Mass1 == fMtau) { // D,DS->taunu; tau->NH (H = pi, rho)
      if ((Dorigin == 0 && PhaseSpaceFactor(fMD, fMtau, 0.) > 0.) || (Dorigin == 1 && PhaseSpaceFactor(fMDS, fMtau, 0.))) {
        if (Mass3 == fMpi) {
          a = U2*life*fGF*fGF*V*V*f*f*Mass1*Mass1*Mass1/(16.*TMath::Pi());
          b = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.) - (1. + Mass2*Mass2/(Mass1*Mass1))*Mass3*Mass3/(Mass1*Mass1);
          c = 1. - ((Mass3 - Mass2)*(Mass3 - Mass2)/(Mass1*Mass1));
          d = 1. - ((Mass3 + Mass2)*(Mass3 + Mass2)/(Mass1*Mass1));
          brt = a*b*TMath::Sqrt(c*d);
        }
        else if (Mass3 == fMrho) {
          a = U2*life*fRho*fRho*fGF*fGF*V*V*Mass1*Mass1*Mass1/(8.*TMath::Pi()*Mass3*Mass3);
          b = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.) + (1. + (Mass2*Mass2 - 2.*Mass3*Mass3)/(Mass1*Mass1))*Mass3*Mass3/(Mass1*Mass1);
          c = 1. - ((Mass3 - Mass2)*(Mass3 - Mass2)/(Mass1*Mass1));
          d = 1. - ((Mass3 + Mass2)*(Mass3 + Mass2)/(Mass1*Mass1));
          brt = a*b*TMath::Sqrt(c*d);
        }
      }
    }
  }

  return brt;
}

// Total BR for 3-body HNL production mode

Double_t HeavyNeutrinoMassScan::ThreeBodyBR(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t Mass4, Int_t Dorigin, Bool_t noU) {

  Double_t br = 0.;
  Double_t U2 = 0.;

  if (Mass1 == fMD) {
    Dorigin = 0;
  }
  else if (Mass1 == fMDS) {
    Dorigin = 1;
  }
	
  if (Mass1 >= (Mass2 + Mass3 + Mass4)) {
    if (Mass3 == fMK || Mass3 == fMK0 || Mass3 == fMpi || Mass3 == fMpi0) { // D,D0->NHl (H = pi, pi0, K, K0)
      if (Mass1 == fMD || Mass1 == fMD0) {
        Double_t ENmin = Mass2; // N at rest, K and e back to back                                     
        Double_t ENmax = (Mass1*Mass1 + Mass2*Mass2 - TMath::Power(Mass4 + Mass3, 2.))/(2.*Mass1); // N one way,K and e other way, their momenta summed equal to the N one                                      
	Double_t q2min = TMath::Power(Mass2 + Mass4, 2.); // sum of masses of lepton pair           
	Double_t q2max = TMath::Power(Mass1 - Mass3, 2.); // sum of 4momenta of lepton pair, when K at rest and N and e back to back                                              
	
	Double_t tau = 0.;
	Double_t V = 0.;
	Double_t f = 0.;
	Double_t a = 0.;
	Double_t b = 0.;
	Double_t g = 0.;

	if (noU == true)
	  U2 = 1.;
	else {
	  if (Mass4 == fMe)
	    U2 = fUeSquared;
	  else if (Mass4 == fMmu)
	    U2 = fUmuSquared;
	  else if (Mass4 == fMtau)
	    U2 = fUtauSquared;
	}

        if (Mass1 == fMD) {
          tau = fDlife;
          if (Mass3 == fMK0) {
            V = fVcs;
            f = fDK0;
            g = fgDK0;
          }
          else if (Mass3 == fMpi0) {
            V = fVcd;
            f = fDpi0;
            g = fgDpi0;
          }
          else {
            cout<<"[ThreeBodyBR] Unknown daughter hadron"<<endl;
            _exit(1);
          }
        }
        else if (Mass1 == fMD0) {
          tau = fD0life;
          if (Mass3 == fMK) {
            V = fVcs;
            f = fD0K;
            g = fgD0K;
          }
          else if (Mass3 == fMpi) {
            V = fVcd;
            f = fD0pi;
            g = fgD0pi;
          }
          else {
            cout<<"[ThreeBodyBR] Unknown daughter hadron"<<endl;
            _exit(1);
          }
        }

        a = U2*tau*V*V*fGF*fGF/(64.*TMath::Power(TMath::Pi(), 3.)*Mass1*Mass1);

	std::string function = ThreeBodyFunction(Mass1, Mass3);
        TF2* func = new TF2("func", function.c_str());

        func->SetParameter(0, f);
        func->SetParameter(1, Mass1);
        func->SetParameter(2, Mass2);
        func->SetParameter(3, Mass3);
        func->SetParameter(4, Mass4);
        func->SetParameter(5, g);

	ROOT::Math::WrappedMultiTF1 wf1(*func, 2);
	ROOT::Math::AdaptiveIntegratorMultiDim ig;
	ig.SetFunction(wf1);
	ig.SetRelTolerance(0.001);
	double xmin[] = {q2min, ENmin};
	double xmax[] = {q2max, ENmax};
	b = ig.Integral(xmin, xmax);
	br = a*b;
	
	delete func;
	func = nullptr;
      }
    }
    else if (Mass3 == fMKStar || Mass3 == fMK0Star) { // D,D0->NVl (V = K*, K0*)  
      if (Mass1 == fMD || Mass1 == fMD0) {
        Double_t ENmin = Mass2; // N at rest, K and e back to back                      
        Double_t ENmax = (Mass1*Mass1 + Mass2*Mass2 - TMath::Power(Mass4 + Mass3, 2.))/(2.*Mass1); // N one way, K and e other way, their momenta summed equal to the N one                                 
	Double_t q2min = TMath::Power(Mass2 + Mass4, 2.); // sum of masses of lepton pair  
        Double_t q2max = TMath::Power(Mass1 - Mass3, 2.); // sum of 4momenta of lepton pair, when K at rest and N and e back to back                                                         
	
	Double_t tau = 0.;
	Double_t V = 0.;
	Double_t f1 = 0.;
	Double_t f2 = 0.;
	Double_t f3 = 0.;
	Double_t f4 = 0.;
	Double_t a = 0.;
	Double_t b = 0.;
	Double_t omega2 = 0.;
	Double_t Omega2 = 0.;

	if (noU == true)
	  U2 = 1.;
	else {
	  if (Mass4 == fMe)
	    U2 = fUeSquared;
	  else if (Mass4 == fMmu)
	    U2 = fUmuSquared;
	  else if (Mass4 == fMtau)
	    U2 = fUtauSquared;
	}

        if (Mass1 == fMD) {
          tau = fDlife;
          V = fVcs;
          f1 = fVD/(Mass1 + Mass3);
          f2 = (Mass1 + Mass3)*fA1D;
          f3 = -fA2D/(Mass1 + Mass3);
          f4 = Mass3*(2.*fA0D - fA1D - fA2D) + Mass1*(fA2D - fA1D); // multiply by 1/x
        }
        else if (Mass1 == fMD0) {
          tau = fD0life;
          V = fVcs;
          f1 = fVD0/(Mass1 + Mass3);
          f2 = (Mass1 + Mass3)*fA1D0;
          f3 = -fA2D0/(Mass1 + Mass3);
          f4 = Mass3*(2.*fA0D0 - fA1D0 - fA2D0) + Mass1*(fA2D0 - fA1D0); // multiply by 1/x 
        }

        omega2 = Mass1*Mass1 - Mass3*Mass3 + Mass2*Mass2 - Mass4*Mass4; // add - 2.*Mass1*y;           
        Omega2 = Mass1*Mass1 - Mass3*Mass3; // add -x                                            
        a = U2*tau*V*V*fGF*fGF/(32.*TMath::Power(TMath::Pi(), 3.)*Mass1*Mass1);

	std::string function = ThreeBodyFunction(Mass1, Mass3);
        TF2* func = new TF2("func", function.c_str());

        func->SetParameter(0, omega2);
        func->SetParameter(1, Omega2);
        func->SetParameter(2, Mass2);
        func->SetParameter(4, Mass4);
        func->SetParameter(3, Mass3);
        func->SetParameter(5, f1);
        func->SetParameter(6, f2);
        func->SetParameter(7, f3);
        func->SetParameter(8, f4);
        func->SetParameter(9, Mass1);

	ROOT::Math::WrappedMultiTF1 wf1(*func, 2);
	ROOT::Math::AdaptiveIntegratorMultiDim ig;
	ig.SetFunction(wf1);
	ig.SetRelTolerance(0.001);
	double xmin[] = {q2min, ENmin};
	double xmax[] = {q2max, ENmax};
	b = ig.Integral(xmin, xmax);
	br = a*b;
	
	delete func;
	func = nullptr;
      }
    }
    else if (Mass1 == fMtau) {
      if ((Dorigin == 0 && PhaseSpaceFactor(fMD, fMtau, 0.) > 0.) || (Dorigin == 1 && PhaseSpaceFactor(fMDS, fMtau, 0.))) {
        Double_t b = 0.;
        Double_t ENmin = 0.;
        Double_t ENmax = 0.;
        Double_t life = ftaulife;

        if (Mass3 == 0.1) { // D,DS->taunu_tau; tau->Nlnu_tau
	  std::string function = ThreeBodyFunction(Mass1, Mass3);

	  if (noU == true)
	    U2 = 1.;
	  else {
	    if (Mass4 == fMe)
	      U2 = fUeSquared;
	    else if (Mass4 == fMmu)
	      U2 = fUmuSquared;
	    else if (Mass4 == fMtau)
	      U2 = fUtauSquared;
	  }

          Mass3 = 0.;
          ENmin = Mass2; // N at rest, l and nu back to back                           
          ENmax = (Mass1*Mass1 + Mass2*Mass2 - TMath::Power(Mass4 + Mass3, 2.))/(2.*Mass1); // N one way, l and nu other way, their momenta summed equal to the N one                                           
	  
	  TF1 func("func", function.c_str());

          func.SetParameter(0, life);
          func.SetParameter(1, Mass1);
          func.SetParameter(2, Mass2);
          func.SetParameter(3, fGF);
          func.SetParameter(4, Mass4);

	  ROOT::Math::WrappedTF1 wf1(func);
	  ROOT::Math::GaussLegendreIntegrator ig;
	  ig.SetFunction(wf1);
	  b = ig.Integral(ENmin, ENmax);
	  br = U2*b;
	}
        else if (Mass3 == 0.01) { // D,DS->taunu_tau; tau->Nlnu_l
	  std::string function = ThreeBodyFunction(Mass1, Mass3);

	  if (noU == true)
	    U2 = 1.;
	  else 
	    U2 = fUtauSquared;

          Mass3 = 0.;
          ENmin = Mass2; // N at rest, l and nu back to back                                       
          ENmax = (Mass1*Mass1 + Mass2*Mass2 - TMath::Power(Mass4 + Mass3, 2.))/(2.*Mass1); // N one way, l and nu other way, their momenta summed equal to the N one                                          
	  
	  TF1 func("func", function.c_str());

          func.SetParameter(0, life);
          func.SetParameter(1, Mass1);
          func.SetParameter(2, Mass2);
          func.SetParameter(3, fGF);
          func.SetParameter(4, Mass4);

	  ROOT::Math::WrappedTF1 wf1(func);
	  ROOT::Math::GaussLegendreIntegrator ig;
	  ig.SetFunction(wf1);
	  b = ig.Integral(ENmin, ENmax);
	  br = U2*b;
        }
        else {
          cout<<"[ThreeBodyBR] Unknown neutrino type in N 3-body production mode"<<endl;
          _exit(1);
        }
      }
    }
    else {
      cout<<"[ThreeBodyBR] Unknown N 3-body production mode"<<endl;
      _exit(1);
    }
  }

  return br;
}

// Create string function for 3-body total BR of HNL production mode

std::string HeavyNeutrinoMassScan::ThreeBodyFunction(Double_t Mass1, Double_t Mass3) {

  std::string function = "";

  if (Mass1 == fMD || Mass1 == fMD0) {
    if (Mass3 == fMK || Mass3 == fMK0 || Mass3 == fMpi || Mass3 == fMpi0) { // D,D0->NHl
      function = "([5]*[5]*(x*([2]*[2] + [4]*[4]) - TMath::Power([2]*[2] - [4]*[4], 2.)) + 2.*[5]*[0]*([2]*[2]*(2.*[1]*[1] - 2.*[3]*[3] -4.*y*[1] - [4]*[4] + [2]*[2] + x) + [4]*[4]*(4.*y*[1] + [4]*[4] - [2]*[2] - x)) + [0]*[0]*((4.*y*[1] + [4]*[4] - [2]*[2] - x)*(2.*[1]*[1] - 2.*[3]*[3] - 4.*y*[1] - [4]*[4] + [2]*[2] + x) - (2.*[1]*[1] + 2.*[3]*[3] - x)*(x - [2]*[2] - [4]*[4])))";
    }
    else if (Mass3 == fMKStar || Mass3 == fMK0Star) { // D,D0->NVl
      function = "(([6]*[6]/2.)*(x - [2]*[2] - [4]*[4] + ([0] - 2.*[9]*y)*([1] - x - ([0] - 2.*[9]*y))/([3]*[3])) + (([7]+[8]*1./x)*([7]+[8]*1./x)/2.)*([2]*[2] + [4]*[4])*(x - [2]*[2] + [4]*[4])*(([1] - x)*([1] - x)/(4.*[3]*[3]) - x) + 2.*[7]*[7]*[3]*[3]*(([1] - x)*([1] - x)/(4.*[3]*[3]) - x)*([2]*[2] + [4]*[4] - x + ([0] - 2.*[9]*y)*([1] - x - ([0] - 2.*[9]*y))/([3]*[3])) + 2.*[7]*([7]+[8]*1./x)*([2]*[2]*([0] - 2.*[9]*y) + ([1] - x - ([0] - 2.*[9]*y))*[4]*[4])*(([1] - x)*([1] - x)/(4.*[3]*[3]) - x) + 2.*[5]*[6]*(x*(2.*([0] - 2.*[9]*y) - [1] + x) + ([1] - x)*([2]*[2] - [4]*[4])) + ([6]*([7]+[8]*1./x)/2.)*(([0] - 2.*[9]*y)*([1] - x)/([3]*[3])*([2]*[2] - [4]*[4]) + ([1] - x)*([1] - x)*[4]*[4]/([3]*[3]) + 2.*TMath::Power([2]*[2] - [4]*[4], 2.) - 2.*x*([2]*[2] + [4]*[4])) + [6]*[7]*(([1] - x)*([0] - 2.*[9]*y)*([1] - x - ([0] - 2.*[9]*y))/([3]*[3]) + 2.*([0] - 2.*[9]*y)*([4]*[4] - [2]*[2]) + ([1] - x)*([2]*[2] - [4]*[4] - x)) + [5]*[5]*(([1] - x)*([1] - x)*(x - [2]*[2] + [4]*[4]) - 2.*[3]*[3]*(x*x - TMath::Power([2]*[2] - [4]*[4], 2.)) + 2.*([0] - 2.*[9]*y)*([1] - x)*([2]*[2] - x - [4]*[4]) + 2.*([0] - 2.*[9]*y)*([0] - 2.*[9]*y)*x))";
    }
  }
  else if (Mass1 == fMtau) { // D,DS->taunu_tau; tau->Nlnu
    if (Mass3 == 0.1) { //nu_tau                                                                      
      function = "([0]*[3]*[3]*[1]*[1]*x/(2.*TMath::Power(TMath::Pi(), 3.)))*(1. + ([2]*[2] - [4]*[4])/([1]*[1]) - 2.*x/[1])*(1. - [4]*[4]/([1]*[1] + [2]*[2] - 2.*x*[1]))*(TMath::Sqrt(x*x - [2]*[2]))";
    }
    else if (Mass3 == 0.01) { //nu_e or nu_mu                                                       
      function = "([0]*[3]*[3]*[1]*[1]/(4.*TMath::Power(TMath::Pi(), 3.)))*(TMath::Power((1. - [4]*[4]/([1]*[1] + [2]*[2] - 2.*x*[1])), 2.)*TMath::Sqrt(x*x - [2]*[2]))*(([1] - x)*(1. - ([2]*[2] + [4]*[4])/([1]*[1])) - (1. - [4]*[4]/([1]*[1] + [2]*[2] - 2.*x*[1]))*(([1] - x)*([1] - x)/[1] + (x*x - [2]*[2])/(3.*[1])))";
    }
  }

  return function;
}

// Decay width for 2-body HNL decay mode 

Double_t HeavyNeutrinoMassScan::Gamma2(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t form, Bool_t noU) {

  Double_t gamma_2 = 0.;
  Double_t V = 0.;
  Double_t a = 0.;
  Double_t b = 0.;
  Double_t c = 0.;
  Double_t d = 0.;
  Double_t f = 0.;
  Double_t g = 0.;
  Double_t U2 = 0.;

  if (noU == true)
    U2 = 1.;
  else {
    if (Mass3 == fMe)
      U2 = fUeSquared;
    else if (Mass3 == fMmu)
      U2 = fUmuSquared;
    else if (Mass3 == fMtau)
      U2 = fUtauSquared;
  }

  if (Mass1 >= (Mass2 + Mass3)) {
    if (Mass2 == fMpi || Mass2 == fMK) {
      if (Mass2 == fMpi)
	V = fVud;
      else if (Mass2 == fMK)
	V = fVus;

      a = (U2*fGF*fGF*V*V*form*form*Mass1*Mass1*Mass1)/(16.*TMath::Pi());
      b = TMath::Power(1. - Mass3*Mass3/(Mass1*Mass1), 2.);
      c = (1. + Mass3*Mass3/(Mass1*Mass1))*(Mass2*Mass2/(Mass1*Mass1));
      d = 1. - (Mass2 - Mass3)*(Mass2 - Mass3)/(Mass1*Mass1);
      f = 1. - (Mass2 + Mass3)*(Mass2 + Mass3)/(Mass1*Mass1);
      gamma_2 = a*(b - c)*TMath::Sqrt(d*f);
    }
    else if (Mass2 == fMrho) {

      V = fVud;
      a = (U2*fGF*fGF*form*form*V*V*Mass1*Mass1*Mass1)/(8.*TMath::Pi()*Mass2*Mass2);
      b = (1. - TMath::Power(Mass2/Mass1 - Mass3/Mass1, 2.))*(1. - TMath::Power(Mass2/Mass1 + Mass3/Mass1, 2.));
      c = (1. + (Mass3*Mass3)/(Mass1*Mass1))*(Mass2*Mass2/(Mass1*Mass1));
      d = 2*Mass2*Mass2*Mass2*Mass2/(Mass1*Mass1*Mass1*Mass1);
      f = TMath::Power(1. - (Mass3*Mass3)/(Mass1*Mass1), 2.);
      g = c - d + f;
      gamma_2 = a*TMath::Sqrt(b)*g;
    }
    else if (Mass2 == fMpi0) {
      a = (U2*fGF*fGF*form*form*Mass1*Mass1*Mass1)/(32.*TMath::Pi());
      b = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.);
      gamma_2 = a*b;
    }
    else if (Mass2 == fMeta || Mass2 == fMetaprime) {
      a = (U2*fGF*fGF*form*form*Mass1*Mass1*Mass1)/(32.*TMath::Pi());
      b = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.);
      gamma_2 = a*b;
    }
    else if (Mass2 == fMrho0) {
      a = (U2*form*form*fGF*fGF*Mass1*Mass1*Mass1)/(16.*TMath::Pi()*Mass2*Mass2);
      b = 1. + 2.*Mass2*Mass2/(Mass1*Mass1);
      c = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.);
      gamma_2 = a*b*c;
    }
    else {
      cout<<"[Gamma2] Unknown N two-body decay mode"<<endl;
      _exit(1);
    }
  }

  return gamma_2;
}

// Decay width for 3-body HNL decay mode 

Double_t HeavyNeutrinoMassScan::GammaLeptonNu3(Double_t Mass1, Double_t Mass2, Double_t Mass3, Bool_t noU) {

  Double_t r = 0.;
  Double_t a = 0.;
  Double_t b = 0.;
  Double_t f = 0.;
  Double_t U2 = 0.;
  Double_t gamma_l_l_nu = 0.;

  if (Mass1 >= (Mass2 + Mass3)) {
    if (Mass2 == Mass3 && Mass2 != 0.) {
      if (noU == true)
	U2 = 1.;
      else 
	U2 = fUSquared;

      r = 4*Mass3*Mass3/(Mass1*Mass1);
      a = TMath::Power(1-r,0.5)*(1./3.-7*r/6.-r*r/24.-r*r*r/16.);
      b = r*r*(1-r*r/16.)*TMath::ATanH(TMath::Power(1-r,0.5));
      f = a+b;
      gamma_l_l_nu = U2*fGF*fGF*TMath::Power(Mass1,5)*f/(64.*TMath::Power(TMath::Pi(),3));
    }
    else if (Mass2 == Mass3 && Mass2 == 0.) {
      if (noU == true)
	U2 = 1.;
      else 
	U2 = fUSquared;

      a = U2*fGF*fGF*TMath::Power(Mass1, 5.)/(192.*TMath::Power(TMath::Pi(), 3.));
      gamma_l_l_nu = a;
    }
    else if (Mass2 != Mass3) {
      if (Mass2 > Mass3) {
        r = Mass2/Mass1;
      }
      else if (Mass2 < Mass3) {
        r = Mass3/Mass1;
      }
      else {
        cout<<"[GammaLeptonNu3] N 3-body decay mode should have equal masses"<<endl;
        _exit(1);
      }

      if (noU == true)
	U2 = 1.;
      else {
	if (Mass2 == fMe)
	  U2 = fUeSquared;
	else if (Mass2 == fMmu)
	  U2 = fUmuSquared;
	else if (Mass2 == fMtau)
	  U2 = fUtauSquared;
      }

      a = (U2*fGF*fGF*TMath::Power(Mass1, 5))/(192*TMath::Power(TMath::Pi(), 3));
      b = 1 - 8.*r*r + 8.*TMath::Power(r, 6) - TMath::Power(r, 8) - 12.*TMath::Power(r, 4)*TMath::Log(r*r);
      gamma_l_l_nu = a*b;
    }
  }

  return gamma_l_l_nu;
}

// Total HNL decay width                                                                              

Double_t HeavyNeutrinoMassScan::GammaTot(Double_t MN) {

  Double_t gammaTot = 0.;

  if (MN < 2*fMe) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false);
  }
  else if (MN >= 2*fMe && MN < (fMe+fMmu)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false);
  }
  else if (MN >= (fMe+fMmu) && MN < (fMpi0)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false);
  }
  else if (MN >= (fMpi0) && MN < (fMe+fMpi)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false);
  }
  else if (MN >= (fMe+fMpi) && MN < 2*fMmu) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false);
  }
  else if (MN >= 2*fMmu && MN < (fMmu+fMpi)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false);
  }
  else if (MN >= (fMmu+fMpi) && MN < (fMK+fMe)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false);
  }
  else if (MN >= (fMK+fMe) && MN < (fMeta)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false);
  }
  else if (MN >= (fMeta) && MN < (fMK+fMmu)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false);
  }
  else if (MN >= (fMK+fMmu) && MN < (fMrho0)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false) + Gamma2(MN, fMK, fMmu, fK, false);
  }
  else if (MN >= (fMrho0) && MN < (fMrho+fMe)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false) + Gamma2(MN, fMK, fMmu, fK, false) + Gamma2(MN, fMrho0, 0., fRho, false);
  }
  else if (MN >= (fMrho+fMe) && MN < (fMrho+fMmu)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false) + Gamma2(MN, fMK, fMmu, fK, false) + Gamma2(MN, fMrho0, 0., fRho, false) + Gamma2(MN, fMrho, fMe, fRho, false);
  }
  else if (MN >= (fMmu+fMrho) && MN < (fMetaprime)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false) + Gamma2(MN, fMK, fMmu, fK, false) + Gamma2(MN, fMrho0, 0., fRho, false) + Gamma2(MN, fMrho, fMe, fRho, false) + Gamma2(MN, fMrho, fMmu, fRho, false);
  }
  else if (MN >= (fMetaprime) && MN < (fMe+fMtau)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false) + Gamma2(MN, fMK, fMmu, fK, false) + Gamma2(MN, fMrho0, 0., fRho, false) + Gamma2(MN, fMrho, fMe, fRho, false) + Gamma2(MN, fMrho, fMmu, fRho, false) + Gamma2(MN, fMetaprime, 0., fEtaprime, false);
  }
  else if (MN >= (fMe+fMtau) && MN < (fMmu+fMtau)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false) + Gamma2(MN, fMK, fMmu, fK, false) + Gamma2(MN, fMrho0, 0., fRho, false) + Gamma2(MN, fMrho, fMe, fRho, false) + Gamma2(MN, fMrho, fMmu, fRho, false) + Gamma2(MN, fMetaprime, 0., fEtaprime, false) + GammaLeptonNu3(MN, fMe, fMtau, false);
  }
  else if (MN >= (fMmu+fMtau) && MN < (fMpi+fMtau)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false) + Gamma2(MN, fMK, fMmu, fK, false) + Gamma2(MN, fMrho0, 0., fRho, false) + Gamma2(MN, fMrho, fMe, fRho, false) + Gamma2(MN, fMrho, fMmu, fRho, false) + Gamma2(MN, fMetaprime, 0., fEtaprime, false) + GammaLeptonNu3(MN, fMe, fMtau, false) + GammaLeptonNu3(MN, fMmu, fMtau, false);
  }
  else if (MN >= (fMpi+fMtau) && MN < (fMK+fMtau)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false) + Gamma2(MN, fMK, fMmu, fK, false) + Gamma2(MN, fMrho0, 0., fRho, false) + Gamma2(MN, fMrho, fMe, fRho, false) + Gamma2(MN, fMrho, fMmu, fRho, false) + Gamma2(MN, fMetaprime, 0., fEtaprime, false) + GammaLeptonNu3(MN, fMe, fMtau, false) + GammaLeptonNu3(MN, fMmu, fMtau, false) + Gamma2(MN, fMpi, fMtau, fPi, false);
  }
  else if (MN >= (fMK+fMtau) && MN < (fMrho+fMtau)) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false) + Gamma2(MN, fMK, fMmu, fK, false) + Gamma2(MN, fMrho0, 0., fRho, false) + Gamma2(MN, fMrho, fMe, fRho, false) + Gamma2(MN, fMrho, fMmu, fRho, false) + Gamma2(MN, fMetaprime, 0., fEtaprime, false) + GammaLeptonNu3(MN, fMe, fMtau, false) + GammaLeptonNu3(MN, fMmu, fMtau, false) + Gamma2(MN, fMpi, fMtau, fPi, false) + Gamma2(MN, fMK, fMtau, fK, false);
  }
  else if (MN >= (fMrho+fMtau) && MN < 2*fMtau) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false) + Gamma2(MN, fMK, fMmu, fK, false) + Gamma2(MN, fMrho0, 0., fRho, false) + Gamma2(MN, fMrho, fMe, fRho, false) + Gamma2(MN, fMrho, fMmu, fRho, false) + Gamma2(MN, fMetaprime, 0., fEtaprime, false) + GammaLeptonNu3(MN, fMe, fMtau, false) + GammaLeptonNu3(MN, fMmu, fMtau, false) + Gamma2(MN, fMpi, fMtau, fPi, false) + Gamma2(MN, fMK, fMtau, fK, false) + Gamma2(MN, fMrho, fMtau, fRho, false);
  }
  else if (MN >= 2*fMtau) {
    gammaTot = GammaLeptonNu3(MN, 0., 0., false) + GammaLeptonNu3(MN, fMe, fMe, false) + GammaLeptonNu3(MN, fMe, fMmu, false) + Gamma2(MN, fMpi0, 0., fPi, false) + Gamma2(MN, fMpi, fMe, fPi, false) + GammaLeptonNu3(MN, fMmu, fMmu, false) + Gamma2(MN, fMpi, fMmu, fPi, false) + Gamma2(MN, fMK, fMe, fK, false) + Gamma2(MN, fMeta, 0., fEta, false) + Gamma2(MN, fMK, fMmu, fK, false) + Gamma2(MN, fMrho0, 0., fRho, false) + Gamma2(MN, fMrho, fMe, fRho, false) + Gamma2(MN, fMrho, fMmu, fRho, false) + Gamma2(MN, fMetaprime, 0., fEtaprime, false) + GammaLeptonNu3(MN, fMe, fMtau, false) + GammaLeptonNu3(MN, fMmu, fMtau, false) + Gamma2(MN, fMpi, fMtau, fPi, false) + Gamma2(MN, fMrho, fMtau, fRho, false) + GammaLeptonNu3(MN, fMtau, fMtau, false);
  }

  return gammaTot;
}

// HNL lifetime

Double_t HeavyNeutrinoMassScan::tauN(Double_t MN) {

  Double_t gammaN = GammaTot(MN);

  return fhc/(gammaN*fcLight);
}

// Lambda function                                                                                   

Double_t HeavyNeutrinoMassScan::lambda(Double_t a, Double_t b, Double_t c) {

  return a*a + b*b + c*c - 2.*a*b - 2.*a*c - 2.*b*c;
}

// Factor related to HNL production

Double_t HeavyNeutrinoMassScan::ComputeProd(KinePart* p, Double_t MN) {

  Double_t br = 0.;
  Double_t Dorigin = -1;

  if (p->GetParticleName().Contains("DS")) 
    Dorigin = 1;
  else if (p->GetParticleName().Contains("D0"))
    Dorigin = -1;
  else
    Dorigin = 0;

  Double_t BR2De     = TwoBodyBR(fMD,   MN, fMe,   Dorigin, true);
  Double_t BR2Dmu    = TwoBodyBR(fMD,   MN, fMmu,  Dorigin, true);
  Double_t BR2DSe    = TwoBodyBR(fMDS,  MN, fMe,   Dorigin, true);
  Double_t BR2DSmu   = TwoBodyBR(fMDS,  MN, fMmu,  Dorigin, true);
  Double_t BR2DStau  = TwoBodyBR(fMDS,  MN, fMtau, Dorigin, true);
  Double_t BR2taupi  = TwoBodyBR(fMtau, MN, fMpi,  Dorigin, true);
  Double_t BR2taurho = TwoBodyBR(fMtau, MN, fMrho, Dorigin, true);

  Double_t BR3DK0e        = ThreeBodyBR(fMD,   MN, fMK0,     fMe,  Dorigin, true);
  Double_t BR3Dpi0e       = ThreeBodyBR(fMD,   MN, fMpi0,    fMe,  Dorigin, true);
  Double_t BR3D0Ke        = ThreeBodyBR(fMD0,  MN, fMK,      fMe,  Dorigin, true);
  Double_t BR3D0pie       = ThreeBodyBR(fMD0,  MN, fMpi,     fMe,  Dorigin, true);
  Double_t BR3DK0mu       = ThreeBodyBR(fMD,   MN, fMK0,     fMmu, Dorigin, true);
  Double_t BR3Dpi0mu      = ThreeBodyBR(fMD,   MN, fMpi0,    fMmu, Dorigin, true);
  Double_t BR3D0Kmu       = ThreeBodyBR(fMD0,  MN, fMK,      fMmu, Dorigin, true);
  Double_t BR3D0pimu      = ThreeBodyBR(fMD0,  MN, fMpi,     fMmu, Dorigin, true);
  Double_t BR3DK0Stare    = ThreeBodyBR(fMD,   MN, fMK0Star, fMe,  Dorigin, true);
  Double_t BR3D0KStare    = ThreeBodyBR(fMD0,  MN, fMKStar,  fMe,  Dorigin, true);
  Double_t BR3DK0Starmu   = ThreeBodyBR(fMD,   MN, fMK0Star, fMmu, Dorigin, true);
  Double_t BR3D0KStarmu   = ThreeBodyBR(fMD0,  MN, fMKStar,  fMmu, Dorigin, true);
  Double_t BR3tauenu_tau  = ThreeBodyBR(fMtau, MN, 0.1,      fMe,  Dorigin, true);
  Double_t BR3tauenu_e    = ThreeBodyBR(fMtau, MN, 0.01,     fMe,  Dorigin, true);
  Double_t BR3taumunu_tau = ThreeBodyBR(fMtau, MN, 0.1,      fMmu, Dorigin, true);
  Double_t BR3taumunu_mu  = ThreeBodyBR(fMtau, MN, 0.01,     fMmu, Dorigin, true);
  Double_t BR2Dtot = BR2De + BR2Dmu + BR2taupi + BR2taurho;
  Double_t BR2DStot = BR2DSe + BR2DSmu + BR2DStau + BR2taupi + BR2taurho;
  Double_t BR3Dtot = BR3DK0e + BR3DK0mu + BR3Dpi0e + BR3Dpi0mu + BR3DK0Stare + BR3DK0Starmu + BR3tauenu_tau + BR3taumunu_tau + BR3tauenu_e + BR3taumunu_mu;
  Double_t BR3DStot = BR3tauenu_tau + BR3taumunu_tau + BR3tauenu_e + BR3taumunu_mu;
  Double_t BR3D0tot = BR3D0Ke + BR3D0Kmu + BR3D0pie + BR3D0pimu + BR3D0KStare + BR3D0KStarmu;
  
  if (p->GetParticleName().Contains("DS")) {
    if ((p->GetParticleName().Contains("taunu") && !(p->GetParticleName().Contains("nu_"))) || !(p->GetParticleName().Contains("taunu"))) // DS->Nl or DS->taunu; tau->NH
      br = BR2DStot;
    else // DS->taunu; tau->Nlnu
      br = BR3DStot;
  }
  else if (p->GetParticleName().Contains("D0")) // D0->HNl
    br = BR3D0tot;
  else {
    if ((p->GetParticleName().Contains("taunu") && !(p->GetParticleName().Contains("nu_"))) || !(p->GetParticleName().Contains("taunu"))) // D->Nl or D->taunu; tau->NH 
      br = BR2Dtot;
    else // D->HNl or D->taunu; tau->Nlnu
      br = BR3Dtot;
  }

  return br;
}

// Factor related to HNL decay

Double_t HeavyNeutrinoMassScan::ComputeDecay(Double_t MN) {

  Double_t br = 0.;
  Double_t Gamma2pimu = Gamma2(MN, fMpi, fMmu, fPi, false);
  Double_t gammaTot = GammaTot(MN);

  if (Gamma2pimu > 0. && gammaTot > 0.) 
    br = Gamma2pimu/gammaTot;

  return br;
}

HeavyNeutrinoMassScan::~HeavyNeutrinoMassScan() {
  
  delete fCDAcomp;
  delete fDistcomp;
  delete fLAVMatching;
  delete fSAVMatching;

  delete fhReach;
  delete fhDecay;
  delete fhWeight;
  delete fgAcc;
  delete fgYield;
  delete fgGammaTot;
  delete fgTau;

  fhReach    = nullptr;
  fhDecay    = nullptr;
  fhWeight   = nullptr;
  fgAcc      = nullptr;
  fgYield    = nullptr;
  fgGammaTot = nullptr;
  fgTau      = nullptr;
}