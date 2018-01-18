#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cmath>
#include <typeinfo>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include "HeavyNeutrinoPiMuSelectionScan.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "BeamParameters.hh"
#include "DownstreamTrack.hh"
#include "GeometricAcceptance.hh"
#include "TriggerConditions.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

#define TRIGGER_L0_PHYSICS_TYPE 1
#define MCTriggerMask 0xFF

/// \class HeavyNeutrinoPiMuSelectionScan

HeavyNeutrinoPiMuSelectionScan::HeavyNeutrinoPiMuSelectionScan(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrinoPiMuSelectionScan") {

  RequestAllMCTrees();
  RequestAllRecoTrees();
  RequestL0Data();
  
  fCDAcomp     = new TwoLinesCDA();
  fDistcomp    = new PointLineDistance();
  fLAVMatching = new LAVMatching();
  fSAVMatching = new SAVMatching();
  fMN = 0.;

  for (Int_t i = 0; i < fN; i++) {
    fSumGood[i] = 0.;
    fSumAll[i] = 0.;
    fNevents[i] = 0;
    fCouplings[i] = 0.;
    fAcc[i] = 0.;
    fGammaTot[i] = 0.;
    fTau[i] = 0.;
    fProb[i] = 0.;
    fYield[i] = 0.;
  }
}

void HeavyNeutrinoPiMuSelectionScan::InitHist() {

  BookHisto("hReach",    new TH2D("Reach", "Probability of HNL reaching the FV vs coupling",    fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1 ));
  BookHisto("hDecay",    new TH2D("Decay", "Probability of HNL decaying in the FV vs coupling", fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1 ));
  BookHisto("hWeight",   new TH2D("Weight", "HNL weight vs coupling",                           fN, fCouplingStart, fCouplingStop, 1000, 0., 1.E-8 ));

  fgGammaTot = new TGraph();
  fgGammaTot->SetNameTitle("GammaTot", "N total decay width vs coupling");
  BookHisto(fgGammaTot);

  fgTau = new TGraph();
  fgTau->SetNameTitle("Tau", "N lifetime vs coupling");
  BookHisto(fgTau);

  fgAcc = new TGraph();
  fgAcc->SetNameTitle("Acc", "Acceptance vs coupling");
  BookHisto(fgAcc);

  fgYield = new TGraph();
  fgYield->SetNameTitle("Yield", "Yield per POT vs coupling");
  BookHisto(fgYield);
}

void HeavyNeutrinoPiMuSelectionScan::Process(Int_t) {

  TRecoLKrEvent* LKrEvent  = (TRecoLKrEvent*)GetEvent("LKr");
  TRecoLAVEvent* LAVEvent  = (TRecoLAVEvent*)GetEvent("LAV");
  TRecoIRCEvent* IRCEvent  = (TRecoIRCEvent*)GetEvent("IRC");
  TRecoSACEvent* SACEvent  = (TRecoSACEvent*)GetEvent("SAC");

  // Plots of KineParts from MC

  Double_t zStraw[4] = {183508.0, 194066.0, 204459.0, 218885.0};
  Double_t xStrawChamberCentre[4] = {101.2, 114.4, 92.4, 52.8};
  Double_t rMinStraw = 60.0;
  Double_t rMaxStraw = 1010.0;
  Double_t rMinCHOD   = 120.0;
  Double_t rMaxCHOD   = 1110.0;
  Double_t zCHOD = 239009.0;
  Double_t zMUV3 = 246800.0;
  TLorentzVector mom1(0., 0., 0., 0.);
  TLorentzVector mom2(0., 0., 0., 0.);
  Double_t mass1 = 0.;
  Double_t mass2 = 0.;
  Double_t HNLTau = 0.;
  Double_t gammaTot = 0.;
  Double_t NDecayProb = 0.;
  Double_t NReachProb = 0.;
  Double_t LFV = 77500.;
  Double_t LReach = 0.;
  Double_t LInitialFV = 102500.;
  Double_t USquared = 0.;
  Double_t UeSquared = 0.;
  Double_t UmuSquared = 0.;
  Double_t UtauSquared = 0.;
  Double_t LeptonUSquared = 0.;
  Double_t BR = 0.;
  Double_t Weight = 0.;
  TVector3 point1(0., 0., 0.);
  TVector3 point2(0., 0., 0.);
  TVector3 momentum1;
  Double_t DBeProdProb = 0.00069;
  Double_t DCuProdProb = DBeProdProb*TMath::Power((29./4.),1./3.); // ACu/ABe
  Double_t DDecayProb = 1.;
  Double_t NPiMu = 0.;
  Double_t DProdProb = 0.;
  Double_t TAXDistance = 24685.;

  // Scan on the coupling

  for(Double_t fCoupling = fCouplingStart; fCoupling < fCouplingStop-fCouplingStep; fCoupling += fCouplingStep) {

    USquared = TMath::Power(10, fCoupling);
    UeSquared = USquared/20.8;
    UmuSquared = 16.*UeSquared;
    UtauSquared = 3.8*UeSquared;
    fIndex = round((std::abs(fCouplingStart-fCoupling))/fCouplingStep);
    fCouplings[fIndex] = fCoupling;

    if (GetWithMC()) {
      Event *evt = GetMCEvent();
      for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
	KinePart *p = evt->GetKinePart(i);
	if (p->GetParentID() == 0) {
	  if (p->GetCharge() == 1.) {
	    mom1 = p->GetInitial4Momentum();
	    mass1 = TMath::Sqrt(mom1.E()*mom1.E() - mom1.P()*mom1.P());
	  }
	  else if (p->GetCharge() == -1.) {
	    mom2 = p->GetInitial4Momentum();
	    mass2 = TMath::Sqrt(mom2.E()*mom2.E() - mom2.P()*mom2.P());
	  }
	}
      }

      // Computation of coupling-related quantities of all HNLs in the MC event

      for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
	KinePart *p = evt->GetKinePart(i);      
	if (p->GetParentID() == -1 && p->GetPDGcode() == 999) {
	  fNevents[fIndex]++;
	  point1.SetXYZ(p->GetProdPos().X(), p->GetProdPos().Y(), p->GetProdPos().Z());
	  point2.SetXYZ(0., 0., LInitialFV);
	  momentum1.SetXYZ(p->GetInitial4Momentum().Px(), p->GetInitial4Momentum().Py(), p->GetInitial4Momentum().Pz());
	  fMN = ComputeHNLMass(p);
	  gammaTot = GammaTot(mass1, mass2, fMN, USquared);
	  HNLTau = Tau(gammaTot);
	  fGammaTot[fIndex] = gammaTot;
	  fTau[fIndex] = HNLTau;
	  LReach = ComputeL(point1, point2, momentum1);
	  BR = ComputeBR(p, fMN);
	  NReachProb = ComputeNReachProb(p, HNLTau, LReach);
	  NDecayProb = ComputeNDecayProb(p, HNLTau, LFV);
	  NPiMu = GammaPiMu(fMN)*UmuSquared/(GammaPiE(fMN)*UeSquared + GammaPiMu(fMN)*UmuSquared);
	  
	  if (p->GetParticleName() == "HNLD+e+" || p->GetParticleName() == "HNLD-e-" || p->GetParticleName() == "HNLDS+e+" || p->GetParticleName() == "HNLDS-e-") 
	    LeptonUSquared = UeSquared;

	  else if (p->GetParticleName() == "HNLD+mu+" || p->GetParticleName() == "HNLD-mu-" || p->GetParticleName() == "HNLDS+mu+" || p->GetParticleName() == "HNLDS-mu-") 
	    LeptonUSquared = UmuSquared;

	  if (p->GetProdPos().Z() < TAXDistance)
	    DProdProb = DBeProdProb;
	  else if (p->GetProdPos().Z() >= TAXDistance)
	    DProdProb = DCuProdProb;

	  Weight = DProdProb*DDecayProb*NReachProb*NDecayProb*NPiMu*BR*LeptonUSquared;

	  if (Weight >= 0. && Weight <= 10.)
	    fSumAll[fIndex] += Weight;

	  FillHisto("hReach",    fCoupling, NReachProb);
	  FillHisto("hDecay",    fCoupling, NDecayProb);
	  FillHisto("hWeight",   fCoupling, Weight);
	}
      }
    }

    // L0 data                                                                       

    Int_t RunNumber = GetWithMC() ? 0 : GetEventHeader()->GetRunID();   
    L0TPData *L0TPData = GetL0Data();
    UChar_t L0DataType = GetWithMC() ? 0x11 : L0TPData->GetDataType();
    UInt_t L0TriggerFlags = GetWithMC() ? 0xFF : L0TPData->GetTriggerFlags();
    Bool_t PhysicsTriggerOK = (L0DataType & TRIGGER_L0_PHYSICS_TYPE);
    Bool_t TriggerFlagsOK = L0TriggerFlags & MCTriggerMask;
    Bool_t TriggerOK = PhysicsTriggerOK && TriggerFlagsOK;

    // Select physics triggers and L0 trigger conditions                                
   
    // IMPORTANTE!!!!!!!!!!!!!!!! aggiungi le trigger condition del 2017!!!!!!!!!!!!!!!!              

    Int_t ID1 = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-Q2-MO1");
    Int_t ID2 = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-Q2-M1");
    Int_t ID3 = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-QX-MO1");
    Int_t ID4 = TriggerConditions::GetInstance()->GetL0TriggerID("RICH-QX-M1");
    Bool_t On1 = TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0TPData, ID1);
    Bool_t On2 = TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0TPData, ID2);
    Bool_t On3 = TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0TPData, ID3);
    Bool_t On4 = TriggerConditions::GetInstance()->L0TriggerOn(RunNumber, L0TPData, ID4);
    Bool_t TriggerStreamOK = (On1 || On2 || On3 || On4);

    if (TriggerOK) {
      if ((TriggerStreamOK && !GetWithMC()) || GetWithMC()) {

	// Select two-track events                                                             
	
	std::vector<DownstreamTrack> Tracks = *(std::vector<DownstreamTrack>*) GetOutput("DownstreamTrackBuilder.Output");
	
	if (Tracks.size() == 2) {
	  
	  // Features of the tracks                                                                 
	  
	  Int_t Charge1 = Tracks[0].GetCharge();
	  Int_t Charge2 = Tracks[1].GetCharge();
	  Double_t ChiSquare1 = Tracks[0].GetChi2();
	  Double_t ChiSquare2 = Tracks[1].GetChi2();
	  TRecoSpectrometerCandidate* SpectrometerCand1 = Tracks[0].GetSpectrometerCandidate();
	  TRecoSpectrometerCandidate* SpectrometerCand2 = Tracks[1].GetSpectrometerCandidate();
	  TVector3 Mom1 = SpectrometerCand1->GetThreeMomentumBeforeMagnet();
	  TVector3 Mom2 = SpectrometerCand2->GetThreeMomentumBeforeMagnet();
	  TVector3 TotMom = Mom1 + Mom2;
	  
	  // Compute CDA of track1,2 wrt kaon axis, between each other and of two-track total momentum wrt beam axis
	  
	  fCDAcomp->SetLine1Point1(0.0, 0.0, 102000.0);     // kaon axis                      
	  fCDAcomp->SetDir1(1.2E-3, 0., 1.);
	  fCDAcomp->SetLine2Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);     // track1                                                             
	  fCDAcomp->SetDir2(Mom1);
	  fCDAcomp->ComputeVertexCDA();
	  
	  Double_t CDA1 = fCDAcomp->GetCDA();
	  Double_t Zvertex1 = fCDAcomp->GetVertex().z();     // kaon axis-track1                      
	  
	  fCDAcomp->SetLine2Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand2->yAtBeforeMagnet(10.0), 10.0);     // track2                                                     
	  fCDAcomp->SetDir2(Mom2);
	  fCDAcomp->ComputeVertexCDA();
	  
	  Double_t CDA2 = fCDAcomp->GetCDA();
	  Double_t Zvertex2 = fCDAcomp->GetVertex().z();     // kaon axis-track2                   
	  
	  fCDAcomp->SetLine1Point1(SpectrometerCand1->xAtBeforeMagnet(10.0), SpectrometerCand1->yAtBeforeMagnet(10.0), 10.0);     // track1                                                                
	  fCDAcomp->SetDir1(Mom1);
	  fCDAcomp->SetLine2Point1(SpectrometerCand2->xAtBeforeMagnet(10.0), SpectrometerCand2->yAtBeforeMagnet(10.0), 10.0);     // track2                                                     
	  fCDAcomp->SetDir2(Mom2);
	  fCDAcomp->ComputeVertexCDA();
	  
	  Double_t CDA = fCDAcomp->GetCDA();
	  TVector3 Vertex = fCDAcomp->GetVertex();
	  Double_t Zvertex = fCDAcomp->GetVertex().z();     // track1-track2                           
	  
	  fCDAcomp->SetLine1Point1(0.0, 0.0, 102000.0);     // beam axis                           
	  fCDAcomp->SetDir1(0., 0., 1.);
	  fCDAcomp->SetLine2Point1(Vertex);     // total momentum                          
	  fCDAcomp->SetDir2(TotMom);
	  fCDAcomp->ComputeVertexCDA();
	  
	  Double_t CDAMom = fCDAcomp->GetCDA();
	  Double_t ZvertexMom = fCDAcomp->GetVertex().z();     // beam axis-total momentum 
	  
	  // Compute distance of two-track momentum wrt target (extrapolation at target)               
	  
	  fDistcomp->SetLineDir(TotMom);     // total momentum                                     
	  fDistcomp->SetLinePoint1(Vertex);
	  fDistcomp->SetPoint(0., 0., 0.);     // target                                             
	  fDistcomp->ComputeDistance();
	  
	  Double_t TargetDist = fDistcomp->GetDistance();
	  
	  // Compute distance of two-track vertex wrt beam axis                                    
	  
	  fDistcomp->SetLinePoint1(0., 0., 102000.);
	  fDistcomp->SetLineDir(0., 0., 1.);
	  fDistcomp->SetPoint(Vertex);
	  fDistcomp->ComputeDistance();
	  
	  Double_t BeamlineDist = fDistcomp->GetDistance();
	  
	  // Track selection, CUT 1: Two tracks in Spectrometer acceptance                  
    
	  for (Int_t i = 0; i < 4; i++) {
	    Double_t x1 = SpectrometerCand1->xAt(zStraw[i]);
	    Double_t y1 = SpectrometerCand1->yAt(zStraw[i]);
	    Double_t r1 = sqrt(x1*x1 + y1*y1);
	    Double_t rShifted1 = sqrt(pow(x1-xStrawChamberCentre[i],2) + y1*y1);
	    Double_t x2 = SpectrometerCand2->xAt(zStraw[i]);
	    Double_t y2 = SpectrometerCand2->yAt(zStraw[i]);
	    Double_t r2 = sqrt(x2*x2 + y2*y2);
	    Double_t rShifted2 = sqrt(pow(x2-xStrawChamberCentre[i],2) + y2*y2);
	    Bool_t inAcc = false;
	    
	    if ((rShifted1 > rMinStraw && r1 < rMaxStraw) && (rShifted2 > rMinStraw && r2 < rMaxStraw))
	      inAcc = true;
	    if (inAcc) {
	      
	      // Track selection, CUT 2: Chi2 and momentum cuts                                        
	      
	      if (ChiSquare1 <= 20. && ChiSquare2 <= 20.) {
		if (SpectrometerCand1->GetNChambers() >= 3 && SpectrometerCand2->GetNChambers() >= 3) {

		  // Track selection, CUT 3: Opposite-charged tracks                                
		  
		  if (Charge1 + Charge2 == 0) {

		    // Downstream track selection, CUT 4: Extrapolation and association to CHOD    
		    
		    Bool_t CHODAssoc = (Tracks[0].CHODAssociationExists() && Tracks[1].CHODAssociationExists());
		    Double_t x1 = SpectrometerCand1->xAtAfterMagnet(zCHOD);
		    Double_t y1 = SpectrometerCand1->yAtAfterMagnet(zCHOD);
		    Double_t r1 = sqrt(x1*x1+y1*y1);
		    Double_t x2 = SpectrometerCand2->xAtAfterMagnet(zCHOD);
		    Double_t y2 = SpectrometerCand2->yAtAfterMagnet(zCHOD);
		    Double_t r2 = sqrt(x2*x2+y2*y2);
		    Bool_t inAcc = false;
		    
		    if ((r1 > rMinCHOD && r1 < rMaxCHOD) && (r2 > rMinCHOD && r2 < rMaxCHOD))
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
			      Int_t Assoc = 0;
			      
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
				
				Double_t EoP1 = Tracks[0].GetLKrEoP();
				Double_t EoP2 = Tracks[1].GetLKrEoP();
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
					fSAVMatching->SetIRCTimeCuts(99999, 99999);
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
    
						// Computation of coupling-related quantities of good HNL in the MC event                      
						
						if (GetWithMC()) {
						  Event *evt = GetMCEvent();
						  for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
						    KinePart *p = evt->GetKinePart(i);
						    if (p->GetParentID() == -1 && p->GetPDGcode() == 999 && p->GetEndProcessName() == "good") {
						      point1.SetXYZ(p->GetProdPos().X(), p->GetProdPos().Y(), p->GetProdPos().Z());
						      point2.SetXYZ(0., 0., LInitialFV);
						      momentum1.SetXYZ(p->GetInitial4Momentum().Px(), p->GetInitial4Momentum().Py(), p->GetInitial4Momentum().Pz());
						      fMN = ComputeHNLMass(p);
						      gammaTot = GammaTot(mass1, mass2, fMN, USquared);
						      HNLTau = Tau(gammaTot);
						      LReach = ComputeL(point1, point2, momentum1);
						      BR = ComputeBR(p, fMN);
						      NReachProb = ComputeNReachProb(p, HNLTau, LReach);
						      NDecayProb = ComputeNDecayProb(p, HNLTau, LFV);
						      NPiMu = GammaPiMu(fMN)*UmuSquared/(GammaPiE(fMN)*UeSquared + GammaPiMu(fMN)*UmuSquared);
						      
						      if (p->GetParticleName() == "HNLD+e+" || p->GetParticleName() == "HNLD-e-" || p->GetParticleName() == "HNLDS+e+" || p->GetParticleName() == "HNLDS-e-") 
							LeptonUSquared = UeSquared;
						      
						      else if (p->GetParticleName() == "HNLD+mu+" || p->GetParticleName() == "HNLD-mu-" || p->GetParticleName() == "HNLDS+mu+" || p->GetParticleName() == "HNLDS-mu-") 
							LeptonUSquared = UmuSquared;

						      if (p->GetProdPos().Z() < TAXDistance)
							DProdProb = DBeProdProb;
						      else if(p->GetProdPos().Z() >= TAXDistance)
							DProdProb = DCuProdProb;
						      
						      Weight = DProdProb*DDecayProb*NReachProb*NDecayProb*NPiMu*BR*LeptonUSquared;
						      
						      if (Weight >= 0. && Weight <= 10.)
							fSumGood[fIndex] += Weight;
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
    }
  }
}

void HeavyNeutrinoPiMuSelectionScan::EndOfJobUser() {

  // Retrieve plots

  fhReach    = (TH2D*)fHisto.GetTH2("hReach");
  fhDecay    = (TH2D*)fHisto.GetTH2("hDecay");
  fhWeight   = (TH2D*)fHisto.GetTH2("hWeight");

  fhReach   ->GetXaxis()->SetTitle("Log of coupling");
  fhDecay   ->GetXaxis()->SetTitle("Log of coupling");
  fhWeight  ->GetXaxis()->SetTitle("Log of coupling");

  fhReach   ->GetYaxis()->SetTitle("Reach probability");
  fhDecay   ->GetYaxis()->SetTitle("Decay probability");
  fhWeight  ->GetYaxis()->SetTitle("Weight");

  // Acceptance computation

  for (Int_t i = 0; i < fN; i++) {

    fAcc[i] = fSumGood[i]/fSumAll[i]; 
    fProb[i] = fSumAll[i]/fNevents[i];
    fYield[i] = fAcc[i]*fProb[i];

    fgGammaTot->SetPoint(i, fCouplings[i], fGammaTot[i]);
    fgTau->SetPoint(i, fCouplings[i], fTau[i]);
    fgAcc->SetPoint(i, fCouplings[i], fAcc[i]);
    fgYield->SetPoint(i, fCouplings[i], fYield[i]);
  }

  Double_t MaxYield = *std::max_element(fYield, fYield+fN);
  Int_t MaxYieldIndex = std::max_element(fYield, fYield+fN) - fYield;

  if (MaxYield*1.E18 > 2.3) {
    ExclusionFile.open("ExclusionFile.txt", ios::app);
    ExclusionFile <<TMath::Power(10,fCouplings[MaxYieldIndex])<<"\t"<<fMN/1000.<<"\n";
    ExclusionFile.close();
  }

  fgGammaTot->GetXaxis()->SetTitle("Log of coupling");
  fgTau     ->GetXaxis()->SetTitle("Log of coupling");
  fgAcc     ->GetXaxis()->SetTitle("Log of coupling");
  fgYield   ->GetXaxis()->SetTitle("Log of coupling");
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

Double_t HeavyNeutrinoPiMuSelectionScan::ComputeHNLMass(KinePart* p) {

  Double_t MN = TMath::Sqrt(p->GetInitial4Momentum().E()*p->GetInitial4Momentum().E() - p->GetInitial4Momentum().P()*p->GetInitial4Momentum().P());

  return MN;
}

Double_t HeavyNeutrinoPiMuSelectionScan::ComputeL(TVector3 p1, TVector3 p2, TVector3 mom1) {
 
  TVector3 r(p1.x() + mom1.Px()/mom1.Pz()*(p2.z()-p1.z()), p1.y() + mom1.Py()/mom1.Pz()*(p2.z()-p1.z()), p2.z());
  Double_t x = r.x()-p1.x();
  Double_t y = r.y()-p1.y();
  Double_t z = r.z()-p1.z();
  Double_t L = TMath::Sqrt(x*x + y*y + z*z);

  return L;
}

Double_t HeavyNeutrinoPiMuSelectionScan::ComputeNDecayProb(KinePart* p, Double_t tau, Double_t l) {

  Double_t cLight = 299.9; // mm/ns
  Double_t NProb = 0.;
  NProb = 1. - TMath::Exp(-l/(p->GetInitial4Momentum().Beta()*p->GetInitial4Momentum().Gamma()*cLight*tau));

  return NProb;
}


Double_t HeavyNeutrinoPiMuSelectionScan::ComputeNReachProb(KinePart* p, Double_t tau, Double_t l) {

  Double_t cLight = 299.9; // mm/ns
  Double_t NProb = 0.;
  NProb = TMath::Exp(-l/(p->GetInitial4Momentum().Beta()*p->GetInitial4Momentum().Gamma()*cLight*tau));

  return NProb;
}

Double_t HeavyNeutrinoPiMuSelectionScan::GammaPiE(Double_t MN) {

  Double_t Mpi = 139.6;
  Double_t Me = 0.511;
  Double_t GF = 1.1E-11; // MeV^-2
  Double_t fPi = 130.41; // MeV
  Double_t gamma_pi_e = 0.;

  gamma_pi_e = GF*GF*fPi*fPi/(8*TMath::Pi()*MN*MN)*(Me*Me + MN*MN)*(MN*MN - Mpi*Mpi - Me*Me)*TMath::Power(MN*MN*MN*MN + Mpi*Mpi*Mpi*Mpi + Me*Me*Me*Me - 2.*MN*Me*MN*Me - 2.*MN*Mpi*MN*Mpi - 2.*Mpi*Me*Mpi*Me,0.5)/(2*MN);
  
  return gamma_pi_e;
}

Double_t HeavyNeutrinoPiMuSelectionScan::GammaPiMu(Double_t MN) {

  Double_t Mpi = 139.6;
  Double_t Mmu = 105.;
  Double_t GF = 1.1E-11; 
  Double_t fPi = 130.41; 
  Double_t gamma_pi_mu = 0.;

  gamma_pi_mu = GF*GF*fPi*fPi/(8*TMath::Pi()*MN*MN)*(Mmu*Mmu + MN*MN)*(MN*MN - Mpi*Mpi - Mmu*Mmu)*TMath::Power(MN*MN*MN*MN + Mpi*Mpi*Mpi*Mpi + Mmu*Mmu*Mmu*Mmu - 2.*MN*Mmu*MN*Mmu - 2.*MN*Mpi*MN*Mpi - 2.*Mpi*Mmu*Mpi*Mmu,0.5)/(2*MN);
  
  return gamma_pi_mu;
}

Double_t HeavyNeutrinoPiMuSelectionScan::GammaENuNu(Double_t MN) {

  Double_t GF = 1.1E-11;
  Double_t Me = 0.511;
  Double_t r  = 4*Me*Me/(MN*MN);
  Double_t a = TMath::Power(1-r,0.5)*(1./3.-7*r/6.-r*r/24.-r*r*r/16.);
  Double_t b = r*r*(1-r*r/16.)*TMath::ATanH(TMath::Power(1-r,0.5));
  Double_t f = a+b;
  Double_t gamma_e_nu_nu = 0.;

  gamma_e_nu_nu = GF*GF*TMath::Power(MN,5)*f/(64.*TMath::Power(TMath::Pi(),3));
  
  return gamma_e_nu_nu;
}

Double_t HeavyNeutrinoPiMuSelectionScan::GammaTot(Double_t mass1, Double_t mass2, Double_t MN, Double_t USquared) {

  Double_t UeSquared = USquared/20.8;
  Double_t UmuSquared = 16.*UeSquared;
  Double_t Mpi = 139.6;
  Double_t Me = 0.511;
  Double_t Mmu = 105.;
  Double_t gammaTot = 0.;

  if ((mass1 < 1. || mass2 < 1.) && MN < Mpi + Me)       // N->eenu                               
    gammaTot = GammaENuNu(MN)*UeSquared;
  else if ((mass1 < 1. || mass2 < 1.) && MN < Mpi + Mmu)  // N->epi                               
    gammaTot = GammaENuNu(MN)*UeSquared + GammaPiE(MN)*UeSquared;
  else if (((mass1 > 1. && mass1 < 120.) || (mass2 > 1. && mass2 < 120.)) && MN >= Mpi + Mmu)  // N->pimu                                                                      
    gammaTot = GammaENuNu(MN)*UeSquared + GammaPiE(MN)*UeSquared + GammaPiMu(MN)*UmuSquared;

  // For now

  gammaTot = GammaPiE(MN)*UeSquared + GammaPiMu(MN)*UmuSquared;

  return gammaTot;
}

Double_t HeavyNeutrinoPiMuSelectionScan::Tau(Double_t gammaTot) {

  Double_t hc = 200.E-12; // MeV*mm
  Double_t cLight = 299.9;
  Double_t tau = 0.;

  tau = hc/(gammaTot*cLight);

  return tau;
}

Double_t HeavyNeutrinoPiMuSelectionScan::NIntoPiMuBR(Double_t USquared, Double_t MN, Double_t tau) {

  Double_t UeSquared = USquared/20.8;
  Double_t UmuSquared = 16.*UeSquared;
  Double_t hc = 200.E-12;
  Double_t cLight = 299.9;
  Double_t nIntoPiMuBR = 0.;

  nIntoPiMuBR = GammaPiMu(MN)*UmuSquared*tau*cLight/hc;

  return nIntoPiMuBR;
}

Double_t HeavyNeutrinoPiMuSelectionScan::PhaseSpace(Double_t mesonMass, Double_t MN, Double_t leptonMass) {

  Double_t phaseSpace = 0.;

  phaseSpace = TMath::Power(mesonMass*mesonMass-MN*MN-leptonMass*leptonMass,2) - 4.*MN*MN*leptonMass*leptonMass;

  return phaseSpace;
}

Double_t HeavyNeutrinoPiMuSelectionScan::PhaseSpaceFactor(Double_t mesonMass, Double_t MN, Double_t leptonMass, Double_t phaseSpace) {

  Double_t factor = 0.;

  if(phaseSpace > 0.)
    factor = (mesonMass*mesonMass*(MN*MN + leptonMass*leptonMass)-TMath::Power(MN*MN-leptonMass*leptonMass,2))*TMath::Power(phaseSpace,0.5)/(leptonMass*leptonMass*TMath::Power(mesonMass*mesonMass-leptonMass*leptonMass,2));

  return factor;
}

Double_t HeavyNeutrinoPiMuSelectionScan::ComputeBR(KinePart* p, Double_t MN) {

  Double_t Me = 0.511;
  Double_t Mmu = 105.;
  Double_t DMass = 1870.;
  Double_t DSMass = 1969.;
  Double_t De2BR = 9.2E-9;
  Double_t Dm2BR = 3.74E-4;
  Double_t DSe2BR = 1.4E-7;
  Double_t DSm2BR = 5.56E-3;
  Double_t mesonMass = 0.;

  if (p->GetParticleName() == "HNLD+e+" || p->GetParticleName() == "HNLD-e-")
    mesonMass = DMass;
 
  else if (p->GetParticleName() == "HNLD+mu+" || p->GetParticleName() == "HNLD-mu-")
    mesonMass = DMass;

  else if(p->GetParticleName() == "HNLDS+e+" || p->GetParticleName() == "HNLDS-e-")
    mesonMass = DSMass;
  
  else if (p->GetParticleName() == "HNLDS+mu+" || p->GetParticleName() == "HNLDS-mu-") 
    mesonMass = DSMass;
   
  Double_t phaseSpaceE = PhaseSpace(mesonMass, MN, Me);
  Double_t phaseSpaceMu = PhaseSpace(mesonMass, MN, Mmu);
  Double_t FactorE = PhaseSpaceFactor(mesonMass, MN, Me, phaseSpaceE);
  Double_t FactorMu = PhaseSpaceFactor(mesonMass, MN, Mmu, phaseSpaceMu);
  Double_t brt = 0.;

  if (mesonMass == DMass)
    brt = De2BR*FactorE + Dm2BR*FactorMu;

  else if (mesonMass == DSMass)
    brt = DSe2BR*FactorE + DSm2BR*FactorMu;

  return brt;
}

Double_t HeavyNeutrinoPiMuSelectionScan::ComputeTotalBR(KinePart* p, Double_t USquared, Double_t MN) {

  Double_t UeSquared = USquared/20.8;
  Double_t UmuSquared = 16.*UeSquared;
  Double_t Me = 0.511;
  Double_t Mmu = 105.;
  Double_t DMass = 1870.;
  Double_t DSMass = 1969.;
  Double_t De2BR = 9.2E-9;
  Double_t Dm2BR = 3.82E-4;
  Double_t DSe2BR = 1.4E-7;
  Double_t DSm2BR = 5.9E-3;
  Double_t mesonMass = 0.;
  Double_t totalBR = 0.;

  if (p->GetParticleName() == "HNLD+e+" || p->GetParticleName() == "HNLD-e-" || p->GetParticleName() == "HNLD+mu+" || p->GetParticleName() == "HNLD-mu-")
    mesonMass = DMass;

  else if(p->GetParticleName() == "HNLDS+e+" || p->GetParticleName() == "HNLDS-e-" || p->GetParticleName() == "HNLDS+mu+" || p->GetParticleName() == "HNLDS-mu-")
    mesonMass = DSMass;

  Double_t phaseSpaceE = PhaseSpace(mesonMass, MN, Me);
  Double_t phaseSpaceMu = PhaseSpace(mesonMass, MN, Mmu);
  Double_t FactorE = PhaseSpaceFactor(mesonMass, MN, Me, phaseSpaceE);
  Double_t FactorMu = PhaseSpaceFactor(mesonMass, MN, Mmu, phaseSpaceMu);

  if (mesonMass == DMass)
    totalBR = De2BR*FactorE*UeSquared + Dm2BR*FactorMu*UmuSquared;

  else if (mesonMass == DSMass)
    totalBR = DSe2BR*FactorE*UeSquared + DSm2BR*FactorMu*UmuSquared;

  return totalBR;
}

HeavyNeutrinoPiMuSelectionScan::~HeavyNeutrinoPiMuSelectionScan() {
  
  delete fCDAcomp;
  delete fDistcomp;
  delete fLAVMatching;
  delete fSAVMatching;
}
