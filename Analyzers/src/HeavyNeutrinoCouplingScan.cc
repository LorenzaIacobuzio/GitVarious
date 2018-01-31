#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPad.h>
#include "HeavyNeutrinoCouplingScan.hh"
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
#include "HNLFunctions.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

#define TRIGGER_L0_PHYSICS_TYPE 1
#define MCTriggerMask 0xFF
#define LabelSize 0.05

/// \class HeavyNeutrinoCouplingScan

HeavyNeutrinoCouplingScan::HeavyNeutrinoCouplingScan(Core::BaseAnalysis *ba) :
  Analyzer(ba, "HeavyNeutrinoCouplingScan") {
  
  RequestAllMCTrees();
  RequestAllRecoTrees();
  RequestL0Data();
  RequestL1Data();

  AddParam("USquared", &fUSquared, 1.E-6);
  AddParam("UeSquaredRatio", &fUeSquaredRatio, 1.);
  AddParam("UmuSquaredRatio", &fUmuSquaredRatio, 16.);
  AddParam("UtauSquaredRatio", &fUtauSquaredRatio, 3.8);

  // Other parameters                                                                                   

  fCDAcomp     = new TwoLinesCDA();
  fDistcomp    = new PointLineDistance();
  fLAVMatching = new LAVMatching();
  fSAVMatching = new SAVMatching();

  // Scan variables

  for (Int_t i = 0; i < fN; i++) {
    fSumGood[i]   = 0.;
    fSumAll[i]    = 0.;
    fNevents[i]   = 0;
    fCouplings[i] = 0.;
    fAcc[i]       = 0.;
    fGammaTot[i]  = 0.;
    fTau[i]       = 0.;
    fProb[i]      = 0.;
    fYield[i]     = 0.;
  }

  // Histos

  fhReach     = nullptr;
  fhDecay     = nullptr;
  fhWeight    = nullptr;
  fgAcc       = nullptr;
  fgYield     = nullptr;
  fgGammaTot  = nullptr;
  fgTau       = nullptr;
}

void HeavyNeutrinoCouplingScan::InitHist() {

  BookHisto("hReach",    new TH2D("Reach", "Probability of N reaching the FV vs coupling",    fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("hDecay",    new TH2D("Decay", "Probability of N decaying in the FV vs coupling", fN, fCouplingStart, fCouplingStop, 1000, -0.1, 1.1));
  BookHisto("hWeight",   new TH2D("Weight", "N weight vs coupling",                           fN, fCouplingStart, fCouplingStop, 1000, 1.E-26, 1.E-21));

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

void HeavyNeutrinoCouplingScan::Process(Int_t) {

  TRecoLAVEvent*          LAVEvent  = (TRecoLAVEvent*)          GetEvent("LAV");
  TRecoIRCEvent*          IRCEvent  = (TRecoIRCEvent*)          GetEvent("IRC");
  TRecoSACEvent*          SACEvent  = (TRecoSACEvent*)          GetEvent("SAC");

  Double_t MN             = 0.;
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

  // Scan on the coupling                                                               

  for(Double_t fCoupling = fCouplingStart; fCoupling < fCouplingStop-fCouplingStep; fCoupling += fCouplingStep) {
    
    fUSquared = TMath::Power(10, fCoupling);
    fUeSquared   = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUeSquaredRatio;
    fUmuSquared  = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUmuSquaredRatio;
    fUtauSquared = fUSquared/(fUeSquaredRatio + fUmuSquaredRatio + fUtauSquaredRatio)*fUtauSquaredRatio;
    fIndex = round((std::abs(fCouplingStart-fCoupling))/fCouplingStep);
    fCouplings[fIndex] = fCoupling;
    
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
      
      // Computation of coupling-related quantities of all HNLs (good and bad)
      
      for (Int_t i = 0; i < evt->GetNKineParts(); i++) {
	KinePart *p = evt->GetKinePart(i);      
	if (p->GetParentID() == -1 && p->GetPDGcode() == 999) {
	  fNevents[fIndex]++;
	  point1.SetXYZ(p->GetProdPos().X(), p->GetProdPos().Y(), p->GetProdPos().Z());
	  point2.SetXYZ(0., 0., fLInitialFV);
	  momentum1.SetXYZ(p->GetInitial4Momentum().Px(), p->GetInitial4Momentum().Py(), p->GetInitial4Momentum().Pz());
	  MN = ComputeHNLMass(p);
	  gammaTot = GammaTot(MN);
	  HNLTau = tauN(MN);
	  fGammaTot[fIndex] = gammaTot;
	  fTau[fIndex] = HNLTau;
	  LReach = ComputeL(point1, point2, momentum1);
	  ProdFactor = ComputeProd(p, MN);
	  DecayFactor = ComputeDecay(MN);
	  NReachProb = ComputeNReachProb(p, HNLTau, LReach);
	  NDecayProb = ComputeNDecayProb(p, HNLTau, fLFV);
	  
	  if (p->GetParticleName().Contains("e") && !p->GetParticleName().Contains("nu_tau"))
	    LeptonUSquared = fUeSquared;
	  else if (p->GetParticleName().Contains("mu") && !p->GetParticleName().Contains("nu_tau"))
	    LeptonUSquared = fUmuSquared;
	  else if (p->GetParticleName() == "DS->Ntau" || p->GetParticleName().Contains("nu_tau") || (p->GetParticleName().Contains("tau") && (p->GetParticleName().Contains("rho") || p->GetParticleName().Contains("pi") || p->GetParticleName().Contains("K"))))
	    LeptonUSquared = fUtauSquared;
	  
	  if (p->GetProdPos().Z() < fTAXDistance)
	    DProdProb = fDBeProdProb;
	  else if (p->GetProdPos().Z() >= fTAXDistance)
	    DProdProb = fDCuProdProb;
	  
	  // Weight to be associated to each HNL
	  
	  Weight = DProdProb*fDDecayProb*NReachProb*NDecayProb*DecayFactor*ProdFactor*LeptonUSquared;
	  fSumAll[fIndex] += Weight;

	  FillHisto("hReach",  fCoupling, NReachProb);
	  FillHisto("hDecay",  fCoupling, NDecayProb);
	  FillHisto("hWeight", fCoupling, Weight);
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
						  MN = ComputeHNLMass(p);
						  gammaTot = GammaTot(MN);
						  HNLTau = tauN(MN);
						  LReach = ComputeL(point1, point2, momentum1);
						  ProdFactor = ComputeProd(p, MN);
						  DecayFactor = ComputeDecay(MN);
						  NReachProb = ComputeNReachProb(p, HNLTau, LReach);
						  NDecayProb = ComputeNDecayProb(p, HNLTau, fLFV);
						  
						  if (p->GetParticleName().Contains("e") && !p->GetParticleName().Contains("nu_tau"))
						    LeptonUSquared = fUeSquared;
						  else if (p->GetParticleName().Contains("mu") && !p->GetParticleName().Contains("nu_tau"))
						    LeptonUSquared = fUmuSquared;
						  else if (p->GetParticleName() == "DS->Ntau" || p->GetParticleName().Contains("nu_tau") || (p->GetParticleName().Contains("tau") && (p->GetParticleName().Contains("rho") || p->GetParticleName().Contains("pi") || p->GetParticleName().Contains("K"))))
						    LeptonUSquared = fUtauSquared;
						  
						  if (p->GetProdPos().Z() < fTAXDistance)
						    DProdProb = fDBeProdProb;
						  else if(p->GetProdPos().Z() >= fTAXDistance)
						    DProdProb = fDCuProdProb;
						  
						  Weight = DProdProb*fDDecayProb*NReachProb*NDecayProb*DecayFactor*ProdFactor*LeptonUSquared;
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

void HeavyNeutrinoCouplingScan::EndOfJobUser() {
  
  // Retrieve histos

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
    fgTau     ->SetPoint(i, fCouplings[i], fTau[i]);
    fgAcc     ->SetPoint(i, fCouplings[i], fAcc[i]);
    fgYield   ->SetPoint(i, fCouplings[i], fYield[i]);
  }

  // Set titles, etc.

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

HeavyNeutrinoCouplingScan::~HeavyNeutrinoCouplingScan() {
  
  delete fCDAcomp;
  delete fDistcomp;
  delete fLAVMatching;
  delete fSAVMatching;

  fhReach    = nullptr;
  fhDecay    = nullptr;
  fhWeight   = nullptr;
  fgAcc      = nullptr;
  fgYield    = nullptr;
  fgGammaTot = nullptr;
  fgTau      = nullptr;
}
