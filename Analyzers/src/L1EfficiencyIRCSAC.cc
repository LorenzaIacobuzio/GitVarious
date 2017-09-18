#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "L1EfficiencyIRCSAC.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "TDCEvent.hh"
#include "TIRCDigi.hh"
#include "TSACDigi.hh"
#include "Persistency.hh"
#include "L0TPData.hh"
#include "L1TPData.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

#define TRIGGER_L0_PHYSICS_TYPE 1

/// \class MyAnalyzer
/// \Brief
/// \EndBrief

/// \Detailed
/// \EndDetailed

L1EfficiencyIRCSAC::L1EfficiencyIRCSAC(Core::BaseAnalysis *ba) : Analyzer(ba, "L1EfficiencyIRCSAC") {

  RequestTree("IRC", new TDCEvent, "Digis");
  RequestTree("SAC", new TDCEvent, "Digis");

  AddParam("IRCSACDen", &fIRCSACDen, 0);
  AddParam("IRCSACNum", &fIRCSACNum, 0);
  AddParam("IRCSACRFDen", &fIRCSACRFDen, 0);
  AddParam("IRCSACRFNum", &fIRCSACRFNum, 0);
  AddParam("IRCSACRVDen", &fIRCSACRVDen, 0);
  AddParam("IRCSACRVNum", &fIRCSACRVNum, 0);

  RequestL0Data();
  RequestL1Data();
}

void L1EfficiencyIRCSAC::InitHist(){
  BookHisto("hC", new TH1D("C", "Control", 10, 0., 10.));
  BookHisto("hC1", new TH2D("C1", "Control1", 10., 0., 10., 10, 0., 10.));
}

void L1EfficiencyIRCSAC::Process(Int_t){

  RawHeader* rhe = GetRawHeader();  
  TDCEvent* IRCEvent = (TDCEvent*) GetEvent("IRC", "Digis");
  TDCEvent* SACEvent = (TDCEvent*) GetEvent("SAC", "Digis");

  // Require physics trigger and mask1 at L0 (pinunu)

  L0TPData *L0TPData = GetL0Data();
  UChar_t L0DataType = L0TPData->GetDataType();
  UInt_t L0TriggerFlags = L0TPData->GetTriggerFlags();
  Bool_t L0PhysicsTrigger = (L0DataType & TRIGGER_L0_PHYSICS_TYPE);
  Bool_t L0MaskID = (L0TriggerFlags & 2);

  if (!L0PhysicsTrigger) return;
  if ((Int_t)L0DataType != 1) return;
  if (!L0MaskID) return;

  // Require flagged events

  Int_t L1Word = (rhe->GetTriggerType() & 0x00FF00) >> 8;
  Bool_t bitG = (L1Word & 64);

  if (!bitG) return;

  // Open L1 data packet, get array of masks; for mask1, get array of enabled algorithms; require for the event to pass KTAG, LAV and STRAW algos at L1

  L1TPData *L1TPData = GetL1Data();
  if (L1TPData == NULL) return;
  std::vector<L1MaskBlock> L1Infos = L1TPData->GetL0Masks();
  Int_t NL0MasksOn = L1Infos.size();
  for (Int_t i=0; i<NL0MasksOn; i++) {
    if ((Int_t)L1Infos[i].GetL0MaskID() == 1) {
      std::vector<L1AlgoBlock> L1Algos = L1Infos[i].GetL1Algorithms();
      Int_t NAlgos = L1Algos.size();
      if (NAlgos != 4) return;
      for (Int_t j=0; j<NAlgos; j++) {
	if ((Int_t)L1Algos[j].GetL1AlgoID() == 2) {
	  UChar_t QualityFlags = L1Algos[j].GetL1QualityFlags();
	  if ((QualityFlags & 64) && !(QualityFlags & 16) && !(QualityFlags & 4) && (QualityFlags & 1)) {
	    for (Int_t k=0; k<NAlgos; k++) {
	      if ((Int_t)L1Algos[k].GetL1AlgoID() == 3) {
		UChar_t QualityFlags = L1Algos[k].GetL1QualityFlags();
		if ((QualityFlags & 64) && !(QualityFlags & 16) && !(QualityFlags & 4) && (QualityFlags & 1)) {
		  for (Int_t l=0; l<NAlgos; l++) {
		    if ((Int_t)L1Algos[l].GetL1AlgoID() == 5) {
		      UChar_t QualityFlags = L1Algos[l].GetL1QualityFlags();
		      if ((QualityFlags & 64) && !(QualityFlags & 16) && !(QualityFlags & 4) && (QualityFlags & 1)) {

			// Denominator for rejection factor

			fIRCSACRFDen++;

			// Loop over IRC and SAC digis to count the number of hits
			
			Int_t NIRCSACHits = 0;
			Double_t AlgoOnlineTimeWindow = 10.;
			Int_t NIRCDigis = IRCEvent->GetNDigis();
			TClonesArray& IRCDigis = (*(IRCEvent->GetDigis()));
			for (Int_t m=0; m<NIRCDigis; m++) {
			  TIRCDigi* IRCDigi = (TIRCDigi*)IRCDigis[m];
			  if (IRCDigi->GetDetectedEdge() == 3) {
			    if (IRCDigi->GetChannelID() < 1000) {
			      Double_t FineTime =  rhe->GetFineTime() * ClockPeriod/256.;
			      Double_t EdgeTime1 = IRCDigi->GetTime();
			      Double_t DtL0_1 = fabs(FineTime - EdgeTime1);
			      if (DtL0_1 < AlgoOnlineTimeWindow) {
				NIRCSACHits++;
			      }
			    }
			  }
			}
			
			Int_t NSACDigis = SACEvent->GetNDigis();
			TClonesArray& SACDigis = (*(SACEvent->GetDigis()));
			for (Int_t m=0; m<NSACDigis; m++) {
			  TSACDigi* SACDigi = (TSACDigi*)SACDigis[m];
			  if (SACDigi->GetDetectedEdge() == 3) {
			    if (SACDigi->GetChannelID() < 1000) {
			      Double_t FineTime =  rhe->GetFineTime() * ClockPeriod/256.;
			      Double_t EdgeTime2 = SACDigi->GetTime();
			      Double_t DtL0_2 = fabs(FineTime - EdgeTime2);
			      if (DtL0_2 < AlgoOnlineTimeWindow)
				NIRCSACHits++;
			    }
			  }
			}
			
			// Select IRCSAC algo at L1
			
			for (Int_t n=0; n<NAlgos; n++) {
			  if ((Int_t)L1Algos[n].GetL1AlgoID() == 4) {
			    UChar_t QualityFlags = L1Algos[n].GetL1QualityFlags();
			    Bool_t IsProcessed = (QualityFlags & 64);
			    Bool_t EmptyPacket = (QualityFlags & 16);
			    Bool_t BadData = (QualityFlags & 4);
			    Bool_t IsPassed = (QualityFlags & 1);
			    if (IsProcessed) {
			      if (!EmptyPacket) {
				if (!BadData) {

				  // Select events with hits in both IRC and SAC; this is the denominator of the algo efficiency
				  
				  if (NIRCSACHits > 0) {
				    fIRCSACDen++;
				    
				    // Require for the event to be vetoed by the IRCSAC algo at L1; this is the numerator of the algo efficiency
				    
				    if (!IsPassed) {
				      fIRCSACNum++;
				    }
				  }
				}
			      }
			    }
			  }
			}
			
			// Numerator for rejection factor
			
			UChar_t L1TrigWord = L1Infos[i].GetL1TriggerWord();

			if (L1TrigWord & 1)
			  fIRCSACRFNum++;
			
			for (Int_t n=0; n<NAlgos; n++) {
			  if ((Int_t)L1Algos[n].GetL1AlgoID() == 4) {
			    UChar_t QualityFlags = L1Algos[n].GetL1QualityFlags();
			    Bool_t IsProcessed = (QualityFlags & 64);
			    Bool_t EmptyPacket = (QualityFlags & 16);
			    Bool_t BadData = (QualityFlags & 4);
			    Bool_t IsPassed = (QualityFlags & 1);
			    if (IsProcessed) {
			      if (!EmptyPacket) {
				if (!BadData) {

				  // Select events with no hits in both IRC and SAC; this is the denominator of the random veto
				  
				  if (NIRCSACHits == 0) {
				    fIRCSACRVDen++;

				    // Require for the event to be vetoed by the IRCSAC algo at L1; this is the numerator of the random veto

				    if (!IsPassed) {
				      fIRCSACRVNum++;
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

void L1EfficiencyIRCSAC::EndOfJobUser() {
  SaveAllPlots();
  cout<<"Algo efficiency = "<<(Double_t)fIRCSACNum/(Double_t)fIRCSACDen<<endl;
  cout<<"Rejection factor = "<<(Double_t)fIRCSACRFNum/(Double_t)fIRCSACRFDen<<endl;
  cout<<"Random veto = "<<(Double_t)fIRCSACRVNum/(Double_t)fIRCSACRVDen<<endl;
  cout<<"Bad offline events = "<<(Double_t)fIRCSACDen<<", good offline events = "<<(Double_t)fIRCSACRVDen<<endl;
}

L1EfficiencyIRCSAC::~L1EfficiencyIRCSAC(){
}
