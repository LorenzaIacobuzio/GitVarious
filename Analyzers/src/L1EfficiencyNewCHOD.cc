#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "L1EfficiencyNewCHOD.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "TDCEvent.hh"
#include "TNewCHODDigi.hh"
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

L1EfficiencyNewCHOD::L1EfficiencyNewCHOD(Core::BaseAnalysis *ba) : Analyzer(ba, "L1EfficiencyNewCHOD") {

  RequestTree("NewCHOD", new TDCEvent, "Digis");

  AddParam("NewCHODDen", &fNewCHODDen, 0);
  AddParam("NewCHODNum", &fNewCHODNum, 0);
  AddParam("NewCHODRFDen", &fNewCHODRFDen, 0);
  AddParam("NewCHODRFNum", &fNewCHODRFNum, 0);
  AddParam("NewCHODRVDen", &fNewCHODRVDen, 0);
  AddParam("NewCHODRVNum", &fNewCHODRVNum, 0);

  RequestL0Data();
  RequestL1Data();
}

void L1EfficiencyNewCHOD::InitHist(){
  BookHisto("hC", new TH1D("C", "Control", 10, 0., 10.));
  BookHisto("hC1", new TH2D("C1", "Control1", 10., 0., 10., 10, 0., 10.));
}

void L1EfficiencyNewCHOD::Process(Int_t){

  RawHeader* rhe = GetRawHeader();
  TDCEvent* NewCHODEvent = (TDCEvent*)GetEvent("NewCHOD", "Digis");  

  // Require physics trigger and mask0 at L0 (newchod) 

  L0TPData *L0TPData = GetL0Data();
  UChar_t L0DataType = L0TPData->GetDataType();
  UInt_t L0TriggerFlags = L0TPData->GetTriggerFlags();
  Bool_t L0PhysicsTrigger = (L0DataType & TRIGGER_L0_PHYSICS_TYPE);
  Bool_t L0MaskID = (L0TriggerFlags & 1);

  if (!L0PhysicsTrigger) return;
  if ((Int_t)L0DataType != 1) return;
  if (!L0MaskID) return;

  // Require flagged events 
  
  Int_t L1Word = (rhe->GetTriggerType() & 0x00FF00) >> 8;
  Bool_t bitG = (L1Word & 64);

  if (!bitG) return;


  // Open L1 data packet, get array of masks; for mask0, get array of enabled algorithms

  L1TPData *L1TPData = GetL1Data();
  if (L1TPData == NULL) return;
  std::vector<L1MaskBlock> L1Infos = L1TPData->GetL0Masks();
  Int_t NL0MasksOn = L1Infos.size();
  for (Int_t i=0; i<NL0MasksOn; i++) {
    if ((Int_t)L1Infos[i].GetL0MaskID() == 0) {
      std::vector<L1AlgoBlock> L1Algos = L1Infos[i].GetL1Algorithms();
      Int_t NAlgos = L1Algos.size();
      if (NAlgos != 1) return;
	
      // Denominator for rejection factor
      
      fNewCHODRFDen++;
      
      // Loop over NewCHOD digis to count the number of hits                       
      
      Int_t NNewCHODHits = 0;
      Double_t AlgoOnlineTimeWindow = 10.;
      Int_t NNewCHODDigis = NewCHODEvent->GetNDigis();
      TClonesArray& NewCHODDigis = (*(NewCHODEvent->GetDigis()));
      for (Int_t i=0; i<NNewCHODDigis; i++) {
	TNewCHODDigi* NewCHODDigi1 = (TNewCHODDigi*)NewCHODDigis[i];
	if (NewCHODDigi1->GetDetectedEdge() & 1) {
	  Int_t pmtID1 = NewCHODDigi1->GetChannelID();
	  if ((pmtID1/10) % 10 >= 0 && (pmtID1/10) % 10 <= 3) {
	    for (Int_t j=0; j<NNewCHODDigis; j++) {
	      TNewCHODDigi* NewCHODDigi2 = (TNewCHODDigi*)NewCHODDigis[j];
	      if ((NewCHODDigi2->GetDetectedEdge() & 1) && i != j) {
		Int_t pmtID2 = NewCHODDigi2->GetChannelID();
		if (fabs(pmtID1-pmtID2) == 50 && (pmtID1/100) % 10 == (pmtID2/100) % 10) {
		  Double_t FineTime =  rhe->GetFineTime() * ClockPeriod/256.;
		  Double_t EdgeTime1 = NewCHODDigi1->GetTime();
		  Double_t EdgeTime2 = NewCHODDigi2->GetTime();
		  Double_t DtL0_1 = fabs(EdgeTime1 - FineTime);
		  Double_t DtL0_2 = fabs(EdgeTime2 - FineTime);
		  if (DtL0_1 < AlgoOnlineTimeWindow && DtL0_2 < AlgoOnlineTimeWindow && fabs(EdgeTime1 - EdgeTime2) < 5.) {
		    NNewCHODHits++;
		  }		
		}
	      }
	    }
	  }
	}
      }
      
      for (Int_t j=0; j<NAlgos; j++) {
	if ((Int_t)L1Algos[j].GetL1AlgoID() == 7) {
	  UChar_t QualityFlags = L1Algos[j].GetL1QualityFlags();
	  Bool_t IsProcessed = (QualityFlags & 64);
	  Bool_t EmptyPacket = (QualityFlags & 16);
	  Bool_t BadData = (QualityFlags & 4);
	  Bool_t IsPassed = (QualityFlags & 1);
	  if (IsProcessed) {
	    if (!EmptyPacket) {
	      if (!BadData) {
		
		// Select events with 0 < nhits < 4; this is the denominator of the algo efficiency
		
		if (NNewCHODHits > 0 && NNewCHODHits < 4) {
		  fNewCHODDen++;
		  
		  // Require for the event to pass the NewCHOD algo at L1; this is the numerator of the algo efficiency
		  
		  if (IsPassed) {
		    fNewCHODNum++;
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
	fNewCHODRFNum++;
      
      for (Int_t j=0; j<NAlgos; j++) {
	if ((Int_t)L1Algos[j].GetL1AlgoID() == 7) {
	  UChar_t QualityFlags = L1Algos[j].GetL1QualityFlags();
	  Bool_t IsProcessed = (QualityFlags & 64);
	  Bool_t EmptyPacket = (QualityFlags & 16);
	  Bool_t BadData = (QualityFlags & 4);
	  Bool_t IsPassed = (QualityFlags & 1);
	  if (IsProcessed) {
	    if (!EmptyPacket) {
	      if (!BadData) {
		
		// Select events with more than 4 hits in NewCHOD; this is the denominator of the mistagging  
		
		if (NNewCHODHits >= 4) {
		  fNewCHODRVDen++;
		  
		  // Require for the event to pass the NewCHOD algo at L1; this is the numerator of the mistagging       

		  if (IsPassed) {
		    fNewCHODRVNum++;
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

void L1EfficiencyNewCHOD::EndOfJobUser() {
    SaveAllPlots();
    cout<<"Algo efficiency = "<<(Double_t)fNewCHODNum/(Double_t)fNewCHODDen<<endl;
    cout<<"Rejection factor = "<<(Double_t)fNewCHODRFNum/(Double_t)fNewCHODRFDen<<endl;
    cout<<"Mistagging = "<<(Double_t)fNewCHODRVNum/(Double_t)fNewCHODRVDen<<endl;
    cout<<"Good offline events = "<<(Double_t)fNewCHODDen<<", bad offline events = "<<(Double_t)fNewCHODRVDen<<endl;
}

L1EfficiencyNewCHOD::~L1EfficiencyNewCHOD(){
}
