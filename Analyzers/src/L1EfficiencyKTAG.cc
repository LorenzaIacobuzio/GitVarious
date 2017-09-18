#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <TChain.h>
#include "L1EfficiencyKTAG.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "TDCEvent.hh"
#include "TCedarDigi.hh"
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

L1EfficiencyKTAG::L1EfficiencyKTAG(Core::BaseAnalysis *ba) : Analyzer(ba, "L1EfficiencyKTAG") {

  RequestTree("Cedar", new TDCEvent, "Digis");

  AddParam("CedarDen", &fCedarDen, 0);
  AddParam("CedarNum", &fCedarNum, 0);

  RequestL0Data();
  RequestL1Data();
}

void L1EfficiencyKTAG::InitHist(){
  BookHisto("hC", new TH1D("C", "Control", 100, 0., 100.));
  BookHisto("hC1", new TH2D("C1", "Control1", 100., 0., 100., 100, 0., 100.));
}

void L1EfficiencyKTAG::Process(Int_t){

  EventHeader* rhe = GetEventHeader();  
  TDCEvent* CedarEvent = (TDCEvent*) GetEvent("Cedar", "Digis");

  // Require physics trigger and mask1 at L0 (pinunu)

  L0TPData *L0TPData = GetL0Data();
  UChar_t L0DataType = L0TPData->GetDataType();
  UInt_t L0TriggerFlags = L0TPData->GetTriggerFlags();
  Bool_t L0PhysicsTrigger = (L0DataType & TRIGGER_L0_PHYSICS_TYPE);
  Bool_t L0MaskID = (L0TriggerFlags & 0xFF);
  Int_t L1TriggerWord = (GetEventHeader()->GetTriggerType()&0x00ff00)>>8;
  Bool_t L1AutoPass = ((L1TriggerWord & 128) == 128);

  if (!L0PhysicsTrigger) return;
  if (!L1AutoPass) return;
  if ((Int_t)L0DataType != 1) return;
  if (!L0MaskID) return;

  FillHisto("hC1", (Int_t)L0TriggerFlags, (Int_t)L0MaskID);
  /*
  Int_t Sectors[8] = {0,0,0,0,0,0,0,0};
  Double_t AlgoOnlineTimeWindow = 5.;
  Int_t NCedarDigis = CedarEvent->GetNDigis();
  TClonesArray& CedarDigis = (*(CedarEvent->GetDigis()));
  for (Int_t m=0; m<NCedarDigis; m++) {
    TCedarDigi* CedarDigi = (TCedarDigi*)CedarDigis[m];
    if (CedarDigi->GetDetectedEdge() & 1) {
      Double_t FineTime =  rhe->GetFineTime() * ClockPeriod/256.;
      Double_t EdgeTime1 = CedarDigi->GetTime();
      Double_t DtL0_1 = fabs(FineTime - EdgeTime1);
      Int_t Channel = CedarDigi->GetChannelID();
      if (DtL0_1 < AlgoOnlineTimeWindow) {
	Int_t pos = Channel/100 - 1;
	Sectors[pos]++;
      }
    }
  }
  
  Int_t counter = 0;

  for (Int_t k=0; k<8; k++) {
    if (Sectors[k] > 0) {
      counter++;
    }
  }
  */
  Int_t counter=5;
  if (counter > 4) {

    L1TPData *L1TPData = GetL1Data();
    if (L1TPData == NULL) return;
    std::vector<L1MaskBlock> L1Infos = L1TPData->GetL0Masks();
    Int_t NL0MasksOn = L1Infos.size();
    for (Int_t i=0; i<NL0MasksOn; i++) {
      if ((Int_t)L1Infos[i].GetL0MaskID() == 1) {
	std::vector<L1AlgoBlock> L1Algos = L1Infos[i].GetL1Algorithms();
	for (UInt_t j=0; j<L1Algos.size(); j++) {
	  if ((Int_t)L1Algos[j].GetL1AlgoID() == 2) {
	    Int_t QualityFlags = L1Algos[j].GetL1QualityFlags();
	    Bool_t IsProcessed = (QualityFlags & 64);
	    Bool_t EmptyPacket = (QualityFlags & 16);
	    Bool_t BadData = (QualityFlags & 4);
	    Bool_t IsPassed = (QualityFlags & 1);
	    if (IsProcessed) {
	      if (!EmptyPacket) {
		if (!BadData) {
		  fCedarDen++;
		  if (IsPassed) {
		    fCedarNum++;
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
  
void L1EfficiencyKTAG::EndOfJobUser() {
  SaveAllPlots();
  cout<<"REJECTION FACTOR = " << (Double_t)fCedarDen << "/" << (Double_t)fCedarNum << " = " << (Double_t)fCedarDen/(Double_t)fCedarNum<<endl;
}

L1EfficiencyKTAG::~L1EfficiencyKTAG(){
}
