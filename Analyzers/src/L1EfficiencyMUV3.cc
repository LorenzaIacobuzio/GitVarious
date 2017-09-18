#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "L1EfficiencyMUV3.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "TDCEvent.hh"
#include "TMUV3Digi.hh"
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

L1EfficiencyMUV3::L1EfficiencyMUV3(Core::BaseAnalysis *ba) : Analyzer(ba, "L1EfficiencyMUV3") {

  RequestTree("MUV3", new TDCEvent, "Digis");

  AddParam("MUV3Den", &fMUV3Den, 0);
  AddParam("MUV3Num", &fMUV3Num, 0);
  AddParam("MUV3RFDen", &fMUV3RFDen, 0);
  AddParam("MUV3RFNum", &fMUV3RFNum, 0);
  AddParam("MUV3RVDen", &fMUV3RVDen, 0);
  AddParam("MUV3RVNum", &fMUV3RVNum, 0);

  RequestL0Data();
  RequestL1Data();
}

void L1EfficiencyMUV3::InitHist(){
  BookHisto("hC", new TH1D("C", "Control", 10, 0., 10.));
  BookHisto("hC1", new TH2D("C1", "Control1", 10., 0., 10., 10, 0., 10.));
}

void L1EfficiencyMUV3::Process(Int_t){

  RawHeader* rhe = GetRawHeader();  
  TDCEvent* MUV3Event = (TDCEvent*) GetEvent("MUV3", "Digis");

  //********* Multiplicity algo ***********
  /*
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

                        fMUV3RFDen++;

                        // Loop over MUV3 digis to count the number of hits
   
			Int_t NMUV3Hits = 0;
			Double_t AlgoOnlineTimeWindow = 6.;
			Int_t NMUV3Digis = MUV3Event->GetNDigis();
			TClonesArray& MUV3Digis = (*(MUV3Event->GetDigis()));
			for (Int_t m=0; m<NMUV3Digis; m++) {
			  TMUV3Digi* MUV3Digi = (TMUV3Digi*)MUV3Digis[m];
			  if (MUV3Digi->GetDetectedEdge() & 1) {
			    Double_t FineTime =  rhe->GetFineTime() * ClockPeriod/256.;
			    Double_t EdgeTime = MUV3Digi->GetTime();
			    Double_t DtL0 = fabs(FineTime - EdgeTime);
			    if (DtL0 < AlgoOnlineTimeWindow)
			      NMUV3Hits++;
			  }
			}
			
			for (Int_t n=0; n<NAlgos; n++) {
			  if ((Int_t)L1Algos[n].GetL1AlgoID() == 6) {
			    UChar_t QualityFlags = L1Algos[n].GetL1QualityFlags();
			    Bool_t IsProcessed = (QualityFlags & 64);
			    Bool_t EmptyPacket = (QualityFlags & 16);
			    Bool_t BadData = (QualityFlags & 4);
			    Bool_t IsPassed = (QualityFlags & 1);
			    if (IsProcessed) {
			      if (!EmptyPacket) {
				if (!BadData) {

				  // Select events with hits in MUV3; this is the denominator of the algo efficiency
				  
				  if (NMUV3Hits > 0) {
				    fMUV3Den++;			  

				    // Require for the event to be vetoed by the MUV3 algo at L1; this is the numerator of the algo efficiency 
			
				    if (!IsPassed) {
				      fMUV3Num++;
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
                          fMUV3RFNum++;
			
			for (Int_t n=0; n<NAlgos; n++) {
			  if ((Int_t)L1Algos[n].GetL1AlgoID() == 6) {
			    UChar_t QualityFlags = L1Algos[n].GetL1QualityFlags();
			    Bool_t IsProcessed = (QualityFlags & 64);
			    Bool_t EmptyPacket = (QualityFlags & 16);
			    Bool_t BadData = (QualityFlags & 4);
			    Bool_t IsPassed = (QualityFlags & 1);
			    if (IsProcessed) {
			      if (!EmptyPacket) {
				if (!BadData) {

				  // Select events with no hits in MUV3; this is the denominator of the random veto
				  
				  if (NMUV3Hits == 0) {
				    fMUV3RVDen++;				    

				    // Require for the event to be vetoed by the MUV3 algo at L1; this is the numerator of the random veto

				    if (!IsPassed) {
				      fMUV3RVNum++;
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
  */
  // ************* Left-right algo **************
  /*  
  // Require physics trigger and mask6 at L0 (exotics)                                                                     

  L0TPData *L0TPData = GetL0Data();
  UChar_t L0DataType = L0TPData->GetDataType();
  UInt_t L0TriggerFlags = L0TPData->GetTriggerFlags();
  Bool_t L0PhysicsTrigger = (L0DataType & TRIGGER_L0_PHYSICS_TYPE);
  Bool_t L0MaskID = (L0TriggerFlags & 64);

  if (!L0PhysicsTrigger) return;
  if ((Int_t)L0DataType != 1) return;
  if (!L0MaskID) return;

  // Require flagged events                                                                                               

  Int_t L1Word = (rhe->GetTriggerType() & 0x00FF00) >> 8;
  Bool_t bitG = (L1Word & 64);

  if (!bitG) return;

  // Open L1 data packet, get array of masks; for mask6, get array of enabled algorithms; require for the event to pass KTAG, LAV and STRAW algos at L1                      

  L1TPData *L1TPData = GetL1Data();
  if (L1TPData == NULL) return;
  std::vector<L1MaskBlock> L1Infos = L1TPData->GetL0Masks();
  Int_t NL0MasksOn = L1Infos.size();
  for (Int_t i=0; i<NL0MasksOn; i++) {
    if ((Int_t)L1Infos[i].GetL0MaskID() == 6) {
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

                        fMUV3RFDen++;

                        // Loop over MUV3 digis to check whether there is a left-right couple    

			Bool_t Couple = false;
			Double_t AlgoOnlineTimeWindow = 6.;
			Int_t NMUV3Digis = MUV3Event->GetNDigis();
			TClonesArray& MUV3Digis = (*(MUV3Event->GetDigis()));
			for (Int_t m=0; m<NMUV3Digis; m++) {
			  TMUV3Digi* MUV3Digi1 = (TMUV3Digi*)MUV3Digis[m];
			  if (MUV3Digi1->GetDetectedEdge() & 1) {
			    Int_t pmtID1 = MUV3Digi1->GetChannelID();
			    if (pmtID1 > 151)
			      pmtID1 = pmtID1 - 200;
			    for (Int_t n=0; n<NMUV3Digis; n++) {
			      if (!Couple) {
				TMUV3Digi* MUV3Digi2 = (TMUV3Digi*)MUV3Digis[n];
				if ((MUV3Digi2->GetDetectedEdge() & 1) && m != n) {
				  Int_t pmtID2 = MUV3Digi2->GetChannelID();
				  if (pmtID2 > 151)
				    pmtID2 = pmtID2 - 200;
				  if (pmtID1 < 144 && pmtID2 < 144) {
				    if (((pmtID1 % 12) <= 5 && (pmtID2 % 12) >= 6) || ((pmtID2 % 12) <= 5 && (pmtID1 % 12) >= 6)) {
				      Double_t FineTime =  rhe->GetFineTime() * ClockPeriod/256.;
				      Double_t EdgeTime1 = MUV3Digi1->GetTime();
				      Double_t DtL0_1 = fabs(FineTime - EdgeTime1);
				      Double_t EdgeTime2 = MUV3Digi2->GetTime();
				      Double_t DtL0_2 = fabs(FineTime - EdgeTime2);
				      if (fabs(EdgeTime1 - EdgeTime2) < 10.) {
					if (DtL0_1 < AlgoOnlineTimeWindow && DtL0_2 < AlgoOnlineTimeWindow) {
					  Couple = true;
					}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		     
			for (Int_t p=0; p<NAlgos; p++) {
			  if ((Int_t)L1Algos[p].GetL1AlgoID() == 6) {
			    UChar_t QualityFlags = L1Algos[p].GetL1QualityFlags();
			    Bool_t IsProcessed = (QualityFlags & 64);
			    Bool_t EmptyPacket = (QualityFlags & 16);
			    Bool_t BadData = (QualityFlags & 4);
			    Bool_t IsPassed = (QualityFlags & 1);
			    if (IsProcessed) {
			      if (!EmptyPacket) {
				if (!BadData) {
				  
				  // Select events with a left-right couple in MUV3; this is the denominator of the algo efficiency    
				  
				  if (Couple) {
				    fMUV3Den++;

				    // Require for the event to pass the MUV3 algo at L1; this is the numerator of the algo efficiency                                                                                                             
				    if (IsPassed) {
				      fMUV3Num++;
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
                          fMUV3RFNum++;

			for (Int_t p=0; p<NAlgos; p++) {
			  if ((Int_t)L1Algos[p].GetL1AlgoID() == 6) {
			    UChar_t QualityFlags = L1Algos[p].GetL1QualityFlags();
			    Bool_t IsProcessed = (QualityFlags & 64);
			    Bool_t EmptyPacket = (QualityFlags & 16);
			    Bool_t BadData = (QualityFlags & 4);
			    Bool_t IsPassed = (QualityFlags & 1);
			    if (IsProcessed) {
			      if (!EmptyPacket) {
				if (!BadData) {

				  // Select events with no couples in MUV3; this is the denominator of the mistagging                                                                                                                   
				  if (!Couple) {
				    fMUV3RVDen++;

				    // Require for the event to pass the MUV3 algo at L1; this is the numerator of the mistagging                                                                                                      
				    if (IsPassed) {
				      fMUV3RVNum++;
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
  */
  // ************ Neighbours algo *************
  
  // Require physics trigger and mask6 at L0 (exotics)                                                                    

  L0TPData *L0TPData = GetL0Data();
  UChar_t L0DataType = L0TPData->GetDataType();
  UInt_t L0TriggerFlags = L0TPData->GetTriggerFlags();
  Bool_t L0PhysicsTrigger = (L0DataType & TRIGGER_L0_PHYSICS_TYPE);
  Bool_t L0MaskID = (L0TriggerFlags & 64);

  if (!L0PhysicsTrigger) return;
  if ((Int_t)L0DataType != 1) return;
  if (!L0MaskID) return;

  // Require flagged events                                                                                               

  Int_t L1Word = (rhe->GetTriggerType() & 0x00FF00) >> 8;
  Bool_t bitG = (L1Word & 64);

  if (!bitG) return;

  // Open L1 data packet, get array of masks; for mask6, get array of enabled algorithms; require for the event to pass KTAG, LAV and STRAW algos at L1                      

  L1TPData *L1TPData = GetL1Data();
  if (L1TPData == NULL) return;
  std::vector<L1MaskBlock> L1Infos = L1TPData->GetL0Masks();
  Int_t NL0MasksOn = L1Infos.size();
  for (Int_t i=0; i<NL0MasksOn; i++) {
    if ((Int_t)L1Infos[i].GetL0MaskID() == 6) {
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
                        fMUV3RFDen++;
					      
                        // Loop over MUV3 digis to count the number of hits    

			Int_t NMUV3Hits = 0;
			Int_t Flag = false;
			Double_t AlgoOnlineTimeWindow = 6.;
			Int_t pmtIDarr[2] = {0, 0};
			Bool_t TileID[152];
			for (Int_t q=0; q<152; q++)
			  TileID[q] = false;
			Int_t NMUV3Digis = MUV3Event->GetNDigis();
			TClonesArray& MUV3Digis = (*(MUV3Event->GetDigis()));
			for (Int_t m=0; m<NMUV3Digis; m++) {
			  TMUV3Digi* MUV3Digi = (TMUV3Digi*)MUV3Digis[m];
			  if (MUV3Digi->GetDetectedEdge() & 1) {
			    Int_t pmtID = MUV3Digi->GetChannelID();
			    if (pmtID > 151)
			      pmtID = pmtID - 200;
			    if (pmtID < 144) {
			      Double_t FineTime =  rhe->GetFineTime() * ClockPeriod/256.;
			      Double_t EdgeTime = MUV3Digi->GetTime();
			      Double_t DtL0 = fabs(FineTime - EdgeTime);
			      if (DtL0 < AlgoOnlineTimeWindow)
				TileID[pmtID] = true;
			    }
			  }
			}

			for (Int_t r=0; r<152; r++) {
			  if (TileID[r]) {
			    NMUV3Hits++;
			    if (NMUV3Hits > 2)
			      Flag = true;
			    else if(NMUV3Hits <= 2)
			      pmtIDarr[NMUV3Hits-1] = r;
			  }
			}

			if (NMUV3Hits <= 1)
			  Flag = false;
			else if (NMUV3Hits == 2) {
			  if (fabs(pmtIDarr[0] - pmtIDarr[1]) == 1 && ((pmtIDarr[0] + pmtIDarr[1]) % 24) != 23)
			    Flag = false;
			  else if (fabs(pmtIDarr[0] - pmtIDarr[1]) == 12)
			    Flag = false;
			  else 
			    Flag = true;
			}
			else
			  Flag = true;

			for (Int_t n=0; n<NAlgos; n++) {
			  if ((Int_t)L1Algos[n].GetL1AlgoID() == 6) {
			    UChar_t QualityFlags = L1Algos[n].GetL1QualityFlags();
			    Bool_t IsProcessed = (QualityFlags & 64);
			    Bool_t EmptyPacket = (QualityFlags & 16);
			    Bool_t BadData = (QualityFlags & 4);
			    Bool_t IsPassed = (QualityFlags & 1);
			    if (IsProcessed) {
			      if (!EmptyPacket) {
				if (!BadData) {

				  // Select events with more than 2 hits in MUV3, or with 2 hits but not with neighbour firing tiles; this is the denominator of the algo efficiency    

				  if (Flag) {
				    fMUV3Den++;

				    // Require for the event to pass the MUV3 algo at L1; this is the numerator of the algo efficiency                                                                                                             
				    if (IsPassed) {
				      fMUV3Num++;
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
                          fMUV3RFNum++;

			for (Int_t n=0; n<NAlgos; n++) {
			  if ((Int_t)L1Algos[n].GetL1AlgoID() == 6) {
			    UChar_t QualityFlags = L1Algos[n].GetL1QualityFlags();
			    Bool_t IsProcessed = (QualityFlags & 64);
			    Bool_t EmptyPacket = (QualityFlags & 16);
			    Bool_t BadData = (QualityFlags & 4);
			    Bool_t IsPassed = (QualityFlags & 1);
			    if (IsProcessed) {
			      if (!EmptyPacket) {
				if (!BadData) {

				  // Select events with less than 2 hits in MUV3 or with 2 hits but with neighbour firing tiles; this is the denominator of the mistagging                                                                                  				  
				  if (!Flag) {                                                                                           
				    fMUV3RVDen++; 

				    // Require for the event to pass the MUV3 algo at L1; this is the numerator of the mistagging
                                                                     
				    if (IsPassed) {
				      fMUV3RVNum++;
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

void L1EfficiencyMUV3::EndOfJobUser() {
  SaveAllPlots();
  cout<<"Algo efficiency = "<<(Double_t)fMUV3Num/(Double_t)fMUV3Den<<endl;
  cout<<"Rejection factor = "<<(Double_t)fMUV3RFNum/(Double_t)fMUV3RFDen<<endl;
  cout<<"Random veto/mistagging = "<<(Double_t)fMUV3RVNum/(Double_t)fMUV3RVDen<<endl;
  cout<<"Bad offline events = "<<(Double_t)fMUV3Den<<", good offline events = "<<(Double_t)fMUV3RVDen<<endl;
}

L1EfficiencyMUV3::~L1EfficiencyMUV3(){
}
