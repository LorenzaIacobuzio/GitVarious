#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "L1DataPacket.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "L0TPData.hh"
#include "L1TPData.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

#define TRIGGER_L0_PHYSICS_TYPE 1
#define TRIGGER_L0_MONITORING_TYPE 20
#define TRIGGER_L0_CALIBRATION_TYPE 4

/// \class L1DataPacket
/// \Brief

/// Tool to retrieve information in the L1 data packet

/// \EndBrief
///
/// \Detailed

/// A list of all the methods available can be found on Doxygen: https://na62-sw.web.cern.ch/sites/na62-sw.web.cern.ch/files/doxygen/da/db2/classL1TPData.html
/// A detailed description of the packet contents can be found in the TDAQ data-format living note (section "Event format, data to storage"): https://twiki.cern.ch/twiki/bin/view/NA62/TdaqSystem
/// A detailed description of L1 algorithms, XML config files, structures, word contents, etc. can be found in the ShifterHandbook2016 ("Subsystem description - L1/L2"): https://twiki.cern.ch/twiki/bin/view/NA62/ShifterHandbook2016DescriptionL12

/// \EndDetailed

L1DataPacket::L1DataPacket(Core::BaseAnalysis *ba) : Analyzer(ba, "L1DataPacket")
{
  RequestL0Data();
  RequestL1Data();
}

void L1DataPacket::Process(Int_t){

  // Example code to retrieve L0 mask1 data

  EventHeader* rhe = GetEventHeader();

  L0TPData *L0TPData = GetL0Data();
  UChar_t L0DataType = L0TPData->GetDataType();
  UInt_t L0TriggerFlags = L0TPData->GetTriggerFlags();
  Bool_t L0PhysicsTrigger = (L0DataType & TRIGGER_L0_PHYSICS_TYPE);   // requiring physics triggers
  Bool_t L0MaskID = (L0TriggerFlags & 2);   // requiring mask1 at L0; L0TriggerFlags = 2^(L0Mask) (eg. for mask 4, triggerflags = 16)

  if (!L0PhysicsTrigger) return;
  if ((Int_t)L0DataType != 1) return;
  if (!L0MaskID) return;

  Int_t L1Word = (rhe->GetTriggerType() & 0x00FF00) >> 8;   // retrieving L1TriggerWord; this is a global, event-related word which contains information about flagging, autopass, bypass, verdict, etc.
  Bool_t bitG = (L1Word & 64);   // asking for flagged events

  if (!bitG) return;

  L1TPData *L1TPData = GetL1Data();   // retrieving L1 data packet
  if (L1TPData == NULL) return;   // every data generally related to the L1 packet of the event can be found at this stage (eg. timestamp, downscaling factor, reference detector, etc.); for documentation, see https://na62-sw.web.cern.ch/sites/na62-sw.web.cern.ch/files/doxygen/da/db2/classL1TPData.html
  std::vector<L1MaskBlock> L1Infos = L1TPData->GetL0Masks();   // retrieving vector of L1 masks; data related to the mask packet can be found here (eg. mask trigger word, number of enabled algos, etc.); for documentation, see https://na62-sw.web.cern.ch/sites/na62-sw.web.cern.ch/files/doxygen/dc/da6/classL1MaskBlock.html
  Int_t NL0MasksOn = L1Infos.size();
  for (Int_t iMask=0; iMask<NL0MasksOn; iMask++) {
    if ((Int_t)L1Infos[iMask].GetL0MaskID() == 1) {   // looping over masks and choosing the mask ID at L1
      std::vector<L1AlgoBlock> L1Algos = L1Infos[iMask].GetL1Algorithms();   // retrieving vector of L1 algorithms; a vector is associated to each mask;  each algo has its own ID (eg. CHOD = 0, RICH = 1, etc.); every data related to the algo packet can be found here (eg. quality flags, algo flags, time window, etc.); for documentation, see https://na62-sw.web.cern.ch/sites/na62-sw.web.cern.ch/files/doxygen/d2/d22/classL1AlgoBlock.html
      Int_t NAlgos = L1Algos.size();
      for (Int_t jAlgo=0; jAlgo<NAlgos; jAlgo++) {  
	if ((Int_t)L1Algos[jAlgo].GetL1AlgoID() == 2) {   // looping over algos and choosing the algo ID
          UChar_t QualityFlags = L1Algos[jAlgo].GetL1QualityFlags();   // retrieving the quality flags for each algo; each quality-flag word contains information about whether:
	  Bool_t IsProcessed = (QualityFlags & 64);   // the algo was processed for the current event
	  Bool_t EmptyPacket = (QualityFlags & 16);   // something wrong at L1 generated an empty packet for the event
	  Bool_t BadData = (QualityFlags & 4);   // some corruption occurred and the packet is bad
	  Bool_t IsPassed = (QualityFlags & 1);   // the current event passed the selection of the algo (eg. if the algo looks for >4 sectors in the KTAG and the event has 6 sectors, bit 0 of the quality-flag word is raised)
	  UChar_t AlgoFlags = L1Algos[jAlgo].GetL1AlgoFlags();   // retrieving the algo flags for each algo; it contains information about whether:
	  Bool_t AlgoLogic = (AlgoFlags & 4);   // the logic was set to positive (1) or negative (0)
	  Bool_t IsEnabled = (AlgoFlags & 16);   // the algo was enabled

	  std::vector<UInt_t> L1AlgoWords = L1Algos[jAlgo].GetL1DataWords();   // retrieving vector of words for a certain algo
	  Int_t NWords = L1AlgoWords.size();
	  for (Int_t kWords=0; kWords<NWords; kWords++) {
	    UInt_t Word = L1AlgoWords[kWords];   // each word contains specific, detector-related information (eg. how many sectors fired in the KTAG, number of NewCHOD hits, etc.)
	  }
	}
      }
    }
  }
}

L1DataPacket::~L1DataPacket(){
}
