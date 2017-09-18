#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "DebugFilter.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "DownstreamTrack.hh"
#include "TriggerConditions.hh"
#include "GeometricAcceptance.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

#define TRIGGER_L0_PHYSICS_TYPE 0x1
#define TRIGGER_L0_CONTROL_TYPE 0x10

/// \class DebugFilter
/// \Brief
/// Filter events for pi-mu, pi-e, e-e and mu-mu exotics triggers
/// \EndBrief

DebugFilter::DebugFilter(Core::BaseAnalysis *ba) :
  Analyzer(ba, "DebugFilter") {

  RequestAllRecoTrees();
  RequestL0Data();

  fCDAcomp = new TwoLinesCDA();
  fDistcomp = new PointLineDistance();
  AddParam("con", &fcon, 0);
}

void DebugFilter::Process(Int_t) {

  L0TPData *L0TPData = GetL0Data();
  UChar_t L0DataType = L0TPData->GetDataType();
  UInt_t L0TriggerFlags = L0TPData->GetTriggerFlags();
  Bool_t L0PhysicsTrigger = (L0DataType & TRIGGER_L0_PHYSICS_TYPE);
  Bool_t L0ControlTrigger = (L0DataType & TRIGGER_L0_CONTROL_TYPE);
  Bool_t L0OK = (L0TriggerFlags & 4) || (L0TriggerFlags & 8) || (L0TriggerFlags & 16);
  Int_t DW = 10;

  if (L0ControlTrigger) {
    fcon++;
    if ((fcon%DW) == 0) {
      FilterAccept();
      return;
    }
  }

  if (!(L0PhysicsTrigger && L0OK)) return;
  std::vector<DownstreamTrack> Tracks =
    *(std::vector<DownstreamTrack>*) GetOutput("DownstreamTrackBuilder.Output");
  if (Tracks.size()<2) return;
  for (UInt_t i = 0; i < Tracks.size(); i++) {
    if (GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[i], kSpectrometer, 0) &&
	GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[i], kSpectrometer, 1) &&
	GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[i], kSpectrometer, 2) &&
	GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[i], kSpectrometer, 3) &&
	GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[i], kCHOD)            &&
	GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[i], kMUV3)            &&
	GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[i], kLKr)) {
      for (UInt_t j = i+1; j < Tracks.size(); j++) {
	if (GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[j], kSpectrometer, 0) &&
	    GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[j], kSpectrometer, 1) &&
	    GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[j], kSpectrometer, 2) &&
	    GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[j], kSpectrometer, 3) &&
	    GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[j], kCHOD)            &&
	    GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[j], kMUV3)            &&
	    GeometricAcceptance::GetInstance()->InAcceptance(&Tracks[j], kLKr)) {
	  fCDAcomp->SetLine1Point1(Tracks[i].GetPositionBeforeMagnet());
	  fCDAcomp->SetLine2Point1(Tracks[j].GetPositionBeforeMagnet());
	  fCDAcomp->SetDir1(Tracks[i].GetMomentumBeforeMagnet());
	  fCDAcomp->SetDir2(Tracks[j].GetMomentumBeforeMagnet());
	  fCDAcomp->ComputeVertexCDA();
	  if (fCDAcomp->GetCDA() < 50.) {
	    fDistcomp->SetLinePoint1(0., 0., 102000.);
	    fDistcomp->SetLineDir(1.2e-3, 0., 1.);
	    fDistcomp->SetPoint(fCDAcomp->GetVertex());
	    fDistcomp->ComputeDistance();
	    if (fDistcomp->GetDistance() > 100.) {
	      ofstream outfile;
              outfile.open("debugfilter.txt", std::ios_base::app);
              outfile << GetEventHeader()->GetTimeStamp() << endl;
	      outfile.close();
	      return;
	    }
	  }
	}
      }
    }
  }
}

DebugFilter::~DebugFilter() {
  delete fCDAcomp;
  delete fDistcomp;
}
