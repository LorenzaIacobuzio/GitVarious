// ---------------------------------------------------------------                                    
// History:                                                                                           
//                                                                                                    
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018                             
//                                                                                                    
// --------------------------------------------------------------- 

#ifndef HEAVYNEUTRINONEG_HH
#define HEAVYNEUTRINONEG_HH

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <TCanvas.h>
#include "Analyzer.hh"
#include "TwoLinesCDA.hh"
#include "PointLineDistance.hh"
#include "LAVMatching.hh"
#include "SAVMatching.hh"
#include "TriggerConditions.hh"

class TH1I;
class TH2F;
class TGraph;
class TTree;

class HeavyNeutrinoNeg : public NA62Analysis::Analyzer {

public:

  HeavyNeutrinoNeg(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrinoNeg() {}
  void InitHist();
  void InitOutput();
  void DefineMCSimple() {}
  void Process(Int_t);
  void StartOfBurstUser() {}
  void EndOfBurstUser();
  void StartOfRunUser() {}
  void EndOfRunUser() {}
  void EndOfJobUser();
  void PostProcess();
  void DrawPlot() {}
  void ProcessEOBEvent() {}

protected:

  TwoLinesCDA *fCDAcomp;
  PointLineDistance *fDistcomp;
  LAVMatching *fLAVMatching;
  SAVMatching *fSAVMatching;
  std::vector<Int_t> fID;
  std::vector<std::string> fStream;
  Bool_t fPassSelection;
  Bool_t fReadingData;
  Bool_t fBlindRegion;
  Bool_t fMCsample;
  Bool_t fEnableChecks;
  Int_t fMode;
  Int_t fBurstCounter;
  Int_t fNK3Pi;
  Double_t fInitialFV;
  Double_t fLFV;
  Double_t fMassForReco;
  Double_t fMassForChecks;
  Double_t fCouplingForChecks;
  Double_t fzCHOD;
  Double_t fzMUV3;
  Double_t fzStraw[4];
  Double_t fxStrawChamberCentre[4];
  Double_t frMinStraw;
  Double_t frMaxStraw;
  Double_t fzCHODPlane;
  Double_t frMinCHOD;
  Double_t frMaxCHOD;
  Double_t fZGTK3;
  TRecoCedarCandidate *KTAGcand;

  // Variables for TTrees

  Double_t Weight;
  Double_t CHODTime1;
  Double_t CHODTime2;
  Double_t CDA;
  Double_t Zvertex;
  Double_t CDALine;
  Double_t ZCDALine;
  Double_t BeamlineDist;
  Double_t xSR;
  Double_t ySR;
  Double_t MuEoP;
  Double_t PiEoP;
  Double_t R;
  Double_t energyPi;
  Double_t energyMu;
  Double_t invMass;
  Double_t trueMass;
  Double_t L0TPTime;
  Double_t xGTK31;
  Double_t yGTK31;
  Double_t xGTK32;
  Double_t yGTK32;
  Double_t BeamCDA1;
  Double_t BeamCDA2;
  Double_t EtotLKr;
  TVector3 Mom1;
  TVector3 Mom2;
  TVector3 threeMomPi;
  TVector3 threeMomMu;
  TVector3 TotMom;
  TVector3 Vertex;
  TVector3 Pos1;
  TVector3 Pos2;
  TVector3 BeamVtx1;
  TVector3 BeamVtx2;
  Bool_t Target;
  Bool_t K3pi;
  Bool_t autoPass;
  Bool_t CHANTIAssoc1;
  Bool_t CHANTIAssoc2;
  Bool_t isControlTrigger;
  Int_t Assoc;
  Int_t Charge1;
  Int_t Charge2;
  Int_t nSec;
  Int_t nCHOD;
  Int_t nNewCHOD;
  Int_t nLKr;
  std::vector<Double_t> RICHInfo1;
  std::vector<Double_t> RICHInfo2;
};

#endif
