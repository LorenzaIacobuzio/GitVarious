// ---------------------------------------------------------------                                    
// History:                                                                                           
//                                                                                                    
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018                             
//                                                                                                    
// --------------------------------------------------------------- 

#ifndef HEAVYNEUTRINO_HH
#define HEAVYNEUTRINO_HH

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <TCanvas.h>
#include "Analyzer.hh"
#include "TwoLinesCDA.hh"
#include "PointLineDistance.hh"
#include "LAVMatching.hh"
#include "SAVMatching.hh"

class TH1I;
class TH2F;
class TGraph;
class TTree;

class HeavyNeutrino : public NA62Analysis::Analyzer {

public:

  HeavyNeutrino(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutrino();
  void InitHist();
  void InitOutput();
  void DefineMCSimple() {}
  void Process(Int_t);
  void StartOfBurstUser() {}
  void EndOfBurstUser();
  void StartOfRunUser() {}
  void EndOfRunUser() {}
  void EndOfJobUser();
  void PostProcess() {}
  void DrawPlot() {}

protected:

  TwoLinesCDA *fCDAcomp;
  PointLineDistance *fDistcomp;
  LAVMatching *fLAVMatching;
  SAVMatching *fSAVMatching;
  std::vector<Int_t> fID;
  std::vector<std::string> fStream;
  Bool_t fPassSelection;
  Double_t fInitialFV;
  Double_t fLFV;

  // Histos

  TH1D *fhNk3pi;
  TH1D *fhNbursts;
  TH1D *fhNEvents;
  TH1D *fhN2tracks;
  TH1D *fhNtracks;
  TH1D *fhMomPi;
  TH1D *fhMomMu;

  TH2D *fhXYSpec0Reco;
  TH2D *fhXYSpec1Reco;
  TH2D *fhXYSpec2Reco;
  TH2D *fhXYSpec3Reco;
  TH2D *fhXYCHODReco;
  TH2D *fhXYCHODTrue;
  TH2D *fhXYMUV3True;

  TH1D *fhCuts;

  TH2D *fhCDAvsZ_In;
  TH2D *fhCDAvsZ_Track;
  TH2D *fhCDAvsZ_Energy;
  TH2D *fhCDAvsZ_Geom;
  TH2D *fhCDAvsZ_Vetoes;
  TH2D *fhCDAvsZ_Fin;

  TH2D *fhZvsBeam_In;
  TH2D *fhZvsBeam_Track;
  TH2D *fhZvsBeam_Energy;
  TH2D *fhZvsBeam_Geom;
  TH2D *fhZvsBeam_Vetoes;
  TH2D *fhZvsBeam_Fin;

  TH2D *fhBeamvsTar_In;
  TH2D *fhBeamvsTar_Track;
  TH2D *fhBeamvsTar_Energy;
  TH2D *fhBeamvsTar_Geom;
  TH2D *fhBeamvsTar_Vetoes;
  TH2D *fhBeamvsTar_Fin;

  TH1D *fhNMUV3Cand;
  TH1D *fhEoP;
  TH2D *fhEoPMuVsPi;
  TH1D *fhInvMassReco;
};

#endif
