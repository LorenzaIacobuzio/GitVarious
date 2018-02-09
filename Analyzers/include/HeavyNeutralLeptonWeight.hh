#ifndef HEAVYNEUTRALLEPTONWEIGHT_HH
#define HEAVYNEUTRALLEPTONWEIGHT_HH

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <TCanvas.h>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include "MCInfo.hh"

class TH1I;
class TH2F;
class TGraph;
class TTree;

class HeavyNeutralLeptonWeight : public NA62Analysis::Analyzer {

public:

  HeavyNeutralLeptonWeight(NA62Analysis::Core::BaseAnalysis *ba);
  ~HeavyNeutralLeptonWeight() {}
  void InitHist() {}
  void InitOutput();
  void DefineMCSimple() {}
  void Process(Int_t);
  void StartOfBurstUser() {}
  void EndOfBurstUser() {}
  void StartOfRunUser() {}
  void EndOfRunUser() {}
  void EndOfJobUser() {} 
  void PostProcess() {}
  void DrawPlot() {}

protected:

  std::vector<std::map<std::string, Double_t>> fWeightContainer;
  Double_t fLInitialFV;
  Double_t fLFV;
};

#endif
