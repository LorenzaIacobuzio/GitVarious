#ifndef L1EfficiencyKTAG_HH
#define L1EfficiencyKTAG_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include "DetectorAcceptance.hh"
#include <TCanvas.h>

class TH1I;
class TH2F;
class TGraph;
class TTree;

class L1EfficiencyKTAG : public NA62Analysis::Analyzer
{
public:
  L1EfficiencyKTAG(NA62Analysis::Core::BaseAnalysis *ba);
  ~L1EfficiencyKTAG();
  void InitHist();
  void InitOutput() {}
  void DefineMCSimple() {}
  void Process(Int_t);
  void StartOfBurstUser() {}
  void EndOfBurstUser() {}
  void StartOfRunUser() {}
  void EndOfRunUser() {}
  void EndOfJobUser();
  void PostProcess() {}
  void DrawPlot() {}
  void PrintStatisticsPerBurst() {}
  
private:
  Int_t fCedarDen;
  Int_t fCedarNum;
};
#endif
