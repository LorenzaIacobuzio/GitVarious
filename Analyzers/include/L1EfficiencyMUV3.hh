#ifndef L1EFFICIENCYMUV3_HH
#define L1EFFICIENCYMUV3_HH

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

class L1EfficiencyMUV3 : public NA62Analysis::Analyzer
{
public:
  L1EfficiencyMUV3(NA62Analysis::Core::BaseAnalysis *ba);
  ~L1EfficiencyMUV3();
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
  Int_t fMUV3Den;
  Int_t fMUV3Num;
  Int_t fMUV3RFDen;
  Int_t fMUV3RFNum;
  Int_t fMUV3RVDen;
  Int_t fMUV3RVNum;
};
#endif
