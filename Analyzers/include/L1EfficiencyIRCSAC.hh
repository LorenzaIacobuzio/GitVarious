#ifndef L1EFFICIENCYIRCSAC_HH
#define L1EFFICIENCYIRCSAC_HH

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

class L1EfficiencyIRCSAC : public NA62Analysis::Analyzer
{
public:
  L1EfficiencyIRCSAC(NA62Analysis::Core::BaseAnalysis *ba);
  ~L1EfficiencyIRCSAC();
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
  Int_t fIRCSACDen;
  Int_t fIRCSACNum;
  Int_t fIRCSACRFDen;
  Int_t fIRCSACRFNum;
  Int_t fIRCSACRVDen;
  Int_t fIRCSACRVNum;
};
#endif
