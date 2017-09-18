#ifndef L1EFFICIENCYNEWCHOD_HH
#define L1EFFICIENCYNEWCHOD_HH

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

class L1EfficiencyNewCHOD : public NA62Analysis::Analyzer
{
public:
  L1EfficiencyNewCHOD(NA62Analysis::Core::BaseAnalysis *ba);
  ~L1EfficiencyNewCHOD();
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
  Int_t fNewCHODDen;
  Int_t fNewCHODNum;
  Int_t fNewCHODRFDen;
  Int_t fNewCHODRFNum;
  Int_t fNewCHODRVDen;
  Int_t fNewCHODRVNum;
};
#endif
