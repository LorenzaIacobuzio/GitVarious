#ifndef ACCEPTANCE_HH
#define ACCEPTANCE_HH

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


class Acceptance : public NA62Analysis::Analyzer
{

public:
  Acceptance(NA62Analysis::Core::BaseAnalysis *ba);
  ~Acceptance();
  void InitHist();
  void InitOutput() {};
  void DefineMCSimple() {};
  void Process(int);
  void StartOfBurstUser() {};
  void EndOfBurstUser() {};
  void StartOfRunUser() {};
  void EndOfRunUser() {};
  void EndOfJobUser();
  void PostProcess() {};
  void DrawPlot() {};

protected:

  TH1D *fhNk3pi;
  TH1D *fhNtot;
};
#endif
