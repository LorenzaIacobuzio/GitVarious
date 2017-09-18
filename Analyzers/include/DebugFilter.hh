#ifndef DEBUGFILTER_HH
#define DEBUGFILTER_HH

#include <stdlib.h>
#include "Analyzer.hh"
#include "TwoLinesCDA.hh"
#include "PointLineDistance.hh"
#include "LAVMatching.hh"
#include "SAVMatching.hh"

class TH1I;
class TH2F;
class TGraph;
class TTree;

class DebugFilter : public NA62Analysis::Analyzer
{
public:
  
  DebugFilter(NA62Analysis::Core::BaseAnalysis *ba);
  ~DebugFilter();
  void InitHist() {};
  void DefineMCSimple() {};
  void InitOutput() {};
  void PreProcess() {};
  void Process(int iEvent);
  void PostProcess() {};
  void StartOfBurstUser() {};
  void EndOfBurstUser() {};
  void StartOfRunUser() {};
  void EndOfRunUser() {};
  void EndOfJobUser() {};
  void DrawPlot() {};

protected:

  long int fcon;
  TwoLinesCDA *fCDAcomp;  
  PointLineDistance *fDistcomp;
};
#endif
