#ifndef L1DATAPACKET_HH
#define L1DATAPACKET_HH

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

class L1DataPacket : public NA62Analysis::Analyzer
{

public:
  
  L1DataPacket(NA62Analysis::Core::BaseAnalysis *ba);
  ~L1DataPacket();
  void InitHist() {} 
  void InitOutput() {}
  void DefineMCSimple() {}
  void Process(int iEvent);
  void StartOfBurstUser() {}
  void EndOfBurstUser() {}
  void StartOfRunUser() {}
  void EndOfRunUser() {}
  void PostProcess() {}
  void DrawPlot() {}

protected:

};
#endif
