#ifndef FILTERMUONELECTRONTHREETRACKVERTEX.CC_HH
#define FILTERMUONELECTRONTHREETRACKVERTEX.CC_HH

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


class FilterMuonElectronThreeTrackVertex.cc : public NA62Analysis::Analyzer
{
	public:
		FilterMuonElectronThreeTrackVertex.cc(NA62Analysis::Core::BaseAnalysis *ba);
		~FilterMuonElectronThreeTrackVertex.cc();
		void InitHist();
		void InitOutput();
		void DefineMCSimple();
		void ProcessSpecialTriggerUser(int iEvent, unsigned int triggerType);
		void Process(int iEvent);
		void PostProcess();
		void StartOfBurstUser();
		void EndOfBurstUser();
		void StartOfRunUser();
		void EndOfRunUser();
        void EndOfJobUser();
		void DrawPlot();
	protected:


};
#endif
