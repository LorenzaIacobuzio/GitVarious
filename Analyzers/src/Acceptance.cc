#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "Acceptance.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "K3piSelection.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

Acceptance::Acceptance(Core::BaseAnalysis *ba) : Analyzer(ba, "Acceptance") {
}

void Acceptance::InitHist() {

  BookHisto("hNk3pi", new TH1D("Nk3pi", "bla", 1, 0., 1.));
  BookHisto("hNtot", new TH1D("Ntot", "bla", 1, 0., 1.));
}

void Acceptance::Process(int){

  FillHisto("hNtot", 0.5);

  Bool_t k3pi = *(Bool_t*)GetOutput("K3piSelection.EventSelected");

  if (k3pi)
    FillHisto("hNk3pi", 0.5);
}

void Acceptance::EndOfJobUser(){

  fhNk3pi = (TH1D*)fHisto.GetTH1("hNk3pi");
  fhNtot = (TH1D*)fHisto.GetTH1("hNtot");

  SaveAllPlots();
}

Acceptance::~Acceptance(){
}
