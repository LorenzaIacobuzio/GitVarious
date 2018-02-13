// ---------------------------------------------------------------
// History:                                                                                 
//                                                                                           
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018            
//                       
// ---------------------------------------------------------------

#ifndef HEAVYNEUTRALLEPTONWEIGHT_HH
#define HEAVYNEUTRALLEPTONWEIGHT_HH

#include <stdlib.h>
#include <iostream>
#include <vector>

class HeavyNeutralLeptonWeight {

public:

  HeavyNeutralLeptonWeight();
  ~HeavyNeutralLeptonWeight() {}
  std::vector<std::map<std::string, Double_t>> ComputeWeight(Event*, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);

protected:

  std::vector<std::map<std::string, Double_t>> fWeightContainer;
};

#endif
