// ---------------------------------------------------------------
// History:                                                                                 
//                                                                                           
// Created by Lorenza Iacobuzio (lorenza.iacobuzio@cern.ch) February 2018            
//                       
// ---------------------------------------------------------------

#ifndef HNLWEIGHT_HH
#define HNLWEIGHT_HH

#include <stdlib.h>
#include <iostream>
#include <vector>

class HNLWeight {

public:
  
  std::vector<std::map<std::string, Double_t>> ComputeWeight(Event*, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);

private:
  
  HNLWeight();
  ~HNLWeight() {}

  std::vector<std::map<std::string, Double_t>> fWeightContainer;
};

#endif
