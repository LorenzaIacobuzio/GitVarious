#ifndef HNLFUNCTIONS_HH
#define HNLFUNCTIONS_HH

#include <stdlib.h>
#include <iostream>
#include <string>
#include <cmath>
#include <TROOT.h>
#include <TChain.h>
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/GaussLegendreIntegrator.h"
#include "TF1.h"
#include "TF2.h"
#include "Rtypes.h"

extern Double_t ComputeHNLMass(KinePart*);
extern Double_t ComputeL(TVector3, TVector3, TVector3);
extern Double_t ComputeNDecayProb(KinePart*, Double_t, Double_t);
extern Double_t ComputeNReachProb(KinePart*, Double_t, Double_t);
extern Double_t PhaseSpace(Double_t, Double_t, Double_t);
extern Double_t PhaseSpaceFactor(Double_t, Double_t, Double_t);
extern Double_t TwoBodyBR(Double_t, Double_t, Double_t, Int_t, Bool_t);
extern Double_t ThreeBodyBR(Double_t, Double_t, Double_t, Double_t, Int_t, Bool_t);
extern std::string ThreeBodyFunction(Double_t, Double_t);
extern Double_t Gamma2(Double_t, Double_t, Double_t, Double_t, Bool_t);
extern Double_t GammaLeptonNu3(Double_t, Double_t, Double_t, Bool_t);
extern Double_t GammaTot(Double_t);
extern Double_t tauN(Double_t);
extern Double_t lambda(Double_t, Double_t, Double_t);
extern Double_t ComputeProd(KinePart*, Double_t);
extern Double_t ComputeDecay(Double_t);

// Masses                                                                                              

extern Double_t fMe;
extern Double_t fMmu;
extern Double_t fMtau;
extern Double_t fMpi;
extern Double_t fMpi0;
extern Double_t fMrho;
extern Double_t fMrho0;
extern Double_t fMeta;
extern Double_t fMetaprime;
extern Double_t fMD;
extern Double_t fMDS;
extern Double_t fMD0;
extern Double_t fMK;
extern Double_t fMK0;
extern Double_t fMp;
extern Double_t fMKStar;
extern Double_t fMK0Star;

// Lifetimes                                                                                          

extern Double_t fDlife;
extern Double_t fDSlife;
extern Double_t fD0life;
extern Double_t ftaulife;

// Constants                                                                                          

extern Double_t fhc;
extern Double_t fcLight;
extern Double_t fGF;
extern Double_t fPi;
extern Double_t fRho;
extern Double_t fD;
extern Double_t fDS;
extern Double_t fK;
extern Double_t fEta;
extern Double_t fEtaprime;
extern Double_t fsigmacc;

// CKM                                                                                                 

extern Double_t fVcs;
extern Double_t fVcd;
extern Double_t fVud;
extern Double_t fVus;

// Form factors, pseudoscalar and vector mesons                                                       

extern Double_t fDK0;
extern Double_t fDpi0;
extern Double_t fD0K;
extern Double_t fD0pi;
extern Double_t fgDK0;
extern Double_t fgDpi0;
extern Double_t fgD0K;
extern Double_t fgD0pi;
extern Double_t fA0D;
extern Double_t fA1D;
extern Double_t fA2D;
extern Double_t fVD;
extern Double_t fA0D0;
extern Double_t fA1D0;
extern Double_t fA2D0;
extern Double_t fVD0;

// Fragmentation fractions                                                                            

extern Double_t ffD;
extern Double_t ffD0;
extern Double_t ffDS;

// Other parameters

extern Double_t fUSquared;
extern Double_t fUeSquared;
extern Double_t fUmuSquared;
extern Double_t fUtauSquared;
extern Double_t fUeSquaredRatio;
extern Double_t fUmuSquaredRatio;
extern Double_t fUtauSquaredRatio;

#endif // HNLFUNCTIONS_HH
