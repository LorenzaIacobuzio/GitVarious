Double_t GF = 1.1E-11; // MeV^-2                             
Double_t fPi = 130.41; // MeV                                   
Double_t fRho = 1.04E5; // MeV^2
Double_t Me = 0.511;
Double_t Mmu = 105.66;
Double_t Mtau = 1776.82;
Double_t Mpi = 139.57;
Double_t DMass = 1869.62;
Double_t DSMass = 1968.28;
Double_t De2BR = 9.2E-9;
Double_t Dm2BR = 3.74E-4;
Double_t DSe2BR = 1.4E-7;
Double_t DSm2BR = 5.56E-3;
Double_t DSt2BR = 5.48E-2;
Double_t TauPi2BR = 10.82E-2;
Double_t D0Mass = 1864.84;
Double_t KMass = 493.68;
Double_t K0Mass = 497.61;
Double_t Mrho = 775.4;
const int Mass = 10000;
const int step = 1;
const int Masses = Mass/step;
TString name = "";
Double_t labelSize = 0.05;
Double_t titleSize = 0.07;
Double_t PiE = 0.;
Double_t PiMu = 0.;



Double_t PhaseSpace2(Double_t Mass1, Double_t Mass2, Double_t Mass3) {

  Double_t phaseSpace = TMath::Power(Mass1*Mass1-Mass2*Mass2-Mass3*Mass3,2) - 4.*Mass2*Mass2*Mass3*Mass3;

  return phaseSpace;
}

Double_t PhaseSpaceFactor2(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t phaseSpace) {

  Double_t factor = 0.;

  if(phaseSpace > 0.) {
    factor = ((Mass1*Mass1*(Mass2*Mass2 + Mass3*Mass3)-TMath::Power(Mass2*Mass2-Mass3*Mass3,2))*TMath::Power(phaseSpace,0.5))/(Mass3*Mass3*TMath::Power(Mass1*Mass1-Mass3*Mass3,2));
  }

  return factor;
}

Double_t GammaLeptonNu3(Double_t Mass1, Double_t Mass2, Double_t Mass3) {

  Double_t r  = 4*Mass2*Mass3/(Mass1*Mass1);
  Double_t a = TMath::Power(1-r,0.5)*(1./3.-7*r/6.-r*r/24.-r*r*r/16.);
  Double_t b = r*r*(1-r*r/16.)*TMath::ATanH(TMath::Power(1-r,0.5));
  Double_t f = a+b;
  Double_t gamma_l_l_nu = 0.;

  gamma_l_l_nu = GF*GF*TMath::Power(Mass1,5)*f/(64.*TMath::Power(TMath::Pi(),3));

  return gamma_l_l_nu;
}

Double_t Gamma2(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t form) {

  Double_t gamma_2 = (GF*GF*form*form/(8*TMath::Pi()*Mass1*Mass1)*(Mass3*Mass3 + Mass1*Mass1)*(Mass1*Mass1 - Mass2*Mass2 - Mass3*Mass3)*TMath::Power(Mass1*Mass1*Mass1*Mass1 + Mass2*Mass2*Mass2*Mass2 + Mass3*Mass3*Mass3*Mass3 - 2.*Mass1*Mass3*Mass1*Mass3 - 2.*Mass1*Mass2*Mass1*Mass2 - 2.*Mass2*Mass3*Mass2*Mass3,0.5))/(2*Mass1);

  return gamma_2;
}
/*
Double_t GammaPiE(Double_t MN) {
  
  Double_t gamma_pi_e = (GF*GF*fPi*fPi/(8*TMath::Pi()*MN*MN)*(Me*Me + MN*MN)*(MN*MN - Mpi*Mpi - Me*Me)*TMath::Power(MN*MN*MN*MN + Mpi*Mpi*Mpi*Mpi + Me*Me*Me*Me - 2.*MN*Me*MN*Me - 2.*MN*Mpi*MN*Mpi - 2.*Mpi*Me*Mpi*Me,0.5))/(2*MN);

  return gamma_pi_e;
}

Double_t GammaPiMu(Double_t MN) {

  Double_t gamma_pi_mu = (GF*GF*fPi*fPi/(8*TMath::Pi()*MN*MN)*(Mmu*Mmu + MN*MN)*(MN*MN - Mpi*Mpi - Mmu*Mmu)*TMath::Power(MN*MN*MN*MN + Mpi*Mpi*Mpi*Mpi + Mmu*Mmu*Mmu*Mmu - 2.*MN*Mmu*MN*Mmu - 2.*MN*Mpi*MN*Mpi - 2.*Mpi*Mmu*Mpi*Mmu,0.5))/(2*MN);

  return gamma_pi_mu;
}
*/
Double_t GammaTot(Double_t MN) {

  Double_t USquared = 1.E-6;
  Double_t UeSquared = USquared/20.8;
  Double_t UmuSquared = 16.*UeSquared;
  Double_t UtauSquared = 3.8*UeSquared;
  Double_t Coupling = 0.;
  Double_t gammaTot = 0.;

  // 3-body decay

  if (MN > 2*Me && MN < (Me+Mmu))
    gammaTot = GammaLeptonNu3(MN, Me, Me)*UeSquared;
  else if (MN > (Me+Mmu) && MN < (Me+Mpi))
    gammaTot = GammaLeptonNu3(MN, Me, Mmu)*TMath::Sqrt(UeSquared)*TMath::Sqrt(UmuSquared);



	   
  return gammaTot;
}

Double_t Compute2BR(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t Factor, Bool_t Prod) {

  Double_t brt = 0.;

  if (Prod == kTRUE) {
    if (Factor > 0.) {
      if (Mass1 == DMass && Mass3 == Me)
	brt = De2BR*Factor;
      else if (Mass1 == DMass && Mass3 == Mmu)
	brt = Dm2BR*Factor;
      else if (Mass1 == DSMass && Mass3 == Me)
	brt = DSe2BR*Factor;
      else if (Mass1 == DSMass && Mass3 == Mmu)
	brt = DSm2BR*Factor;
      else if (Mass1 == DSMass && Mass3 == Mtau)
	brt = DSt2BR*Factor;
      else if (Mass1 == Mtau && Mass3 == Mpi) 
	brt = TauPi2BR*Factor*DSt2BR*PhaseSpaceFactor2(DSMass, 0., Mtau, PhaseSpace2(DSMass, 0., Mtau));
    }
  }
  else if (Prod == kFALSE) {}

  return brt;
}
/*
Double_t ComputeTotalBR(Double_t mesonMass, Double_t MN) {

  Double_t brt = 0.;

  Double_t phaseSpaceE = PhaseSpace2(mesonMass, MN, Me);
  Double_t phaseSpaceMu = PhaseSpace2(mesonMass, MN, Mmu);
  Double_t FactorE = PhaseSpaceFactor2(mesonMass, MN, Me, phaseSpaceE);
  Double_t FactorMu = PhaseSpaceFactor2(mesonMass, MN, Mmu, phaseSpaceMu);
  
  if (FactorE > 0. && FactorMu > 0.) {
    if (mesonMass == DMass)
      brt = De2BR*FactorE/(De2BR*FactorE+Dm2BR*FactorMu);
    
    else if (mesonMass == DSMass)
      brt = DSe2BR*FactorE/(DSe2BR*FactorE+DSm2BR*FactorMu);
  }
  
  return brt;
}
*/
void MassScan2(Double_t Mass1, Double_t Mass2, Bool_t Prod) {

  Double_t PS = 0.;
  Double_t PSF = 0.;
  Double_t BR = 0.;
  Double_t xMN[Masses];
  Double_t yBR[Masses];

  if (Prod == kTRUE) {
    for (Int_t MN = 0; MN < Mass; MN += step) {
      PS = PhaseSpace2(Mass1, MN, Mass2);
      PSF = PhaseSpaceFactor2(Mass1, MN, Mass2, PS);
      BR = Compute2BR(Mass1, MN, Mass2, PSF, Prod);
      xMN[MN] = MN;
      yBR[MN] = BR;
    }
  }

  else if (Prod == kFALSE) {}
  
  TGraph* gr = new TGraph(Mass, xMN, yBR);
  Int_t color = (Int_t)(rand() % 8 + 2);
  gr->SetLineColor(color);
  gr->SetLineWidth(3);
  gr->GetYaxis()->SetTitle("F*BR");
  gr->SetTitle("Phase space factor*BR vs N mass");
  gr->GetXaxis()->SetTitle("N mass [MeV]");
  gr->GetXaxis()->SetTitleOffset(0.9);
  gr->GetXaxis()->SetTitleSize(labelSize);
  gr->GetYaxis()->SetTitleSize(labelSize);
  gr->GetXaxis()->SetLabelSize(labelSize);
  gr->GetYaxis()->SetLabelSize(labelSize);
  gr->Draw("AC");
  gr->Write();  
}

void GeneralPlotsNew() {

  TFile *f = new TFile("/home/li/Desktop/Plots.root","RECREATE");
  Double_t yGammaPiE[Masses];
  Double_t yGammaPiMu[Masses];

  // HNL production via two-body decays
  
  MassScan2(DMass, Me, 1);
  MassScan2(DMass, Mmu, 1);
  MassScan2(DSMass, Me, 1);
  MassScan2(DSMass, Mmu, 1);
  MassScan2(DSMass, Mtau, 1);
  MassScan2(Mtau, Mpi, 1);
  /*
  for (Int_t MN = 0; MN < Mass; MN += step) {
    if (MN > Mpi + Me)
      PiE = GammaPiE(MN)*1.E13;
    else
      PiE = 0.;
    if (MN > Mpi + Mmu)
      PiMu = GammaPiMu(MN)*1.E13;
    else
      PiMu = 0.;

    xMN[MN] = MN;
    if (PiE >= 0.)
      yGammaPiE[MN] = PiE;
    if (PiMu >= 0.)
      yGammaPiMu[MN] = PiMu;
  }
  
  TGraph* gr5 = new TGraph(Mass, xMN, yGammaPiE);
  gr5->SetName("gr5");

  TGraph* gr6 = new TGraph(Mass, xMN, yGammaPiMu);
  gr6->SetName("gr6");
  */

  f->Write();
}
