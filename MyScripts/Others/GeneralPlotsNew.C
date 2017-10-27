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
const int InitialMass = 100;
const int Mass = 100000;
const int step = 1;
const int Masses = Mass/step;
TString name = "";
Double_t labelSize = 0.05;
Double_t titleSize = 0.07;
Double_t PiE = 0.;
Double_t PiMu = 0.;



Double_t PhaseSpace(Double_t Mass1, Double_t Mass2, Double_t Mass3, Bool_t TwoBody) {

  Double_t phaseSpace = 0.;

  if (TwoBody == kTRUE)
    phaseSpace = TMath::Power(Mass1*Mass1-Mass2*Mass2-Mass3*Mass3,2) - 4.*Mass2*Mass2*Mass3*Mass3;

  return phaseSpace;
}

Double_t PhaseSpaceFactor(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t phaseSpace, Bool_t TwoBody) {

  Double_t factor = 0.;

  if (TwoBody == kTRUE) {
    if(phaseSpace > 0.) {
      factor = (Mass1*Mass1*(Mass2*Mass2 + Mass3*Mass3)-TMath::Power(Mass2*Mass2-Mass3*Mass3,2))*TMath::Power(phaseSpace,0.5)/(Mass3*Mass3*TMath::Power(Mass1*Mass1-Mass3*Mass3,2));
    }
  }

  return factor;
}

Double_t GammaLeptonNu3(Double_t Mass1, Double_t Mass2, Double_t Mass3) {

  Double_t r  = 4*Mass2*Mass3/(Mass1*Mass1);
  Double_t a = TMath::Power(1-r,0.5)*(1./3.-7*r/6.-r*r/24.-r*r*r/16.);
  Double_t b = r*r*(1-r*r/16.)*TMath::ATanH(TMath::Power(1-r,0.5));
  Double_t f = a+b;
  Double_t gamma_l_l_nu = 0.;

  if (Mass1 >= (Mass2+Mass3))
    gamma_l_l_nu = GF*GF*TMath::Power(Mass1,5)*f/(64.*TMath::Power(TMath::Pi(),3));

  return gamma_l_l_nu;
}

Double_t Gamma2(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t form) {

  Double_t gamma_2 = 0.;
  
  if (Mass1 >= (Mass2+Mass3)) 
    gamma_2 = (GF*GF*form*form/(8*TMath::Pi()*Mass1*Mass1)*(Mass3*Mass3 + Mass1*Mass1)*(Mass1*Mass1 - Mass2*Mass2 - Mass3*Mass3)*TMath::Power(Mass1*Mass1*Mass1*Mass1 + Mass2*Mass2*Mass2*Mass2 + Mass3*Mass3*Mass3*Mass3 - 2.*Mass1*Mass3*Mass1*Mass3 - 2.*Mass1*Mass2*Mass1*Mass2 - 2.*Mass2*Mass3*Mass2*Mass3,0.5))/(2*Mass1);

  return gamma_2;
}

Double_t GammaTot(Double_t MN) {
  /*
  Double_t USquared = 1.E-6;
  Double_t UeSquared = USquared/20.8;
  Double_t UmuSquared = 16.*UeSquared;
  Double_t UtauSquared = 3.8*UeSquared;
  */
  Double_t gammaTot = 0.;

  if (MN < 2*Me) 
    gammaTot = 0.;
  else if (MN >= 2*Me && MN < (Me+Mmu)) 
    gammaTot = GammaLeptonNu3(MN, Me, Me)/**UeSquared*/;
  else if (MN >= (Me+Mmu) && MN < (Me+Mpi)) 
    gammaTot = GammaLeptonNu3(MN, Me, Me) + GammaLeptonNu3(MN, Me, Mmu)/**TMath::Sqrt(UeSquared)*TMath::Sqrt(UmuSquared)*/;
  else if (MN >= (Me+Mpi) && MN < 2*Mmu) 
    gammaTot = GammaLeptonNu3(MN, Me, Me) + GammaLeptonNu3(MN, Me, Mmu) + Gamma2(MN, Me, Mpi, fPi)/**UeSquared*/;
  else if (MN >= 2*Mmu && MN < (Mmu+Mpi)) 
    gammaTot = GammaLeptonNu3(MN, Me, Me) + GammaLeptonNu3(MN, Me, Mmu) + Gamma2(MN, Me, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mmu)/**UmuSquared*/;
  else if (MN >= (Mmu+Mpi) && MN < (Mrho+Me)) 
    gammaTot = GammaLeptonNu3(MN, Me, Me) + GammaLeptonNu3(MN, Me, Mmu) + Gamma2(MN, Me, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mmu) + Gamma2(MN, Mmu, Mpi, fPi)/**UmuSquared*/;
  else if (MN >= (Mrho+Me) && MN < (Mrho+Mmu)) 
    gammaTot = GammaLeptonNu3(MN, Me, Me) + GammaLeptonNu3(MN, Me, Mmu) + Gamma2(MN, Me, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mmu) + Gamma2(MN, Mmu, Mpi, fPi) + Gamma2(MN, Mrho, Me, fRho)/**UeSquared*/;
  else if (MN >= (Mmu+Mrho) && MN < (Me+Mtau)) 
    gammaTot = GammaLeptonNu3(MN, Me, Me) + GammaLeptonNu3(MN, Me, Mmu) + Gamma2(MN, Me, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mmu) + Gamma2(MN, Mmu, Mpi, fPi) + Gamma2(MN, Mrho, Me, fRho) + Gamma2(MN, Mrho, Mmu, fRho)/**UmuSquared*/;
  else if (MN >= (Me+Mtau) && MN < (Mpi+Mtau)) 
    gammaTot = GammaLeptonNu3(MN, Me, Me) + GammaLeptonNu3(MN, Me, Mmu) + Gamma2(MN, Me, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mmu) + Gamma2(MN, Mmu, Mpi, fPi) + Gamma2(MN, Mrho, Me, fRho) + Gamma2(MN, Mrho, Mmu, fRho) + GammaLeptonNu3(MN, Me, Mtau)/**TMath::Sqrt(UeSquared)*TMath::Sqrt(UtauSquared)*/;
  else if (MN >= (Mpi+Mtau) && MN < (Mmu+Mtau)) 
    gammaTot = GammaLeptonNu3(MN, Me, Me) + GammaLeptonNu3(MN, Me, Mmu) + Gamma2(MN, Me, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mmu) + Gamma2(MN, Mmu, Mpi, fPi) + Gamma2(MN, Mrho, Me, fRho) + Gamma2(MN, Mrho, Mmu, fRho) + GammaLeptonNu3(MN, Me, Mtau) + Gamma2(MN, Mtau, Mpi, fPi)/**UtauSquared*/;
  else if (MN >= (Mmu+Mtau) && MN < (Mrho+Mtau)) 
    gammaTot = GammaLeptonNu3(MN, Me, Me) + GammaLeptonNu3(MN, Me, Mmu) + Gamma2(MN, Me, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mmu) + Gamma2(MN, Mmu, Mpi, fPi) + Gamma2(MN, Mrho, Me, fRho) + Gamma2(MN, Mrho, Mmu, fRho) + GammaLeptonNu3(MN, Me, Mtau) + Gamma2(MN, Mtau, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mtau)/**TMath::Sqrt(UmuSquared)*TMath::Sqrt(UtauSquared)*/;
  else if (MN >= (Mrho+Mtau) && MN < 2*Mtau) 
    gammaTot = GammaLeptonNu3(MN, Me, Me) + GammaLeptonNu3(MN, Me, Mmu) + Gamma2(MN, Me, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mmu) + Gamma2(MN, Mmu, Mpi, fPi) + Gamma2(MN, Mrho, Me, fRho) + Gamma2(MN, Mrho, Mmu, fRho) + GammaLeptonNu3(MN, Me, Mtau) + Gamma2(MN, Mtau, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mtau) + Gamma2(MN, Mrho, Mtau, fRho)/**UtauSquared*/;
  else if (MN >= 2*Mtau) 
    gammaTot = GammaLeptonNu3(MN, Me, Me) + GammaLeptonNu3(MN, Me, Mmu) + Gamma2(MN, Me, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mmu) + Gamma2(MN, Mmu, Mpi, fPi) + Gamma2(MN, Mrho, Me, fRho) + Gamma2(MN, Mrho, Mmu, fRho) + GammaLeptonNu3(MN, Me, Mtau) + Gamma2(MN, Mtau, Mpi, fPi) + GammaLeptonNu3(MN, Mmu, Mtau) + Gamma2(MN, Mrho, Mtau, fRho) + GammaLeptonNu3(MN, Mtau, Mtau)/**UtauSquared*/;

  return gammaTot;
}

Double_t ComputeBR(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t Factor, Bool_t Prod, Bool_t TwoBody) {

  Double_t brt = 0.;
  
  if (Prod == kTRUE) {
    if (Factor > 0. && Mass1 >= (Mass2+Mass3)) {
      if (TwoBody == kTRUE) {
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
	  brt = TauPi2BR*Factor*DSt2BR*PhaseSpaceFactor(DSMass, 0., Mtau, PhaseSpace(DSMass, 0., Mtau, TwoBody), TwoBody);
      }
      else if (TwoBody == kFALSE)
	brt = 0.;
    }
    else
      brt = 0.;
  }      
  else if (Prod == kFALSE) {
    if (TwoBody == kTRUE) {
      if ((Mass2 == Mpi || Mass3 == Mpi) && Gamma2(Mass1, Mass2, Mass3, fPi) > 0.)
	brt = Gamma2(Mass1, Mass2, Mass3, fPi)/GammaTot(Mass1);
      else if ((Mass2 == Mpi || Mass3 == Mpi) && Gamma2(Mass1, Mass2, Mass3, fPi) == 0.)
	brt = 0.;
      if ((Mass2 == Mrho || Mass3 == Mrho) && Gamma2(Mass1, Mass2, Mass3, fRho) > 0.)
	brt = Gamma2(Mass1, Mass2, Mass3, fRho)/GammaTot(Mass1);
      else if ((Mass2 == Mrho || Mass3 == Mrho) && Gamma2(Mass1, Mass2, Mass3, fRho) == 0.)
	brt = 0.;
    }
    else if (TwoBody == kFALSE) {
      if (GammaLeptonNu3(Mass1, Mass2, Mass3) > 0. && GammaTot(Mass1) > 0.)
	brt = GammaLeptonNu3(Mass1, Mass2, Mass3)/GammaTot(Mass1);
      else if (GammaLeptonNu3(Mass1, Mass2, Mass3) == 0. || GammaTot(Mass1) == 0.)
	brt = 0.;
    }
  }

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
void MassScan(Double_t Mass1, Double_t Mass2, Bool_t Prod, Bool_t TwoBody) {

  Double_t PS = 0.;
  Double_t PSF = 0.;
  Double_t BR = 0.;
  Double_t xMN[Masses];
  Double_t yBR[Masses];

  if (Prod == kTRUE) {
    for (Int_t MN = InitialMass; MN < Mass; MN += step) {
      PS = PhaseSpace(Mass1, MN, Mass2, TwoBody);
      PSF = PhaseSpaceFactor(Mass1, MN, Mass2, PS, TwoBody);
      BR = ComputeBR(Mass1, MN, Mass2, PSF, Prod, TwoBody);
      xMN[MN] = MN/1000.;
      yBR[MN] = BR;
    }
  }

  else if (Prod == kFALSE) {
    for (Int_t MN = InitialMass; MN < Mass; MN += step) {
      BR = ComputeBR(MN, Mass1, Mass2, PSF, Prod, TwoBody);
      xMN[MN] = MN/1000.;
      yBR[MN] = BR;
    }
  }
    
  TGraph* gr = new TGraph(Mass, xMN, yBR);
  Int_t color = (Int_t)(rand() % 8 + 2);
  gr->SetLineColor(color);
  gr->SetLineWidth(3);

  if (Prod == kTRUE) {
    gr->GetYaxis()->SetTitle("F*BR(D->NX)");
    gr->SetTitle("N production BR vs N mass");
  }
  else {
    gr->GetYaxis()->SetTitle("BR(N->X)");
    gr->SetTitle("N decay BR vs N mass");
  }

  gr->GetXaxis()->SetTitle("N mass [GeV]");
  //gr->GetXaxis()->SetTitleOffset(0.9);
  //gr->GetXaxis()->SetTitleSize(labelSize);
  //gr->GetYaxis()->SetTitleSize(labelSize);
  //gr->GetXaxis()->SetLabelSize(labelSize);
  //gr->GetYaxis()->SetLabelSize(labelSize);
  TGaxis::SetMaxDigits(2);
  gr->Draw("AC");
  gPad->SetLogx();
  gPad->SetLogy();
  gr->SetMinimum(1.E-15);
  gr->SetMaximum(1.5);
  gPad->Update();
  gPad->Modified();
  gPad->Write();
  gr->Write();  
}

void GeneralPlotsNew() {

  TFile *f = new TFile("/home/li/Desktop/Plots.root","RECREATE");

  // HNL production via two-body decay

  MassScan(DMass, Me, 1, 1);
  MassScan(DMass, Mmu, 1, 1);
  MassScan(DSMass, Me, 1, 1);
  MassScan(DSMass, Mmu, 1, 1);
  MassScan(DSMass, Mtau, 1, 1);
  MassScan(Mtau, Mpi, 1, 1);

  // HNL decay via two-body decay

  MassScan(Me, Me, 0, 0);
  MassScan(Me, Mmu, 0, 0);
  MassScan(Mpi, Me, 0, 1);
  MassScan(Mmu, Mmu, 0, 0);
  MassScan(Mpi, Mmu, 0, 1);
  MassScan(Mrho, Me, 0, 1);
  MassScan(Mrho, Mmu, 0, 1);
  MassScan(Me, Mtau, 0, 0);
  MassScan(Mpi, Mtau, 0, 1);
  MassScan(Mmu, Mtau, 0, 0);
  MassScan(Mrho, Mtau, 0, 1);
  MassScan(Mtau, Mtau, 0, 0);

  f->Write();
}
