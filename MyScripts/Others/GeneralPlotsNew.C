Double_t GF = 1.1E-11; // MeV^-2                             
Double_t fPi = 130.41; // MeV                                   
Double_t fRho = 1.04E5; // MeV^2
Double_t e = 0.511;
Double_t mu = 105.66;
Double_t tau = 1776.82;
Double_t pi = 139.57;
Double_t rho = 775.4;
Double_t D = 1869.62;
Double_t DS = 1968.28;
Double_t De2BR = 9.2E-9;
Double_t Dm2BR = 3.74E-4;
Double_t DSe2BR = 1.4E-7;
Double_t DSm2BR = 5.56E-3;
Double_t DSt2BR = 5.48E-2;
Double_t TauPi2BR = 10.82E-2;
Double_t D0 = 1864.84;
Double_t K = 493.68;
Double_t K0 = 497.61;
const int InitialMass = 100;
const int Mass = 10000;
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
  else if (TwoBody == kFALSE)
    phaseSpace = 0.;
  
  return phaseSpace;
}

Double_t PhaseSpaceFactor(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t phaseSpace, Bool_t TwoBody) {

  Double_t factor = 0.;

  if (TwoBody == kTRUE) {
    if(phaseSpace > 0.) {
      factor = (Mass1*Mass1*(Mass2*Mass2 + Mass3*Mass3)-TMath::Power(Mass2*Mass2-Mass3*Mass3,2))*TMath::Power(phaseSpace,0.5)/(Mass3*Mass3*TMath::Power(Mass1*Mass1-Mass3*Mass3,2));
    }
  }
  else if (TwoBody == kFALSE) {
    if(phaseSpace > 0.) {
      factor = 0.;
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

  if (MN < 2*e) 
    gammaTot = 0.;
  else if (MN >= 2*e && MN < (e+mu)) 
    gammaTot = GammaLeptonNu3(MN, e, e)/**UeSquared*/;
  else if (MN >= (e+mu) && MN < (e+pi)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu)/**TMath::Sqrt(UeSquared)*TMath::Sqrt(UmuSquared)*/;
  else if (MN >= (e+pi) && MN < 2*mu) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, e, pi, fPi)/**UeSquared*/;
  else if (MN >= 2*mu && MN < (mu+pi)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, e, pi, fPi) + GammaLeptonNu3(MN, mu, mu)/**UmuSquared*/;
  else if (MN >= (mu+pi) && MN < (rho+e)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, e, pi, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, mu, pi, fPi)/**UmuSquared*/;
  else if (MN >= (rho+e) && MN < (rho+mu)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, e, pi, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, mu, pi, fPi) + Gamma2(MN, rho, e, fRho)/**UeSquared*/;
  else if (MN >= (mu+rho) && MN < (e+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, e, pi, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, mu, pi, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho)/**UmuSquared*/;
  else if (MN >= (e+tau) && MN < (pi+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, e, pi, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, mu, pi, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + GammaLeptonNu3(MN, e, tau)/**TMath::Sqrt(UeSquared)*TMath::Sqrt(UtauSquared)*/;
  else if (MN >= (pi+tau) && MN < (mu+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, e, pi, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, mu, pi, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + GammaLeptonNu3(MN, e, tau) + Gamma2(MN, tau, pi, fPi)/**UtauSquared*/;
  else if (MN >= (mu+tau) && MN < (rho+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, e, pi, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, mu, pi, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + GammaLeptonNu3(MN, e, tau) + Gamma2(MN, tau, pi, fPi) + GammaLeptonNu3(MN, mu, tau)/**TMath::Sqrt(UmuSquared)*TMath::Sqrt(UtauSquared)*/;
  else if (MN >= (rho+tau) && MN < 2*tau) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, e, pi, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, mu, pi, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + GammaLeptonNu3(MN, e, tau) + Gamma2(MN, tau, pi, fPi) + GammaLeptonNu3(MN, mu, tau) + Gamma2(MN, rho, tau, fRho)/**UtauSquared*/;
  else if (MN >= 2*tau) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, e, pi, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, mu, pi, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + GammaLeptonNu3(MN, e, tau) + Gamma2(MN, tau, pi, fPi) + GammaLeptonNu3(MN, mu, tau) + Gamma2(MN, rho, tau, fRho) + GammaLeptonNu3(MN, tau, tau)/**UtauSquared*/;

  return gammaTot;
}

Double_t ComputeBR(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t Factor, Bool_t Prod, Bool_t TwoBody) {

  Double_t brt = 0.;
  
  if (Prod == kTRUE) {
    if (Factor > 0. && Mass1 >= (Mass2+Mass3)) {
      if (TwoBody == kTRUE) {
	if (Mass1 == D && Mass3 == e)
	  brt = De2BR*Factor;
	else if (Mass1 == D && Mass3 == mu)
	  brt = Dm2BR*Factor;
	else if (Mass1 == DS && Mass3 == e)
	  brt = DSe2BR*Factor;
	else if (Mass1 == DS && Mass3 == mu)
	  brt = DSm2BR*Factor;
	else if (Mass1 == DS && Mass3 == tau)
	  brt = DSt2BR*Factor;
	else if (Mass1 == tau && Mass3 == pi)
	  brt = TauPi2BR*Factor*DSt2BR*PhaseSpaceFactor(DS, 0., tau, PhaseSpace(DS, 0., tau, TwoBody), TwoBody);
      }
      else if (TwoBody == kFALSE)
	brt = 0.;
    }
    else
      brt = 0.;
  }      
  else if (Prod == kFALSE) {
    if (TwoBody == kTRUE) {
      if ((Mass2 == pi || Mass3 == pi) && Gamma2(Mass1, Mass2, Mass3, fPi) > 0.)
	brt = Gamma2(Mass1, Mass2, Mass3, fPi)/GammaTot(Mass1);
      else if ((Mass2 == pi || Mass3 == pi) && Gamma2(Mass1, Mass2, Mass3, fPi) == 0.)
	brt = 0.;
      if ((Mass2 == rho || Mass3 == rho) && Gamma2(Mass1, Mass2, Mass3, fRho) > 0.)
	brt = Gamma2(Mass1, Mass2, Mass3, fRho)/GammaTot(Mass1);
      else if ((Mass2 == rho || Mass3 == rho) && Gamma2(Mass1, Mass2, Mass3, fRho) == 0.)
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

  Double_t phaseSpaceE = PhaseSpace2(mesonMass, MN, e);
  Double_t phaseSpaceMu = PhaseSpace2(mesonMass, MN, mu);
  Double_t FactorE = PhaseSpaceFactor2(mesonMass, MN, e, phaseSpaceE);
  Double_t FactorMu = PhaseSpaceFactor2(mesonMass, MN, mu, phaseSpaceMu);
  
  if (FactorE > 0. && FactorMu > 0.) {
    if (mesonMass == D)
      brt = De2BR*FactorE/(De2BR*FactorE+Dm2BR*FactorMu);
    
    else if (mesonMass == DS)
      brt = DSe2BR*FactorE/(DSe2BR*FactorE+DSm2BR*FactorMu);
  }
  
  return brt;
}
*/
void MassScan(Double_t Mass1, Double_t Mass2, Bool_t Prod, Bool_t TwoBody, std::string Title) {

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

  TCanvas *c1 = new TCanvas();
  TGraph* gr = new TGraph(Mass, xMN, yBR);
  Int_t color = (Int_t)(rand() % 8 + 2);
  
  gr->SetLineColor(color);
  gr->SetLineWidth(3);

  if (Prod == kTRUE) {
    string title = string("F*BR(") + Title + string(")");
    gr->GetYaxis()->SetTitle(title.c_str());
    gr->SetTitle("N production BR vs N mass");
  }
  else {
    string title = string("BR(") + Title + string(")");
    gr->GetYaxis()->SetTitle(title.c_str());
    gr->SetTitle("N decay BR vs N mass");
  }

  gr->GetXaxis()->SetTitle("N mass [GeV]");
  gr->GetXaxis()->SetTitleOffset(1.2);
  gr->GetXaxis()->SetTitleSize(labelSize);
  gr->GetYaxis()->SetTitleSize(labelSize);
  gr->GetXaxis()->SetLabelSize(labelSize);
  gr->GetYaxis()->SetLabelSize(labelSize);
  TGaxis::SetMaxDigits(2);
  gr->Draw("AC");
  gPad->SetLogx();
  gPad->SetLogy();
  gr->SetMinimum(1.E-15);
  gr->SetMaximum(1.5);
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);
  gPad->Update();
  gPad->Modified();
  c1->Write();
}

void GeneralPlotsNew() {

  TFile *f = new TFile("/home/li/Desktop/Plots.root","RECREATE");

  // HNL production via two-body decay
  /*
  MassScan(D, e, 1, 1, "D->Ne");
  MassScan(D, mu, 1, 1, "D->Nmu");
  MassScan(DS, e, 1, 1, "DS->Ne");
  MassScan(DS, mu, 1, 1, "DS->Nmu");
  MassScan(DS, tau, 1, 1, "DS->Ntau");
  MassScan(tau, pi, 1, 1, "DS->taunu; tau->Npi");
  */

  // HNL production via three-body decay

  
  
  // HNL decay via two- and three-body decay
  /*
  MassScan(e, e, 0, 0, "N->eenu");
  MassScan(e, mu, 0, 0, "N->emunu");
  MassScan(pi, e, 0, 1, "N->pie");
  MassScan(mu, mu, 0, 0, "N->mumunu");
  MassScan(pi, mu, 0, 1, "N->pimu");
  MassScan(rho, e, 0, 1, "N->rhoe");
  MassScan(rho, mu, 0, 1, "N->rhomu");
  MassScan(e, tau, 0, 0, "N->etaunu");
  MassScan(pi, tau, 0, 1, "N->pitau");
  MassScan(mu, tau, 0, 0, "N->mutaunu");
  MassScan(rho, tau, 0, 1, "N->rhotau");
  MassScan(tau, tau, 0, 0, "N->tautaunu");
  */
  f->Write();
}
