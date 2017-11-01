Double_t GF = 1.1E-11; // MeV^-2                             
Double_t cos2ThetaC = 0.9471;
Double_t fPi = 130.41; // MeV                                   
Double_t fRho = 1.04E5; // MeV^2
Double_t e = 0.511;
Double_t mu = 105.66;
Double_t tau = 1776.82;
Double_t pi = 139.57;
Double_t pi0 = 134.98;
Double_t rho = 775.4;
Double_t D = 1869.62;
Double_t DS = 1968.28;
Double_t D0 = 1864.84;
Double_t K = 493.68;
Double_t K0 = 497.61;
Double_t De2BR = 9.2E-9;
Double_t Dm2BR = 3.74E-4;
Double_t DSe2BR = 1.4E-7;
Double_t DSm2BR = 5.56E-3;
Double_t DSt2BR = 5.48E-2;
Double_t TauPi2BR = 10.82E-2;
const int InitialMass = 100;
const int Mass = 10000;
const int step = 1;
const int Masses = Mass/step;
TString name = "";
Double_t labelSize = 0.05;
Double_t titleSize = 0.07;
Double_t PiE = 0.;
Double_t PiMu = 0.;
Int_t counterProd = 0;
Int_t counterDecay = 0;
Int_t colors[10] = {434, 600, 882, 880, 632, 634, 797, 402, 419}; // cyan+2, blue, violet+2, violet, red, red+2, orange-3, yellow+2, green+3



Double_t PhaseSpace(Double_t Mass1, Double_t Mass2, Double_t Mass3) {

  Double_t phaseSpace = 0.;
  
  phaseSpace = TMath::Power(Mass1*Mass1-Mass2*Mass2-Mass3*Mass3,2) - 4.*Mass2*Mass2*Mass3*Mass3;
  
  return phaseSpace;
}

Double_t PhaseSpaceFactor(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t phaseSpace) {

  Double_t factor = 0.;

  if(phaseSpace > 0.) {
    factor = (Mass1*Mass1*(Mass2*Mass2 + Mass3*Mass3)-TMath::Power(Mass2*Mass2-Mass3*Mass3,2))*TMath::Power(phaseSpace,0.5)/(Mass3*Mass3*TMath::Power(Mass1*Mass1-Mass3*Mass3,2));
  }
  
  return factor;
}

Double_t GammaLeptonNu3(Double_t Mass1, Double_t Mass2, Double_t Mass3) {

  Double_t r, a, b, f;
  Double_t gamma_l_l_nu = 0.;

  if (Mass1 >= (Mass2+Mass3)) {
    if (Mass2 == Mass3) {
      r  = 4*Mass2*Mass3/(Mass1*Mass1);
      a = TMath::Power(1-r,0.5)*(1./3.-7*r/6.-r*r/24.-r*r*r/16.);
      b = r*r*(1-r*r/16.)*TMath::ATanH(TMath::Power(1-r,0.5));
      f = a+b;
      gamma_l_l_nu = GF*GF*TMath::Power(Mass1,5)*f/(64.*TMath::Power(TMath::Pi(),3));
    }
    else if (Mass2 != Mass3) {
      if (Mass2 > Mass3)
	r = Mass2/Mass1;
      else if (Mass2 < Mass3)
	r = Mass3/Mass1;
      else {
	cout<<"[GammaLeptonNu3] N 3-body decay mode should have equal masses"<<endl;
	exit(1);
      }
      a = (GF*GF*TMath::Power(Mass1, 5))/(192*TMath::Power(TMath::Pi(), 3));
      b = 1 - 8.*r*r + 8.*TMath::Power(r, 6) - TMath::Power(r, 8) -12.*TMath::Power(r, 4)*TMath::Log(r*r);
      gamma_l_l_nu = a*b;
    }
  }
  
  return gamma_l_l_nu;
 }

Double_t Gamma2(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t form) {

  Double_t gamma_2 = 0.;
  Double_t a, b, c, d, e, f;

  if (Mass1 >= (Mass2+Mass3)) {
    if (Mass2 == pi || Mass3 == pi) {
      a = (cos2ThetaC*GF*GF*form*form*Mass1*Mass1*Mass1)/(16.*TMath::Pi());
      b = (1. - TMath::Power(Mass2/Mass1 - Mass3/Mass1, 2.))*(1. - TMath::Power(Mass2/Mass1 + Mass3/Mass1, 2.));
      c = 1. - (Mass3*Mass3)/(Mass1*Mass1);
      d = (1. + (Mass3*Mass3)/(Mass1*Mass1))*(Mass2*Mass2)/(Mass1*Mass1);
      e = c*c - d;
      gamma_2 = a*TMath::Sqrt(b)*e;
    }
    else if (Mass2 == rho || Mass3 == rho) {
      a = (cos2ThetaC*GF*GF*form*form*Mass1*Mass1*Mass1)/(8.*TMath::Pi()*Mass2*Mass2);
      b = (1. - TMath::Power(Mass2/Mass1 - Mass3/Mass1, 2.))*(1. - TMath::Power(Mass2/Mass1 + Mass3/Mass1, 2.));
      c = (1. + (Mass3*Mass3)/(Mass1*Mass1))*(Mass2*Mass2/(Mass1*Mass1));
      d = 2*Mass2*Mass2*Mass2*Mass2/(Mass1*Mass1*Mass1*Mass1);
      e = TMath::Power(1. - (Mass3*Mass3)/(Mass1*Mass1), 2.);
      f = c - d + e;
      gamma_2 = a*TMath::Sqrt(b)*e;
    }
    else {
      cout<<"[Gamma2] Unknown N two-body decay mode"<<endl;
      exit(1);
    }
  }
  
  return gamma_2;
}

Double_t GammaTot(Double_t MN) {

  Double_t gammaTot = 0.;

  if (MN < 2*e) 
    gammaTot = 0.;
  else if (MN >= 2*e && MN < (e+mu)) 
    gammaTot = GammaLeptonNu3(MN, e, e);
  else if (MN >= (e+mu) && MN < (e+pi)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu);
  else if (MN >= (e+pi) && MN < 2*mu) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi, e, fPi);
  else if (MN >= 2*mu && MN < (mu+pi)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu);
  else if (MN >= (mu+pi) && MN < (rho+e)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, mu, pi, fPi);
  else if (MN >= (rho+e) && MN < (rho+mu)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, rho, e, fRho);
  else if (MN >= (mu+rho) && MN < (e+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho);
  else if (MN >= (e+tau) && MN < (pi+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + GammaLeptonNu3(MN, e, tau);
  else if (MN >= (pi+tau) && MN < (mu+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + GammaLeptonNu3(MN, e, tau) + Gamma2(MN, pi, tau, fPi);
  else if (MN >= (mu+tau) && MN < (rho+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + GammaLeptonNu3(MN, e, tau) + Gamma2(MN, pi, tau, fPi) + GammaLeptonNu3(MN, mu, tau);
  else if (MN >= (rho+tau) && MN < 2*tau) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + GammaLeptonNu3(MN, e, tau) + Gamma2(MN, pi, tau, fPi) + GammaLeptonNu3(MN, mu, tau) + Gamma2(MN, rho, tau, fRho);
  else if (MN >= 2*tau) 
    gammaTot = GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + GammaLeptonNu3(MN, e, tau) + Gamma2(MN, pi, tau, fPi) + GammaLeptonNu3(MN, mu, tau) + Gamma2(MN, rho, tau, fRho) + GammaLeptonNu3(MN, tau, tau);

  return gammaTot;
}
/*
Double_t GammaTotUSquared(Double_t MN) {

  Double_t USquared = 1.E-6;
  Double_t UeSquared = USquared/20.8;
  Double_t UmuSquared = 16.*UeSquared;
  Double_t UtauSquared = 3.8*UeSquared;
  Double_t gammaTot = 0.;
  
  if (MN < 2*e) 
    gammaTot = 0.;
  else if (MN >= 2*e && MN < (e+mu))
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared;
  else if (MN >= (e+mu) && MN < (e+pi)) 
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared + GammaLeptonNu3(MN, e, mu)*(UeSquared + UmuSquared);
  else if (MN >= (e+pi) && MN < 2*mu) 
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared + GammaLeptonNu3(MN, e, mu)*(UeSquared + UmuSquared) + Gamma2(MN, pi, e, fPi)*UeSquared;
  else if (MN >= 2*mu && MN < (mu+pi)) 
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared + GammaLeptonNu3(MN, e, mu)*(UeSquared + UmuSquared) + Gamma2(MN, pi, e, fPi)*UeSquared + GammaLeptonNu3(MN, mu, mu)*USquared;
  else if (MN >= (mu+pi) && MN < (rho+e)) 
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared + GammaLeptonNu3(MN, e, mu)*(UeSquared + UmuSquared) + Gamma2(MN, pi, e, fPi)*UeSquared + GammaLeptonNu3(MN, mu, mu)*USquared + Gamma2(MN, pi, mu, fPi)*UmuSquared;
  else if (MN >= (rho+e) && MN < (rho+mu)) 
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared + GammaLeptonNu3(MN, e, mu)*(UeSquared + UmuSquared) + Gamma2(MN, pi, e, fPi)*UeSquared + GammaLeptonNu3(MN, mu, mu)*USquared + Gamma2(MN, pi, mu, fPi)*UmuSquared + Gamma2(MN, rho, e, fRho)*UeSquared;
  else if (MN >= (mu+rho) && MN < (e+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared + GammaLeptonNu3(MN, e, mu)*(UeSquared + UmuSquared) + Gamma2(MN, pi, e, fPi)*UeSquared + GammaLeptonNu3(MN, mu, mu)*USquared + Gamma2(MN, pi, mu, fPi)*UmuSquared + Gamma2(MN, rho, e, fRho)*UeSquared + Gamma2(MN, rho, mu, fRho)*UmuSquared;
  else if (MN >= (e+tau) && MN < (pi+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared + GammaLeptonNu3(MN, e, mu)*(UeSquared + UmuSquared) + Gamma2(MN, pi, e, fPi)*UeSquared + GammaLeptonNu3(MN, mu, mu)*USquared + Gamma2(MN, pi, mu, fPi)*UmuSquared + Gamma2(MN, rho, e, fRho)*UeSquared + Gamma2(MN, rho, mu, fRho)*UmuSquared + GammaLeptonNu3(MN, e, tau)*(UeSquared + UtauSquared);
  else if (MN >= (pi+tau) && MN < (mu+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared + GammaLeptonNu3(MN, e, mu)*(UeSquared + UmuSquared) + Gamma2(MN, pi, e, fPi)*UeSquared + GammaLeptonNu3(MN, mu, mu)*USquared + Gamma2(MN, pi, mu, fPi)*UmuSquared + Gamma2(MN, rho, e, fRho)*UeSquared + Gamma2(MN, rho, mu, fRho)*UmuSquared + GammaLeptonNu3(MN, e, tau)*(UeSquared + UtauSquared) + Gamma2(MN, pi, tau, fPi)*UtauSquared;
  else if (MN >= (mu+tau) && MN < (rho+tau)) 
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared + GammaLeptonNu3(MN, e, mu)*(UeSquared + UmuSquared) + Gamma2(MN, pi, e, fPi)*UeSquared + GammaLeptonNu3(MN, mu, mu)*USquared + Gamma2(MN, pi, mu, fPi)*UmuSquared + Gamma2(MN, rho, e, fRho)*UeSquared + Gamma2(MN, rho, mu, fRho)*UmuSquared + GammaLeptonNu3(MN, e, tau)*(UeSquared + UtauSquared) + Gamma2(MN, pi, tau, fPi)*UtauSquared + GammaLeptonNu3(MN, mu, tau)*(UmuSquared + UtauSquared);
  else if (MN >= (rho+tau) && MN < 2*tau) 
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared + GammaLeptonNu3(MN, e, mu)*(UeSquared + UmuSquared) + Gamma2(MN, pi, e, fPi)*UeSquared + GammaLeptonNu3(MN, mu, mu)*USquared + Gamma2(MN, pi, mu, fPi)*UmuSquared + Gamma2(MN, rho, e, fRho)*UeSquared + Gamma2(MN, rho, mu, fRho)*UmuSquared + GammaLeptonNu3(MN, e, tau)*(UeSquared + UtauSquared) + Gamma2(MN, pi, tau, fPi)*UtauSquared + GammaLeptonNu3(MN, mu, tau)*(UmuSquared + UtauSquared) + Gamma2(MN, rho, tau, fRho)*UtauSquared;
  else if (MN >= 2*tau) 
    gammaTot = GammaLeptonNu3(MN, e, e)*USquared + GammaLeptonNu3(MN, e, mu)*(UeSquared + UmuSquared) + Gamma2(MN, pi, e, fPi)*UeSquared + GammaLeptonNu3(MN, mu, mu)*USquared + Gamma2(MN, pi, mu, fPi)*UmuSquared + Gamma2(MN, rho, e, fRho)*UeSquared + Gamma2(MN, rho, mu, fRho)*UmuSquared + GammaLeptonNu3(MN, e, tau)*(UeSquared + UtauSquared) + Gamma2(MN, pi, tau, fPi)*UtauSquared + GammaLeptonNu3(MN, mu, tau)*(UmuSquared + UtauSquared) + Gamma2(MN, rho, tau, fRho)*UtauSquared; + GammaLeptonNu3(MN, tau, tau)*USquared;

  return gammaTot;
}
*/
Double_t ComputeBR(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t Mass4, Double_t Factor, Bool_t Prod, Bool_t TwoBody) {

  Double_t brt = 0.;
  
  if (Prod == kTRUE) {
    if (TwoBody == kTRUE) {
      if (Factor > 0. && Mass1 >= (Mass2+Mass3)) {
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
	  brt = TauPi2BR*Factor*DSt2BR*PhaseSpaceFactor(DS, 0., tau, PhaseSpace(DS, 0., tau));
	else {
	  cout<<"[ComputeBR] Unknown N 2-body production mode"<<endl;
	  exit(1);
	}
      }
    }
    else if (TwoBody == kFALSE) {
      if (Factor > 0. && Mass1 >= (Mass2+Mass3+Mass4)) {
	if (Mass1 == D && Mass3 == K0) {
	  if (Mass4 == e)
	    brt = 0.;
	  else if (Mass4 == mu)
	    brt = 0.;
	  else {
	    cout<<"[ComputeBR] Unknown N 3-body production mode"<<endl;
	    exit(1);
	  }
	}
	else if (Mass1 == D0 && Mass3 == K) {
	  if (Mass4 == e)
	    brt = 0.;
	  else if (Mass4 == mu)
	    brt = 0.;
	  else {
	    cout<<"[ComputeBR] Unknown N 3-body production mode"<<endl;
	    exit(1);
	  }
	}
      }
    }
  }      
  else if (Prod == kFALSE) {
    if (TwoBody == kTRUE) {
      if (Mass2 == pi || Mass3 == pi) {
	if (Gamma2(Mass1, Mass2, Mass3, fPi) > 0.)
	  brt = Gamma2(Mass1, Mass2, Mass3, fPi)/GammaTot(Mass1);
	else
	  brt = 0.;
      }
      else if (Mass2 == rho || Mass3 == rho) {
	if (Gamma2(Mass1, Mass2, Mass3, fRho) > 0.)
	  brt = Gamma2(Mass1, Mass2, Mass3, fRho)/GammaTot(Mass1);
	else
	  brt = 0.;
      }
      else {
	cout<<"[ComputeBR] Unknown N 2-body decay mode"<<endl;
	exit(1);
      }
    }
    else if (TwoBody == kFALSE) {
      if (GammaLeptonNu3(Mass1, Mass2, Mass3) > 0. && GammaTot(Mass1) > 0.)
	brt = GammaLeptonNu3(Mass1, Mass2, Mass3)/GammaTot(Mass1);
      else {
	GammaLeptonNu3(Mass1, Mass2, Mass3);
	GammaTot(Mass1);
	brt = 0.;
      }
    }
  }

  return brt;
}
/*
Double_t ComputeTotalBR(Double_t Mass1, Double_t MN) {

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
void MassScan(Double_t Mass1, Double_t Mass2, Double_t Mass3, Bool_t Prod, Bool_t TwoBody, std::string Title, TMultiGraph* M, TLegend* L) {

  Double_t PS = 0.;
  Double_t PSF = 0.;
  Double_t BR = 0.;
  Double_t xMN[Masses];
  Double_t yBR[Masses];

  if (Prod == kTRUE) {
    if (TwoBody == kTRUE) {
      for (Int_t MN = InitialMass; MN < Mass; MN += step) {
	PS = PhaseSpace(Mass1, MN, Mass2);
	PSF = PhaseSpaceFactor(Mass1, MN, Mass2, PS);
	BR = ComputeBR(Mass1, MN, Mass2, 0, PSF, Prod, TwoBody);
	xMN[MN] = MN/1000.;
	yBR[MN] = BR;
      }
    }
    else if (TwoBody == kFALSE) {
      for (Int_t MN = InitialMass; MN < Mass; MN += step) {
        BR = ComputeBR(Mass1, MN, Mass2, Mass3, 0., Prod, TwoBody);
        xMN[MN] = MN/1000.;
        yBR[MN] = BR;
      }
    }
  }
  else if (Prod == kFALSE) {
    if (TwoBody == kTRUE) {
      for (Int_t MN = InitialMass; MN < Mass; MN += step) {
	BR = ComputeBR(MN, Mass1, Mass2, 0., 0., Prod, TwoBody);
	xMN[MN] = MN/1000.;
	yBR[MN] = BR;
      }
    }
    else if (TwoBody == kFALSE) {
      for (Int_t MN = InitialMass; MN < Mass; MN += step) {
        BR = ComputeBR(MN, Mass1, Mass2, 0., 0., Prod, TwoBody);
        xMN[MN] = MN/1000.;
        yBR[MN] = BR;
      }
    }
  }

  TGraph* gr = new TGraph(Mass, xMN, yBR);
  gr->SetLineWidth(3);
  if (Prod == kTRUE) {
    gr->SetLineColor(colors[counterProd]);
    counterProd++;
  }
  else if (Prod == kFALSE) {
    gr->SetLineColor(colors[counterDecay]);
    counterDecay++;
  }
  L->AddEntry(gr, Title.c_str(), "l");
  M->Add(gr);
}

void GeneralPlotsNewMultigraph() {

  TFile *f = new TFile("/home/li/Desktop/Plots.root","RECREATE");
  TMultiGraph* Mprod = new TMultiGraph("Mprod", "N production modes vs N mass");
  Mprod->SetName("Mprod");
  TMultiGraph* Mdecay = new TMultiGraph("Mdecay", "N decay modes vs N mass");
  Mdecay->SetName("Mdecay");
  auto Lprod = new TLegend(0.1,0.65,0.4,0.9);
  auto Ldecay = new TLegend(0.1,0.65,0.4,0.9);
  TCanvas *cProd = new TCanvas();
  TCanvas *cDecay = new TCanvas();
  
  // HNL production via two-body decay

  MassScan(D,   e,   0., 1, 1, "D->Ne",               Mprod, Lprod);
  MassScan(D,   mu,  0., 1, 1, "D->Nmu",              Mprod, Lprod);
  MassScan(DS,  e,   0., 1, 1, "DS->Ne",              Mprod, Lprod);
  MassScan(DS,  mu,  0., 1, 1, "DS->Nmu",             Mprod, Lprod);
  MassScan(DS,  tau, 0., 1, 1, "DS->Ntau",            Mprod, Lprod);
  MassScan(tau, pi,  0., 1, 1, "DS->taunu; tau->Npi", Mprod, Lprod);

  // HNL production via three-body decay

  //MassScan(D, K0, e, 1, 0, "D->K0eN");
  
  // HNL decay via two- and three-body decay
    
  MassScan(e,   e,   0., 0, 0, "N->eenu",     Mdecay, Ldecay);
  //MassScan(e,   mu,  0., 0, 0, "N->emunu",    Mdecay, Ldecay);
  //MassScan(pi,  e,   0., 0, 1, "N->pie",      Mdecay, Ldecay);
  //MassScan(mu,  mu,  0., 0, 0, "N->mumunu",   Mdecay, Ldecay);
  //MassScan(pi,  mu,  0., 0, 1, "N->pimu",     Mdecay, Ldecay);
  //MassScan(rho, e,   0., 0, 1, "N->rhoe",     Mdecay, Ldecay);
  //MassScan(rho, mu,  0., 0, 1, "N->rhomu",    Mdecay, Ldecay);
  //MassScan(e,   tau, 0., 0, 0, "N->etaunu",   Mdecay, Ldecay);
  //MassScan(pi,  tau, 0., 0, 1, "N->pitau",    Mdecay, Ldecay);
  //MassScan(mu,  tau, 0., 0, 0, "N->mutaunu",  Mdecay, Ldecay);
  //MassScan(rho, tau, 0., 0, 1, "N->rhotau",   Mdecay, Ldecay);
  //MassScan(tau, tau, 0., 0, 0, "N->tautaunu", Mdecay, Ldecay);
  
  TGaxis::SetMaxDigits(2);
  
  Mprod->Draw("AC");
  Mprod->GetXaxis()->SetTitle("N mass [GeV]");
  Mprod->GetYaxis()->SetTitle("BR*F");
  Mprod->GetXaxis()->SetTitleOffset(1.2);
  Mprod->GetXaxis()->SetTitleSize(labelSize);
  Mprod->GetYaxis()->SetTitleSize(labelSize);
  Mprod->GetXaxis()->SetLabelSize(labelSize);
  Mprod->GetYaxis()->SetLabelSize(labelSize);
  gPad->SetLogx();
  gPad->SetLogy();
  Mprod->SetMinimum(1.E-15);
  Mprod->SetMaximum(1.5);
  cProd->SetLeftMargin(0.2);
  cProd->SetBottomMargin(0.2);
  gPad->Update();
  gPad->Modified();
  Lprod->Draw();
  cProd->Write();
  
  Mdecay->Draw("AC");
  Mdecay->GetXaxis()->SetTitle("N mass [GeV]");
  Mdecay->GetYaxis()->SetTitle("BR*F");
  Mdecay->GetXaxis()->SetTitleOffset(1.2);
  Mdecay->GetXaxis()->SetTitleSize(labelSize);
  Mdecay->GetYaxis()->SetTitleSize(labelSize);
  Mdecay->GetXaxis()->SetLabelSize(labelSize);
  Mdecay->GetYaxis()->SetLabelSize(labelSize);
  gPad->SetLogx();
  gPad->SetLogy();
  Mdecay->SetMinimum(1.E-15);
  Mdecay->SetMaximum(1.5);
  cDecay->SetLeftMargin(0.2);
  cDecay->SetBottomMargin(0.2);
  gPad->Update();
  gPad->Modified();
  Ldecay->Draw();
  cDecay->Write();

  //Mprod->Write();
  //Mdecay->Write();
  f->Write();
}
