// physical constants

Double_t GF = 1.17E-11; // MeV^-2                             
Double_t cos2ThetaC = 0.9471;
Double_t fPi = 130.41; // MeV                                   
Double_t fRho = 1.04E5; // MeV^2
Double_t fD = 222.6;
Double_t fDS = 280.1;

// masses

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
Double_t KStar = 891.76;
Double_t K0Star = 895.55;

// lifetimes

Double_t Dlife = 1.04E-3;
Double_t D0life = 4.101E-4;
Double_t DSlife = 5.E-4;
Double_t taulife = 2.91E-4;

// CKM

Double_t Vcs = 0.9734;
Double_t Vcd = 0.2252;
Double_t Vud = 0.9743;

//form factors, pseudoscalar meson

Double_t fDK0 = 0.745; // f+
Double_t fDpi0 = 0.648;
Double_t fD0K = 0.736;
Double_t fD0pi = 0.637;
Double_t gDK0 = -0.495; // f-
Double_t gDpi0 = -0.435;
Double_t gD0K = gDK0;
Double_t gD0pi = gDpi0;

//form factors, vector meson

Double_t fA0D = 0.398;
Double_t fA1D = 0.47;
Double_t fA2D = -1.24;
Double_t fVD = 0.66;
Double_t fA0D0 = 0.4;
Double_t fA1D0 = 0.47;
Double_t fA2D0 = -1.24;
Double_t fVD0 = 0.66;

// BRs

Double_t De2BR = 9.2E-9;
Double_t Dm2BR = 3.74E-4;
Double_t DSe2BR = 1.4E-7;
Double_t DSm2BR = 5.56E-3;
Double_t DSt2BR = 5.48E-2;
Double_t TauPi2BR = 10.82E-2;

// variables for macro

const int InitialMass = 100;
const int Mass = 10000;
const int step = 1;
const int Masses = Mass/step;
TString name = "";
Double_t labelSize = 0.05;
Double_t titleSize = 0.07;
Int_t counterProd = 0;
Int_t counterDecay = 0;
Int_t colors[14] = {602, 434, 887, 861, 623, 632, 797, 803, 402, 419, 416, 1, 922, 886}; //blue+2, cyan+2, violet+7, azure+1, red-9, red, orange-3, orange+3, yellow+2, green+3, green, black, grey+2, violet+6





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
  else {
    //cout<<"[PhaseSpaceFactor] Warning: phasespace < 0."<<endl;
    factor = 0.;
  }
  
  return factor;
}

Double_t TwoBodyBR(Double_t Mass1, Double_t Mass2, Double_t Mass3) {

  Double_t brt, life, V, f, a, b, c, d;

  if (Mass1 == D) {
    life = Dlife;
    V = Vcd;
    f = fD;
  }
  else if (Mass1 == DS) {
    life = DSlife;
    V = Vcs;
    f = fDS;
  }
  else if (Mass1 == tau) {
    life = taulife;
    f = fPi;
    V = Vud;
  }
  else {
    cout<<"[TwoBodyBR] Unknown mother hadron"<<endl;
    exit(1);
  }

  if (Mass3 != e && Mass3 != mu && Mass3 != tau && Mass3 != pi && Mass3 != rho) {
    cout<<"[TwoBodyBR] Unknown 2-body decay"<<endl;
    exit(1);
  }

  if (Mass1 != tau) {
    a = life*GF*GF*f*f*V*V*Mass1*Mass2*Mass2/(8.*TMath::Pi());
    b = 1. - Mass2*Mass2/(Mass1*Mass1) + 2.*Mass3*Mass3/(Mass1*Mass1);
    c = (1. - Mass3*Mass3/(Mass1*Mass1))*Mass3*Mass3/(Mass2*Mass2);
    d = TMath::Power(1. + Mass2*Mass2/(Mass1*Mass1) - Mass3*Mass3/(Mass1*Mass1), 2.) - 4.*Mass2*Mass2/(Mass1*Mass1);
    brt = a*(b+c)*TMath::Sqrt(d);
  }
  else if (Mass1 == tau && Mass3 == pi) {
    a = life*GF*GF*V*V*f*f*Mass1*Mass1*Mass1/(16.*TMath::Pi());
    b = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.) - (1. + Mass2*Mass2/(Mass1*Mass1))*Mass3*Mass3/(Mass1*Mass1);
    c = 1. - ((Mass3 - Mass2)*(Mass3 - Mass2)/(Mass1*Mass1));
    d = 1. - ((Mass3 + Mass2)*(Mass3 + Mass2)/(Mass1*Mass1));
    brt = a*b*TMath::Sqrt(c*d);
  }
  else if (Mass1 == tau && Mass3 == rho) {
    a = life*fRho*fRho*GF*GF*V*V*Mass1*Mass1*Mass1/(8.*TMath::Pi()*Mass3*Mass3);
    b = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.) + (1. + (Mass2*Mass2 - 2.*Mass3*Mass3)/(Mass1*Mass1))*Mass3*Mass3/(Mass1*Mass1);
    c = 1. - ((Mass3 - Mass2)*(Mass3 - Mass2)/(Mass1*Mass1));
    d = 1. - ((Mass3 + Mass2)*(Mass3 + Mass2)/(Mass1*Mass1));
    brt = a*b*TMath::Sqrt(c*d);
  }
    
  return brt;
}

Double_t ThreeBodyBR(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t Mass4) {

  Double_t br = 0.;

  if (Mass3 == K || Mass3 == K0 || Mass3 == pi || Mass3 == pi0) {
    if (Mass1 == D || Mass1 == D0) { 
      Double_t ENmin = Mass2; // N at rest, K and e back to back
      Double_t ENmax = (Mass1*Mass1+Mass2*Mass2-TMath::Power(Mass4+Mass3, 2.))/(2.*Mass1); // N one way, K and e other way, their momenta summed equal to the N one
      Double_t q2min = TMath::Power(Mass2+Mass4, 2.); // sum of masses of lepton pair
      Double_t q2max = TMath::Power(Mass1-Mass3, 2.); // sum of 4momenta of lepton pair, when K at rest and N and e back to back
      Double_t tau, V, f, g, a, b;
    
      if (Mass1 >= (Mass2+Mass3+Mass4)) {
	if (Mass1 == D) {
	  tau = Dlife;
	  if (Mass3 == K0) {
	    V = Vcs;
	    f = fDK0;
	    g = gDK0;
	  }
	  else if (Mass3 == pi0) {
	    V = Vcd;
	    f = fDpi0;
	    g =	gDpi0;
	  }
	  else {
	    cout<<"[ThreeBodyBR] Unknown daughter hadron"<<endl;
	    exit(1);
	  }
	}
	else if (Mass1 == D0) {
	  tau = D0life;
	  if (Mass3 == K) {
	    V = Vcs;
	    f = fD0K;
	    g = gD0K;
	  }
	  else if (Mass3 == pi) {
	    V = Vcd;
	    f = fD0pi;
	    g = gD0pi;
	  }
	  else {
	    cout<<"[ThreeBodyBR] Unknown daughter hadron"<<endl;
	    exit(1);
	  }
	}
      
	a = tau*V*V*GF*GF/(64.*TMath::Power(TMath::Pi(), 3.)*Mass1*Mass1);
      
	//(g*g*(x*(Mass2*Mass2 + Mass4*Mass4) - TMath::Power(Mass2*Mass2 - Mass4*Mass4, 2.)) + 2.*f*g*(Mass2*Mass2*(2.*Mass1*Mass1 - 2.*Mass3*Mass3 -4.*EN*Mass1 - Mass4*Mass4 + Mass2*Mass2 + x) + Mass4*Mass4*(4.*EN*Mass1 + Mass4*Mass4 - Mass2*Mass2 - x)) + f*f*((4.*EN*Mass1 + Mass4*Mass4 - Mass2*Mass2 - x)*(2.*Mass1*Mass1 - 2.*Mass3*Mass3 - 4.*EN*Mass1 - Mass4*Mass4 + Mass2*Mass2 + x) - (2.*Mass1*Mass1 + 2.*Mass3*Mass3 - x)*(x - Mass2*Mass2 - Mass4*Mass4)));
      
	TF2 func("func", "([5]*[5]*(x*([2]*[2] + [4]*[4]) - TMath::Power([2]*[2] - [4]*[4], 2.)) + 2.*[5]*[0]*([2]*[2]*(2.*[1]*[1] - 2.*[3]*[3] -4.*y*[1] - [4]*[4] + [2]*[2] + x) + [4]*[4]*(4.*y*[1] + [4]*[4] - [2]*[2] - x)) + [0]*[0]*((4.*y*[1] + [4]*[4] - [2]*[2] - x)*(2.*[1]*[1] - 2.*[3]*[3] - 4.*y*[1] - [4]*[4] + [2]*[2] + x) - (2.*[1]*[1] + 2.*[3]*[3] - x)*(x - [2]*[2] - [4]*[4])))");
      
	func.SetParameter(0, f);
	func.SetParameter(1, Mass1);
	func.SetParameter(2, Mass2);
	func.SetParameter(3, Mass3);
	func.SetParameter(4, Mass4);
	func.SetParameter(5, g);
	
	ROOT::Math::WrappedMultiTF1 wf1(func, 2);
	ROOT::Math::AdaptiveIntegratorMultiDim ig;
	ig.SetFunction(wf1);
	ig.SetRelTolerance(0.001);
	double xmin[] = {q2min, ENmin};
	double xmax[] = {q2max, ENmax};
	b = ig.Integral(xmin, xmax);
	br = a*b;
      }
      else {
	//cout<<"[ThreeBodyBR] Warning: mother mass smaller than daughter masses"<<endl;
	br = 0.;
      }
    }
  }
  else if (Mass3 == KStar || Mass3 == K0Star) {
    if (Mass1 == D || Mass1 == D0) { 
      Double_t ENmin = Mass2; // N at rest, K and e back to back
      Double_t ENmax = (Mass1*Mass1+Mass2*Mass2-TMath::Power(Mass4+Mass3, 2.))/(2.*Mass1); // N one way, K and e other way, their momenta summed equal to the N one
      Double_t q2min = TMath::Power(Mass2+Mass4, 2.); // sum of masses of lepton pair
      Double_t q2max = TMath::Power(Mass1-Mass3, 2.); // sum of 4momenta of lepton pair, when K at rest and N and e back to back
      Double_t tau, V, f, f1, f2, f3, f4, omega2, Omega2, a, b;
      
      if (Mass1 >= (Mass2+Mass3+Mass4)) {
	if (Mass1 == D) {
	  tau = Dlife;
	  V = Vcs;
	  f1 = fVD/(Mass1+Mass3);
	  f2 = (Mass1+Mass3)*fA1D;
	  f3 = -fA2D/(Mass1+Mass3);
	  f4 = Mass3*(2.*fA0D-fA1D-fA2D) + Mass1*(fA2D-fA1D); // to be multiplied by 1/x
	}
	else if (Mass1 == D0) {
	  tau = D0life;
	  V = Vcs;
	  f1 = fVD0/(Mass1+Mass3);
	  f2 = (Mass1+Mass3)*fA1D0;
	  f3 = -fA2D0/(Mass1+Mass3);
	  f4 = Mass3*(2.*fA0D0-fA1D0-fA2D0) + Mass1*(fA2D0-fA1D0); // to be multiplied by 1/x
	}
	
	omega2 = Mass1*Mass1 - Mass3*Mass3 + Mass2*Mass2 - Mass4*Mass4; // add - 2.*Mass1*y;
	Omega2 = Mass1*Mass1 - Mass3*Mass3; // add -x
	a = tau*V*V*GF*GF/(32.*TMath::Power(TMath::Pi(), 3.)*Mass1*Mass1);
      
	// (f2*f2/2.)*(x - Mass2*Mass2 - Mass4*Mass4 + (omega2 - 2.*Mass1*y)*(Omega2 - x - (omega2 - 2.*Mass1*y))/(Mass3*Mass3))
	// + ((f3+f4*1./x)*(f3+f4*1./x)/2.)*(Mass2*Mass2 + Mass4*Mass4)*(x - Mass2*Mass2 + Mass4*Mass4)*((Omega2 - x)*(Omega2 - x)/(4.*Mass3*Mass3) - x)
	// + 2.*f3*f3*Mass3*Mass3*((Omega2 - x)*(Omega2 - x)/(4.*Mass3*Mass3) - x)*(Mass2*Mass2 + Mass4*Mass4 - x + (omega2 - 2.*Mass1*y)*(Omega2 - x - (omega2 - 2.*Mass1*y))/(Mass3*Mass3))
	// + 2.*f3*(f3+f4*1./x)*(Mass2*Mass2*(omega2 - 2.*Mass1*y) + (Omega2 - x - (omega2 - 2.*Mass1*y))*Mass4*Mass4)*((Omega2 - x)*(Omega2 - x)/(4.*Mass3*Mass3) - x)
	// + 2.*f1*f2*(x*(2.*(omega2 - 2.*Mass1*y) - Omega2 - x) + (Omega2 - x)*(Mass2*Mass2 - Mass4*Mass4))
	// + (f2*(f3+f4*1./x)/2.)*((omega2 - 2.*Mass1*y)*(Omega2 - x)/(Mass3*Mass3)*(Mass2*Mass2 - Mass4*Mass4) + (Omega2 - x)*(Omega2 - x)*Mass4*Mass4/(Mass3*Mass3) + 2.*TMath::Power(Mass2*Mass2 - Mass4*Mass4, 2.) - 2.*x*(Mass2*Mass2 + Mass4*Mass4))
	// + f2*f3*((Omega2 - x)*(omega2 - 2.*Mass1*y)*(Omega2 - x - (omega2 - 2.*Mass1*y))/(Mass3*Mass3) + 2.*(omega2 - 2.*Mass1*y)*(Mass4*Mass4 - Mass2*Mass2) + (Omega2 - x)*(Mass2*Mass2 - Mass4*Mass4 - x))
	// + f1*f1*((Omega2 - x)*(Omega2 - x)*(x - Mass2*Mass2 + Mass4*Mass4) - 2.*Mass3*Mass3*(x*x - TMath::Power(Mass2*Mass2 - Mass4*Mass4, 2.)) + 2.*(omega2 - 2.*Mass1*y)*(Omega2 - x)*(Mass2*Mass2 - x - Mass4*Mass4) + 2.*(omega2 - 2.*Mass1*y)*(omega2 - 2.*Mass1*y)*x)
	
	TF2 func("func", "(([6]*[6]/2.)*(x - [2]*[2] - [4]*[4] + ([0] - 2.*[9]*y)*([1] - x - ([0] - 2.*[9]*y))/([3]*[3])) + (([7]+[8]*1./x)*([7]+[8]*1./x)/2.)*([2]*[2] + [4]*[4])*(x - [2]*[2] + [4]*[4])*(([1] - x)*([1] - x)/(4.*[3]*[3]) - x) + 2.*[7]*[7]*[3]*[3]*(([1] - x)*([1] - x)/(4.*[3]*[3]) - x)*([2]*[2] + [4]*[4] - x + ([0] - 2.*[9]*y)*([1] - x - ([0] - 2.*[9]*y))/([3]*[3])) + 2.*[7]*([7]+[8]*1./x)*([2]*[2]*([0] - 2.*[9]*y) + ([1] - x - ([0] - 2.*[9]*y))*[4]*[4])*(([1] - x)*([1] - x)/(4.*[3]*[3]) - x) + 2.*[5]*[6]*(x*(2.*([0] - 2.*[9]*y) - [1] + x) + ([1] - x)*([2]*[2] - [4]*[4])) + ([6]*([7]+[8]*1./x)/2.)*(([0] - 2.*[9]*y)*([1] - x)/([3]*[3])*([2]*[2] - [4]*[4]) + ([1] - x)*([1] - x)*[4]*[4]/([3]*[3]) + 2.*TMath::Power([2]*[2] - [4]*[4], 2.) - 2.*x*([2]*[2] + [4]*[4])) + [6]*[7]*(([1] - x)*([0] - 2.*[9]*y)*([1] - x - ([0] - 2.*[9]*y))/([3]*[3]) + 2.*([0] - 2.*[9]*y)*([4]*[4] - [2]*[2]) + ([1] - x)*([2]*[2] - [4]*[4] - x)) + [5]*[5]*(([1] - x)*([1] - x)*(x - [2]*[2] + [4]*[4]) - 2.*[3]*[3]*(x*x - TMath::Power([2]*[2] - [4]*[4], 2.)) + 2.*([0] - 2.*[9]*y)*([1] - x)*([2]*[2] - x - [4]*[4]) + 2.*([0] - 2.*[9]*y)*([0] - 2.*[9]*y)*x))");
      
	func.SetParameter(0, omega2);
	func.SetParameter(1, Omega2);
	func.SetParameter(2, Mass2);
	func.SetParameter(4, Mass4);
	func.SetParameter(3, Mass3);
	func.SetParameter(5, f1);
	func.SetParameter(6, f2);
	func.SetParameter(7, f3);
	func.SetParameter(8, f4);
	func.SetParameter(9, Mass1);
	
	ROOT::Math::WrappedMultiTF1 wf1(func, 2);
	ROOT::Math::AdaptiveIntegratorMultiDim ig;
	ig.SetFunction(wf1);
	ig.SetRelTolerance(0.001);
        double xmin[] = {q2min, ENmin};
        double xmax[] = {q2max, ENmax};
        b = ig.Integral(xmin, xmax);
        br = a*b;
      }
      else {
	//cout<<"[ThreeBodyBR] Warning: mother mass smaller than daughter masses"<<endl;
	br = 0.;
      }
    }
  }
  else if (Mass1 == tau) {
    if (Mass1 >= (Mass2+Mass3+Mass4)) {
      Double_t a, b, c, d, ENmin, ENmax, EN;
      Double_t life = taulife;
	
      if (Mass3 == 0.1) {
	Mass3 = 0.;
	ENmin = Mass2; // N at rest, l and nu back to back
	ENmax = (Mass1*Mass1+Mass2*Mass2-TMath::Power(Mass4+Mass3, 2.))/(2.*Mass1); // N one way, l and nu other way, their momenta summed equal to the N one

	TF1 *func = new TF1("func", "([0]*[3]*[3]*[1]*[1]*x/(2.*TMath::Power(TMath::Pi(), 3.)))*(1. + ([2]*[2] - [4]*[4])/([1]*[1]) - 2.*x/[1])*(1. - [4]*[4]/([1]*[1] + [2]*[2] - 2.*x*[1]))*(TMath::Sqrt(x*x - [2]*[2]))");
      
	func->SetParameter(0, life);
	func->SetParameter(1, Mass1);
	func->SetParameter(2, Mass2);
	func->SetParameter(3, GF);
	func->SetParameter(4, Mass4);
      
	ROOT::Math::WrappedTF1 wf1(*func);
	ROOT::Math::GaussLegendreIntegrator ig;
	ig.SetFunction(wf1);
	b = ig.Integral(ENmin, ENmax);
	br = b;
      }
      else if (Mass3 == 0.01) {
	Mass3 = 0.;
	ENmin = Mass2; // N at rest, l and nu back to back
	ENmax = (Mass1*Mass1+Mass2*Mass2-TMath::Power(Mass4+Mass3, 2.))/(2.*Mass1); // N one way, l and nu other way, their momenta summed equal to the N one
	
	TF1 *func = new TF1("func", "([0]*[3]*[3]*[1]*[1]/(4.*TMath::Power(TMath::Pi(), 3.)))*(TMath::Power((1. - [4]*[4]/([1]*[1] + [2]*[2] - 2.*x*[1])), 2.)*TMath::Sqrt(x*x - [2]*[2]))*(([1] - x)*(1. - ([2]*[2] + [4]*[4])/([1]*[1])) - (1. - [4]*[4]/([1]*[1] + [2]*[2] - 2.*x*[1]))*(([1] - x)*([1] - x)/[1] + (x*x - [2]*[2])/(3.*[1])))");
      
	func->SetParameter(0, life);
	func->SetParameter(1, Mass1);
	func->SetParameter(2, Mass2);
	func->SetParameter(3, GF);
	func->SetParameter(4, Mass4);
      
	ROOT::Math::WrappedTF1 wf1(*func);
	ROOT::Math::GaussLegendreIntegrator ig;
	ig.SetFunction(wf1);
	b = ig.Integral(ENmin, ENmax);
	br = b;
      }
      else {
	cout<<"[ThreeBodyBR] Unknown neutrino type in N 3-body production mode"<<endl;
	exit(1);
      }
    }
    else {
      //cout<<"[ThreeBodyBR] Warning: mother mass smaller than daughter masses"<<endl;               
      br = 0.;
    }
  }
  else {
    cout<<"[ThreeBodyBR] Unknown N 3-body production mode"<<endl;
    exit(1);
  }

  return br;
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
  else {
    //cout<<"[GammaLeptonNu3] Warning: mother mass smaller than daughter masses"<<endl;
    gamma_l_l_nu = 0.;
  }
  
  return gamma_l_l_nu;
 }

Double_t Gamma2(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t form) {

  Double_t gamma_2 = 0.;
  Double_t a, b, c, d, f, g;

  if (Mass1 >= (Mass2+Mass3)) {
    if (Mass2 == pi || Mass3 == pi) {
      a = (cos2ThetaC*GF*GF*form*form*Mass1*Mass1*Mass1)/(16.*TMath::Pi());
      b = (1. - TMath::Power(Mass2/Mass1 - Mass3/Mass1, 2.))*(1. - TMath::Power(Mass2/Mass1 + Mass3/Mass1, 2.));
      c = 1. - (Mass3*Mass3)/(Mass1*Mass1);
      d = (1. + (Mass3*Mass3)/(Mass1*Mass1))*(Mass2*Mass2)/(Mass1*Mass1);
      f = c*c - d;
      gamma_2 = a*TMath::Sqrt(b)*f;
    }
    else if (Mass2 == rho || Mass3 == rho) {
      a = (cos2ThetaC*GF*GF*form*form*Mass1*Mass1*Mass1)/(8.*TMath::Pi()*Mass2*Mass2);
      b = (1. - TMath::Power(Mass2/Mass1 - Mass3/Mass1, 2.))*(1. - TMath::Power(Mass2/Mass1 + Mass3/Mass1, 2.));
      c = (1. + (Mass3*Mass3)/(Mass1*Mass1))*(Mass2*Mass2/(Mass1*Mass1));
      d = 2*Mass2*Mass2*Mass2*Mass2/(Mass1*Mass1*Mass1*Mass1);
      f = TMath::Power(1. - (Mass3*Mass3)/(Mass1*Mass1), 2.);
      g = c - d + f;
      gamma_2 = a*TMath::Sqrt(b)*f;
    }
    else {
      cout<<"[Gamma2] Unknown N two-body decay mode"<<endl;
      exit(1);
    }
  }
  else {
    //cout<<"[Gamma2] Warning: mother mass smaller than daughter masses"<<endl;
    gamma_2 = 0.;
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
/*
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
      }
    }
  }
  
  return brt;
}
*/
Double_t ComputeBR(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t Mass4, Double_t Factor, Bool_t Prod, Bool_t TwoBody) {

  Double_t brt = 0.;
  
  if (Prod == kTRUE) {
    if (TwoBody == kTRUE) {
      if (Mass1 != tau && Mass3 != pi) {
	if (Factor > 0. && Mass1 >= (Mass2+Mass3)) 
	  brt = TwoBodyBR(Mass1, Mass2, Mass3);
	else
	  brt = 0.;
      }
      else if (Mass1 == tau && (Mass3 == pi || Mass3 == rho)) {
	if (PhaseSpaceFactor(DS, 0., tau, PhaseSpace(DS, 0., tau)) > 0. && Factor > 0.)
	  brt = TwoBodyBR(Mass1, Mass2, Mass3);
      }
      else {
	cout << "[ComputeBR] Unknown N 2-body production mode"<<endl;
	exit(1);
      }
    }
    else if (TwoBody == kFALSE) 
      brt = ThreeBodyBR(Mass1, Mass2, Mass3, Mass4);
  }
  else if (Prod == kFALSE) {
    if (TwoBody == kTRUE) {
      if (Mass2 == pi || Mass3 == pi) {
	if (Gamma2(Mass1, Mass2, Mass3, fPi) > 0.)
	  brt = Gamma2(Mass1, Mass2, Mass3, fPi)/GammaTot(Mass1);
	else {
	  //cout<<"[ComputeBR] Warning: two-body decay width or total decay width might be 0 or negative"<<endl;
	  brt = 0.;
	}
      }
      else if (Mass2 == rho || Mass3 == rho) {
	if (Gamma2(Mass1, Mass2, Mass3, fRho) > 0.)
	  brt = Gamma2(Mass1, Mass2, Mass3, fRho)/GammaTot(Mass1);
	else {
	  //cout<<"[ComputeBR] Warning: two-body decay width or total decay width might be 0 or negative"<<endl;
	  brt = 0.;
	}
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
	//cout<<"[ComputeBR] Warning: three-body decay width or total decay width might be 0 or negative"<<endl;
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
void MassScan(Double_t Mass1, Double_t Mass2, Double_t Mass3, Bool_t Prod, Bool_t TwoBody, std::string Title, TMultiGraph* M) {

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
	BR = ComputeBR(Mass1, MN, Mass2, 0., PSF, Prod, TwoBody);
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

  gr->SetNameTitle(Title.c_str(), Title.c_str());
  gr->SetLineWidth(5);

  if (Prod == kTRUE) {
    if (counterProd < 14) {
      gr->SetLineColor(colors[counterProd]);
      counterProd++;
    }
    else {
      gr->SetLineColor(colors[counterProd-14]);
      gr->SetLineStyle(7);
      counterProd++;
    }
  }
  else if (Prod == kFALSE) {
    if (counterDecay < 14) {
      gr->SetLineColor(colors[counterDecay]);
      counterDecay++;
    }
    else {
      gr->SetLineColor(colors[counterDecay-14]);
      gr->SetLineStyle(7);
      counterDecay++;
    }
  }
  
  gr->Draw();
  M->Add(gr);
}

void PhaseSpace(Double_t Mass1, Double_t Mass3, Double_t Mass4, std::string Title) {

  Mass1 = Mass1/1000.;
  Mass3 = Mass3/1000.;
  Mass4 = Mass4/1000.;

  for (Double_t Mass2 = 0.1; Mass2 < (Mass1-Mass3-Mass4); Mass2 += 0.1) {

    Int_t N = 100;
    Double_t xmin = (Mass2+Mass3)*(Mass2+Mass3);
    Double_t xmax = (Mass1-Mass4)*(Mass1-Mass4);
    Double_t ymin = (Mass2+Mass4)*(Mass2+Mass4);
    Double_t ymax = (Mass1-Mass3)*(Mass1-Mass3);
    
    TH2D* h2 = new TH2D(Form("%s_%.1f",Title.c_str(), Mass2), "", N, xmin, xmax, N, ymin, ymax);
    TLorentzVector Mother(0., 0., 0., Mass1);
    Double_t Daughters[3] = {Mass2, Mass3, Mass4} ;  
    TGenPhaseSpace Event;
    
    Event.SetDecay(Mother, 3, Daughters);
    
    for (Int_t n = 0; n < 100000; n++) {
      Double_t W = Event.Generate();
      TLorentzVector *p1 = Event.GetDecay(0);
      TLorentzVector *p2 = Event.GetDecay(1);
      TLorentzVector *p3 = Event.GetDecay(2);
      TLorentzVector p12 = *p1 + *p2;
      TLorentzVector p13 = *p1 + *p3;

      h2->Fill(p12.M2(), p13.M2(), W);
    }

    h2->Draw("colz");
    h2->SetTitle(Form("%s_%.1f",Title.c_str(), Mass2));
    h2->GetXaxis()->SetTitle("Leptonic invariant mass [GeV^2]");
    h2->GetYaxis()->SetTitle("Neutrino-meson invariant mass [GeV^2]");
    h2->GetXaxis()->SetRange(h2->FindFirstBinAbove(0., 1)-2, h2->FindLastBinAbove(0., 1)+2);
    h2->GetYaxis()->SetRange(h2->FindFirstBinAbove(0., 2)-2, h2->FindLastBinAbove(0., 2)+2);
    gPad->SaveAs(Form("/home/li/Desktop/HeavyNeutrino/DalitzPlots/%s_%.1f.png", Title.c_str(), Mass2));
  }
}

void GeneralPlots() {

  TMultiGraph* Mprod = new TMultiGraph("Mprod", "N production modes vs N mass");
  Mprod->SetName("Mprod");
  TMultiGraph* Mdecay = new TMultiGraph("Mdecay", "N decay modes vs N mass");
  Mdecay->SetName("Mdecay");
  TCanvas *c = new TCanvas();

  TGaxis::SetMaxDigits(2);

  // Dalitz plots
  /*    
  PhaseSpace(D,  K0,     e,  "D->K0eN");  // 3-body HNL production: Dalitz plots
  PhaseSpace(D,  pi0,    e,  "D->pi0eN");
  PhaseSpace(D0, K,      e,  "D0->KeN");
  PhaseSpace(D0, pi,     e,  "D0->pieN");
  PhaseSpace(D,  K0,     mu, "D->K0muN");
  PhaseSpace(D,  pi0,    mu, "D->pi0muN");
  PhaseSpace(D0, K,      mu, "D0->KmuN");
  PhaseSpace(D0, pi,     mu, "D0->pimuN");
  PhaseSpace(D,  K0Star, e,  "D->K0*eN");
  PhaseSpace(D0, KStar,  e,  "D0->K*eN");
  PhaseSpace(D,  K0Star, mu, "D->K0*muN");
  PhaseSpace(D0, KStar,  mu, "D0->K*muN");
  */
  // HNL production
  /*  
  MassScan(D,   e,   0., 1, 1, "D->Ne",                Mprod);  // HNL production via two-body decay
  MassScan(D,   mu,  0., 1, 1, "D->Nmu",               Mprod);
  MassScan(DS,  e,   0., 1, 1, "DS->Ne",               Mprod);
  MassScan(DS,  mu,  0., 1, 1, "DS->Nmu",              Mprod);
  MassScan(DS,  tau, 0., 1, 1, "DS->Ntau",             Mprod);
  MassScan(tau, pi,  0., 1, 1, "DS->taunu; tau->Npi",  Mprod);
  MassScan(tau, rho, 0., 1, 1, "DS->taunu; tau->Nrho", Mprod);
  MassScan(D,  K0,     e,  1, 0, "D->K0eN",   Mprod);           // HNL production via three-body decay
  MassScan(D,  pi0,    e,  1, 0, "D->pi0eN",  Mprod);
  MassScan(D0, K,      e,  1, 0, "D0->KeN",   Mprod);
  MassScan(D0, pi,     e,  1, 0, "D0->pieN",  Mprod);
  MassScan(D,  K0,     mu, 1, 0, "D->K0muN",  Mprod);
  MassScan(D,  pi0,    mu, 1, 0, "D->pi0muN", Mprod);
  MassScan(D0, K,      mu, 1, 0, "D0->KmuN",  Mprod);
  MassScan(D0, pi,     mu, 1, 0, "D0->pimuN", Mprod);
  MassScan(D,  K0Star, e,  1, 0, "D->K0*eN",   Mprod);
  MassScan(D0, KStar,  e,  1, 0, "D0->K*eN",   Mprod);
  MassScan(D,  K0Star, mu, 1, 0, "D->K0*muN",  Mprod);
  MassScan(D0, KStar,  mu, 1, 0, "D0->K*muN",  Mprod);
  MassScan(tau, 0.1,   e,  1, 0, "DS->taunu; tau->Nenu_tau",  Mprod);
  MassScan(tau, 0.01,  e,  1, 0, "DS->taunu; tau->Nenu_e",    Mprod);
  MassScan(tau, 0.1,   mu, 1, 0, "DS->taunu; tau->Nmunu_tau", Mprod);
  MassScan(tau, 0.01,  mu, 1, 0, "DS->taunu; tau->Nmunu_mu",  Mprod);

  Mprod->Draw("AC");
  Mprod->GetXaxis()->SetTitle("N mass [GeV]");
  Mprod->GetYaxis()->SetTitle("BR");
  Mprod->GetXaxis()->SetTitleOffset(1.2);
  Mprod->GetXaxis()->SetTitleSize(labelSize);
  Mprod->GetYaxis()->SetTitleSize(labelSize);
  Mprod->GetXaxis()->SetLabelSize(labelSize);
  Mprod->GetYaxis()->SetLabelSize(labelSize);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  Mprod->SetMinimum(1.E-20);
  Mprod->SetMaximum(1.E-12);
  Mprod->GetXaxis()->SetLimits(0.1, 5.);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.2);
  c->SetWindowSize(20000., 12000.);
  gPad->BuildLegend(0.818, 0.223, 0.984, 0.881);
  gPad->Update();
  gPad->Modified();
  gPad->Write();

  TImage *img = TImage::Create();
  img->FromPad(gPad);
  img->WriteImage("/home/li/Desktop/HeavyNeutrino/BRs/NProdGraph.png");

  // HNL decay

  MassScan(e,   e,   0., 0, 0, "N->eenu",     Mdecay);  // HNL decay via two- and three-body decay
  MassScan(e,   mu,  0., 0, 0, "N->emunu",    Mdecay);
  MassScan(pi,  e,   0., 0, 1, "N->pie",      Mdecay);
  MassScan(mu,  mu,  0., 0, 0, "N->mumunu",   Mdecay);
  MassScan(pi,  mu,  0., 0, 1, "N->pimu",     Mdecay);
  MassScan(rho, e,   0., 0, 1, "N->rhoe",     Mdecay);
  MassScan(rho, mu,  0., 0, 1, "N->rhomu",    Mdecay);
  MassScan(e,   tau, 0., 0, 0, "N->etaunu",   Mdecay);
  MassScan(pi,  tau, 0., 0, 1, "N->pitau",    Mdecay);
  MassScan(mu,  tau, 0., 0, 0, "N->mutaunu",  Mdecay);
  MassScan(rho, tau, 0., 0, 1, "N->rhotau",   Mdecay);
  MassScan(tau, tau, 0., 0, 0, "N->tautaunu", Mdecay);
  
  Mdecay->Draw("AC");
  Mdecay->GetXaxis()->SetTitle("N mass [GeV]");
  Mdecay->GetYaxis()->SetTitle("BR");
  Mdecay->GetXaxis()->SetTitleOffset(1.2);
  Mdecay->GetXaxis()->SetTitleSize(labelSize);
  Mdecay->GetYaxis()->SetTitleSize(labelSize);
  Mdecay->GetXaxis()->SetLabelSize(labelSize);
  Mdecay->GetYaxis()->SetLabelSize(labelSize);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  Mdecay->SetMinimum(1.E-10);
  Mdecay->SetMaximum(1.5);
  Mdecay->GetXaxis()->SetLimits(0.1, 10.);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.2);
  c->SetWindowSize(20000., 12000.);
  gPad->BuildLegend(0.24, 0.25, 0.42, 0.56);
  gPad->Update();
  gPad->Modified();
  gPad->Write();

  TImage *img1 = TImage::Create();
  img1->FromPad(gPad);
  img1->WriteImage("/home/li/Desktop/HeavyNeutrino/BRs/NDecayGraph.png");
  */
}
