// couplings according to model I, II, III

TRandom3 *R = new TRandom3();
Double_t USquared = 1.;
Double_t UeSquared = 0.;
Double_t UmuSquared = 0.;
Double_t UtauSquared = 0.;

void SetModel(Int_t model) {
  if (model == 1) {
    UmuSquared = USquared/54.;
    UtauSquared = UmuSquared;
    UeSquared = 52.*UmuSquared;
  }
  else if (model == 2) {
    UeSquared = USquared/20.8;
    UmuSquared = 16.*UeSquared;
    UtauSquared = 3.8*UeSquared;
  }
  else if (model == 3) {
    UmuSquared = USquared/5.361;
    UeSquared = 0.061*UmuSquared;
    UtauSquared = 4.3*UmuSquared;
  }
}

// physical constants

Double_t GF = 1.17E-11; // MeV^-2                             
Double_t fPi = 130.41; // MeV                                   
Double_t fRho = 1.04E5; // MeV^2
Double_t fD = 222.6;
Double_t fDS = 280.1;
Double_t fK = 159.8;
Double_t fEta = 1.2*fPi;
Double_t fEtaprime = -0.45*fPi;
//Double_t sigmacc = 2.3*9.; // Taken from Gaia's note. sigmacc = f*sigma, where sigma is in mubarn at sqrt(s) = 82 GeV (400 GeV proton on p/n (m_p = 1 GeV) and f is fragmentation fraction
Double_t sigmacc = 10.; // Taken from Gaia's note. sigma is in mubarn at sqrt(s) = 27.4 GeV (400 GeV proton on p/n (m_p = 1 GeV)
Double_t sigmapp = 40000.; // sigma is in mubarn at sqrt(s) = 27.4 GeV
Double_t sin2thetaW = 0.2223;

// masses

Double_t e = 0.511;
Double_t mu = 105.66;
Double_t tau = 1776.82;
Double_t pi = 139.57;
Double_t pi0 = 134.98;
Double_t rho = 775.4;
Double_t rho0 = 775.49;
Double_t eta = 547.86;
Double_t etaprime = 957.78;  
Double_t D = 1869.62;
Double_t DS = 1968.28;
Double_t D0 = 1864.84;
Double_t K = 493.68;
Double_t K0 = 497.61;
Double_t KStar = 891.76;
Double_t K0Star = 895.55;

// lifetimes

Double_t Dlife = 1.579*1.E9; // MeV^-1, right?                                                          
Double_t D0life = 6.227*1.E8;
Double_t DSlife = 7.595*1.E8;
Double_t taulife = 4.42*1.E8;

// CKM

Double_t Vcs = 0.9734;
Double_t Vcd = 0.2252;
Double_t Vud = 0.9743;
Double_t Vus = 0.2253;

// form factors, pseudoscalar meson

Double_t fDK0 = 0.745; // f+
Double_t fDpi0 = 0.648;
Double_t fD0K = 0.736;
Double_t fD0pi = 0.637;
Double_t gDK0 = -0.495; // f-
Double_t gDpi0 = -0.435;
Double_t gD0K = gDK0;
Double_t gD0pi = gDpi0;

// form factors, vector meson

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

// fragmentation fractions

Double_t ffD = 0.246;
Double_t ffD0 = 0.565;
Double_t ffDS = 0.08;

// variables for macro

const int InitialMass = 1;
const int Mass = 10000;
const int step = 1;
const int Masses = Mass/step;
TString name = "";
Double_t labelSize = 0.05;
Double_t titleSize = 0.07;
Int_t counterProd = 0;
Int_t counterDecay = 0;
Int_t colors[14] = {602, 434, 887, 861, 623, 632, 797, 803, 402, 419, 416, 1, 922, 886}; //blue+2, cyan+2, violet+7, azure+1, red-9, red, orange-3, orange+3, yellow+2, green+3, green, black, grey+2, violet+6





// phasespace for 2-body N production

Double_t PhaseSpace(Double_t Mass1, Double_t Mass2, Double_t Mass3) {

  Double_t phaseSpace = 0.;
  
  phaseSpace = TMath::Power(Mass1*Mass1-Mass2*Mass2-Mass3*Mass3,2) - 4.*Mass2*Mass2*Mass3*Mass3;
  
  return phaseSpace;
}

// phasespace factor for 2-body N production

Double_t PhaseSpaceFactor(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t phaseSpace) {

  Double_t factor = 0.;

  if(phaseSpace > 0.) {
    factor = (Mass1*Mass1*(Mass2*Mass2 + Mass3*Mass3)-TMath::Power(Mass2*Mass2-Mass3*Mass3,2))*TMath::Power(phaseSpace,0.5)/(Mass3*Mass3*TMath::Power(Mass1*Mass1-Mass3*Mass3,2));
  }
  else {
    factor = 0.;
  }
  
  return factor;
}

// BR for 2-body N production

Double_t TwoBodyBR(Double_t Mass1, Double_t Mass2, Double_t Mass3) {

  Double_t brt, life, U, V, f, a, b, c, d;

  if (Mass3 == e)
    U = UeSquared;
  else if (Mass3 == mu)
    U = UmuSquared;
  else if (Mass3 == tau || Mass3 == pi || Mass3 == rho || Mass3 == K)
    U = UtauSquared;
  
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
    V = Vud;
    if (Mass3 == pi)
      f = fPi;
    else if (Mass3 == K)
      f = fK;
    else if (Mass3 == rho)
      f = fRho;
  }
  else {
    cout<<"[TwoBodyBR] Unknown mother hadron"<<endl;
    exit(1);
  }

  if (Mass3 != e && Mass3 != mu && Mass3 != tau && Mass3 != pi && Mass3 != rho && Mass3 != K) {
    cout<<"[TwoBodyBR] Unknown 2-body decay"<<endl;
    exit(1);
  }

  if (Mass1 != tau) {
    a = U*life*GF*GF*f*f*V*V*Mass1*Mass2*Mass2/(8.*TMath::Pi());
    b = 1. - Mass2*Mass2/(Mass1*Mass1) + 2.*Mass3*Mass3/(Mass1*Mass1);
    c = (1. - Mass3*Mass3/(Mass1*Mass1))*Mass3*Mass3/(Mass2*Mass2);
    d = TMath::Power(1. + Mass2*Mass2/(Mass1*Mass1) - Mass3*Mass3/(Mass1*Mass1), 2.) - 4.*Mass2*Mass2/(Mass1*Mass1);
    brt = a*(b+c)*TMath::Sqrt(d);
  }
  else if (Mass1 == tau && (Mass3 == pi || Mass3 == K)) {
    a = U*life*GF*GF*V*V*f*f*Mass1*Mass1*Mass1/(16.*TMath::Pi());
    b = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.) - (1. + Mass2*Mass2/(Mass1*Mass1))*Mass3*Mass3/(Mass1*Mass1);
    c = 1. - ((Mass3 - Mass2)*(Mass3 - Mass2)/(Mass1*Mass1));
    d = 1. - ((Mass3 + Mass2)*(Mass3 + Mass2)/(Mass1*Mass1));
    brt = a*b*TMath::Sqrt(c*d);
  }
  else if (Mass1 == tau && Mass3 == rho) {
    a = U*life*fRho*fRho*GF*GF*V*V*Mass1*Mass1*Mass1/(8.*TMath::Pi()*Mass3*Mass3);
    b = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.) + (1. + (Mass2*Mass2 - 2.*Mass3*Mass3)/(Mass1*Mass1))*Mass3*Mass3/(Mass1*Mass1);
    c = 1. - ((Mass3 - Mass2)*(Mass3 - Mass2)/(Mass1*Mass1));
    d = 1. - ((Mass3 + Mass2)*(Mass3 + Mass2)/(Mass1*Mass1));
    brt = a*b*TMath::Sqrt(c*d);
  }
    
  return brt;
}

// BR for 3-body N production

Double_t ThreeBodyBR(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t Mass4) {

  Double_t br = 0.;
  Double_t U = 0.;
  
  if (Mass4 == e)
    U = UeSquared;
  else if (Mass4 == mu)
    U = UmuSquared;
  else if (Mass4 == tau)
    U = UtauSquared;

  if (Mass3 == K || Mass3 == K0 || Mass3 == pi || Mass3 == pi0) {
    if (Mass1 == D || Mass1 == D0) { 
      Double_t mKl = -1.;
      Double_t mNl = -1.;
      Double_t minKl = TMath::Power(Mass3+Mass4, 2.);
      Double_t maxKl = TMath::Power(Mass1-Mass2, 2.);
      Double_t minNl = TMath::Power(Mass2+Mass4, 2.);
      Double_t maxNl = TMath::Power(Mass1-Mass3, 2.);
      Double_t EN = 0.;
      Double_t El = 0.;
      Double_t PN = 0.;
      Double_t Pl = 0.;

      while(mNl < minNl || mNl > maxNl){
  
	mNl = R->Rndm();
	mNl = minNl + mNl*(maxNl-minNl);
	mKl = R->Rndm();
	mKl = minKl + mKl*(maxKl-minKl);
	EN = (Mass1*Mass1 - mKl + Mass2*Mass2)/(2.*Mass1);
	El = (mKl - Mass3*Mass3 + Mass4*Mass4)/(2.*TMath::Sqrt(mKl));
	PN = TMath::Sqrt(EN*EN - Mass2*Mass2); 
	Pl = TMath::Sqrt(El*El - Mass4*Mass4);
	maxNl = (EN+El)*(EN+El) - (PN-Pl)*(PN-Pl); 
	minNl = (EN+El)*(EN+El) - (PN+Pl)*(PN+Pl); 
      }

      Double_t q2min = TMath::Power(Mass2+Mass4, 2.);
      Double_t q2max = TMath::Power(Mass1-Mass3, 2.);
      Double_t ENmin = (Mass1*Mass1 - maxKl + Mass2*Mass2)/(2.*Mass1);
      Double_t ENmax = (Mass1*Mass1 - minKl + Mass2*Mass2)/(2.*Mass1);
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
      
	a = U*tau*V*V*GF*GF/(64.*TMath::Power(TMath::Pi(), 3.)*Mass1*Mass1);
      
	//(g*g*(x*(Mass2*Mass2 + Mass4*Mass4) - TMath::Power(Mass2*Mass2 - Mass4*Mass4, 2.)) + 2.*f*g*(Mass2*Mass2*(2.*Mass1*Mass1 - 2.*Mass3*Mass3 -4.*EN*Mass1 - Mass4*Mass4 + Mass2*Mass2 + x) + Mass4*Mass4*(4.*EN*Mass1 + Mass4*Mass4 - Mass2*Mass2 - x)) + f*f*((4.*EN*Mass1 + Mass4*Mass4 - Mass2*Mass2 - x)*(2.*Mass1*Mass1 - 2.*Mass3*Mass3 - 4.*EN*Mass1 - Mass4*Mass4 + Mass2*Mass2 + x) - (2.*Mass1*Mass1 + 2.*Mass3*Mass3 - x)*(x - Mass2*Mass2 - Mass4*Mass4)));
      
	TF2 func("func", "([5]*[5]*(x*([2]*[2] + [4]*[4]) - TMath::Power([2]*[2] - [4]*[4], 2.)) + 2.*[5]*[0]*([2]*[2]*(2.*[1]*[1] - 2.*[3]*[3] -4.*y*[1] - [4]*[4] + [2]*[2] + x) + [4]*[4]*(4.*y*[1] + [4]*[4] - [2]*[2] - x)) + [0]*[0]*((4.*y*[1] + [4]*[4] - [2]*[2] - x)*(2.*[1]*[1] - 2.*[3]*[3] - 4.*y*[1] - [4]*[4] + [2]*[2] + x) + (2.*[1]*[1] + 2.*[3]*[3] - x)*(x - [2]*[2] - [4]*[4])))");
      
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
	br = 0.;
      }
    }
  }
  else if (Mass3 == KStar || Mass3 == K0Star) {
    if (Mass1 == D || Mass1 == D0) { 
      Double_t mKl = -1.;
      Double_t mNl = -1.;
      Double_t minKl = TMath::Power(Mass3+Mass4, 2.);
      Double_t maxKl = TMath::Power(Mass1-Mass2, 2.);
      Double_t minNl = TMath::Power(Mass2+Mass4, 2.);
      Double_t maxNl = TMath::Power(Mass1-Mass3, 2.);
      Double_t EN = 0.;
      Double_t El = 0.;
      Double_t PN = 0.;
      Double_t Pl = 0.;

      while(mNl < minNl || mNl > maxNl){
  
	mNl = R->Rndm();
	mNl = minNl + mNl*(maxNl-minNl);
	mKl = R->Rndm();
	mKl = minKl + mKl*(maxKl-minKl);
	EN = (Mass1*Mass1 - mKl + Mass2*Mass2)/(2.*Mass1);
	El = (mKl - Mass3*Mass3 + Mass4*Mass4)/(2.*TMath::Sqrt(mKl));
	PN = TMath::Sqrt(EN*EN - Mass2*Mass2); 
	Pl = TMath::Sqrt(El*El - Mass4*Mass4);
	maxNl = (EN+El)*(EN+El) - (PN-Pl)*(PN-Pl); 
	minNl = (EN+El)*(EN+El) - (PN+Pl)*(PN+Pl); 
      }

      Double_t q2min = TMath::Power(Mass2+Mass4, 2.);
      Double_t q2max = TMath::Power(Mass1-Mass3, 2.);
      Double_t ENmin = (Mass1*Mass1 - maxKl + Mass2*Mass2)/(2.*Mass1);
      Double_t ENmax = (Mass1*Mass1 - minKl + Mass2*Mass2)/(2.*Mass1);
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
	a = U*tau*V*V*GF*GF/(32.*TMath::Power(TMath::Pi(), 3.)*Mass1*Mass1);
      
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
	br = 0.;
      }
    }
  }
  else if (Mass1 == tau) {
    if (Mass1 >= (Mass2+Mass3+Mass4)) {
      Double_t a, b, c, d, ENmin, ENmax, EN;
      Double_t life = taulife;
	
      if (Mass3 == 0.1) { //nu_tau
	if (Mass4 == e)
	  U = UeSquared;
	else if (Mass4 == mu)
	  U = UmuSquared;
	else if (Mass4 == tau)
	  U = UtauSquared;

	Mass3 = 0.;
	ENmin = Mass2; // N at rest, l and nu back to back
	ENmax = (Mass1*Mass1+Mass2*Mass2-TMath::Power(Mass4+Mass3, 2.))/(2.*Mass1); // N one way, l and nu other way, their momenta summed equal to the N one

	TF1 *func = new TF1("func", "([5]*[0]*[3]*[3]*[1]*[1]*x/(2.*TMath::Power(TMath::Pi(), 3.)))*(1. + ([2]*[2] - [4]*[4])/([1]*[1]) - 2.*x/[1])*(1. - [4]*[4]/([1]*[1] + [2]*[2] - 2.*x*[1]))*(TMath::Sqrt(x*x - [2]*[2]))");
      
	func->SetParameter(0, life);
	func->SetParameter(1, Mass1);
	func->SetParameter(2, Mass2);
	func->SetParameter(3, GF);
	func->SetParameter(4, Mass4);
	func->SetParameter(5, U);
      
	ROOT::Math::WrappedTF1 wf1(*func);
	ROOT::Math::GaussLegendreIntegrator ig;
	ig.SetFunction(wf1);
	b = ig.Integral(ENmin, ENmax);
	br = b;
      }
      else if (Mass3 == 0.01) { //nu_e or nu_mu
	U = UtauSquared;
	Mass3 = 0.;
	ENmin = Mass2; // N at rest, l and nu back to back
	ENmax = (Mass1*Mass1+Mass2*Mass2-TMath::Power(Mass4+Mass3, 2.))/(2.*Mass1); // N one way, l and nu other way, their momenta summed equal to the N one
	
	TF1 *func = new TF1("func", "([5]*[0]*[3]*[3]*[1]*[1]/(4.*TMath::Power(TMath::Pi(), 3.)))*(TMath::Power((1. - [4]*[4]/([1]*[1] + [2]*[2] - 2.*x*[1])), 2.)*TMath::Sqrt(x*x - [2]*[2]))*(([1] - x)*(1. - ([2]*[2] + [4]*[4])/([1]*[1])) - (1. - [4]*[4]/([1]*[1] + [2]*[2] - 2.*x*[1]))*(([1] - x)*([1] - x)/[1] + (x*x - [2]*[2])/(3.*[1])))");
      
	func->SetParameter(0, life);
	func->SetParameter(1, Mass1);
	func->SetParameter(2, Mass2);
	func->SetParameter(3, GF);
	func->SetParameter(4, Mass4);
      	func->SetParameter(5, U);
	
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
      br = 0.;
    }
  }
  else {
    cout<<"[ThreeBodyBR] Unknown N 3-body production mode"<<endl;
    exit(1);
  }

  return br;
}

// Decay width for 2-body N decay

Double_t Gamma2(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t form) {

  Double_t gamma_2 = 0.;
  Double_t U, V, a, b, c, d, f, g;

  if (Mass1 >= (Mass2+Mass3)) {
    if (Mass2 == pi || Mass2 == K) {
      if (Mass2 == pi)
	V = Vud;
      else if (Mass2 == K)
	V = Vus;
      
      if (Mass3 == e)
	U = UeSquared;
      else if (Mass3 == mu)
	U = UmuSquared;
      else if (Mass3 == tau)
	U = UtauSquared;
      
      a = (U*GF*GF*V*V*form*form*Mass1*Mass1*Mass1)/(16.*TMath::Pi());
      b = TMath::Power(1. - Mass3*Mass3/(Mass1*Mass1), 2.);
      c = (1. + Mass3*Mass3/(Mass1*Mass1))*(Mass2*Mass2/(Mass1*Mass1));
      d = 1. - (Mass2 - Mass3)*(Mass2 - Mass3)/(Mass1*Mass1);
      f = 1. - (Mass2 + Mass3)*(Mass2 + Mass3)/(Mass1*Mass1);
      gamma_2 = a*(b - c)*TMath::Sqrt(d*f);
    }
    else if (Mass2 == rho) {
      V = Vud;
      
      if (Mass3 == e)
	U = UeSquared;
      else if (Mass3 == mu)
	U = UmuSquared;
      else if (Mass3 == tau)
	U = UtauSquared;
      
      a = (U*GF*GF*V*V*form*form*Mass1*Mass1*Mass1)/(8.*TMath::Pi()*Mass2*Mass2);
      b = (1. - TMath::Power(Mass2/Mass1 - Mass3/Mass1, 2.))*(1. - TMath::Power(Mass2/Mass1 + Mass3/Mass1, 2.));
      c = (1. + (Mass3*Mass3)/(Mass1*Mass1))*(Mass2*Mass2/(Mass1*Mass1));
      d = 2*Mass2*Mass2*Mass2*Mass2/(Mass1*Mass1*Mass1*Mass1);
      f = TMath::Power(1. - (Mass3*Mass3)/(Mass1*Mass1), 2.);
      g = c - d + f;
      gamma_2 = a*TMath::Sqrt(b)*g;
    }
    else if (Mass2 == pi0) {
      U = USquared;
      a = (U*GF*GF*form*form*Mass1*Mass1*Mass1)/(32.*TMath::Pi());
      b = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.);
      gamma_2 = a*b;
    }
    else if (Mass2 == eta || Mass2 == etaprime) {
      U = USquared;
      a = (U*GF*GF*form*form*Mass1*Mass1*Mass1)/(32.*TMath::Pi());
      b = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.);
      gamma_2 = a*b;
    }
    else if (Mass2 == rho0) {
      U= USquared;
      a = (U*form*form*GF*GF*Mass1*Mass1*Mass1)/(16.*TMath::Pi()*Mass2*Mass2);
      b = 1. + 2.*Mass2*Mass2/(Mass1*Mass1);
      c = TMath::Power(1. - Mass2*Mass2/(Mass1*Mass1), 2.);
      gamma_2 = a*b*c;
    }
    else {
      cout<<"[Gamma2] Unknown N two-body decay mode"<<endl;
      exit(1);
    }
  }
  else {
    gamma_2 = 0.;
  }
  
  return 2.*gamma_2;
}

// Decay width for 3-body N decay

Double_t GammaLeptonNu3(Double_t Mass1, Double_t Mass2, Double_t Mass3) {

  Double_t U, U1, r, a, b, b1, c, c1, d, d1, f, f1, g, g1, j, j1, L, C1, C2, C3, C4;
  Double_t gamma_l_l_nu = 0.;
  Double_t gamma_l_l_nu0 = 0.;
  Double_t gamma_l_l_nu1 = 0.;

  if (Mass1 >= (Mass2+Mass3)) {
    if (Mass2 == Mass3 && Mass2 != 0.) {
      r = Mass2/Mass1;
      L = TMath::Log((1. - 3.*r*r - (1. - r*r)*(TMath::Sqrt(1. - 4.*r*r)))/(r*r*(1. + TMath::Sqrt(1. - 4.*r*r))));
      if (Mass2 == e) {
	U = UeSquared;
	U1 = UmuSquared + UtauSquared;
	L = -100.;
      }
      else if (Mass2 == mu) {
	U = UmuSquared;
	U1 = UeSquared + UtauSquared;
      }
      else if (Mass2 == tau) {
	U = UtauSquared;
	U1 = UeSquared + UmuSquared;
      }
      a = GF*GF*TMath::Power(Mass1, 5.)/(192.*TMath::Power(TMath::Pi(), 3.));
      C1 = 0.25*(1. - 4.*sin2thetaW + 8.*sin2thetaW*sin2thetaW);
      C2 = 0.5*sin2thetaW*(2.*sin2thetaW - 1.);
      C3 = 0.25*(1. + 4.*sin2thetaW + 8.*sin2thetaW*sin2thetaW);
      C4 = 0.5*sin2thetaW*(2.*sin2thetaW + 1.);
      b = C3;
      c = (1. - 14.*r*r - 2.*TMath::Power(r, 4.) - 12.*TMath::Power(r, 6.))*(TMath::Sqrt(1. - 4.*r*r));
      d = 12.*TMath::Power(r, 4.)*(TMath::Power(r, 4.) - 1.)*L;
      f = 4.*C4;
      g = r*r*(2. + 10.*r*r - 12.*TMath::Power(r, 4.))*(TMath::Sqrt(1. - 4.*r*r));
      j = 6.*TMath::Power(r, 4.)*(1. - 2.*r*r + 2.*TMath::Power(r, 4.))*L;
      gamma_l_l_nu0 = U*a*(b*(c+d) + f*(g+j));
      b1 = C1;
      c1 = c;
      d1 = d;
      f1 = 4.*C2;
      g1 = g;
      j1 = j;
      gamma_l_l_nu1 = U1*a*(b1*(c1+d1) + f1*(g1+j1));
      gamma_l_l_nu = gamma_l_l_nu0 + gamma_l_l_nu1;
    }
    else if (Mass2 == Mass3 && Mass2 == 0.) {
      U = USquared;
      a = U*GF*GF*TMath::Power(Mass1, 5.)/(192.*TMath::Power(TMath::Pi(), 3.));
      gamma_l_l_nu = a;
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
      if ((Mass2 == e || Mass3 == e) && (Mass2 == mu || Mass3 == mu))
	U = UeSquared + UmuSquared;
      else if ((Mass2 == e || Mass3 == e) && (Mass2 == tau || Mass3 == tau))
	U = UeSquared + UtauSquared;
      else if ((Mass2 == mu || Mass3 == mu) && (Mass2 == tau || Mass3 == tau))
	U = UmuSquared + UtauSquared;
      
      a = (U*GF*GF*TMath::Power(Mass1, 5))/(192*TMath::Power(TMath::Pi(), 3));
      b = 1 - 8.*r*r + 8.*TMath::Power(r, 6) - TMath::Power(r, 8) -12.*TMath::Power(r, 4)*TMath::Log(r*r);
      gamma_l_l_nu = a*b;
    }
  }
  else {
    gamma_l_l_nu = 0.;
  }
  
  return 2.*gamma_l_l_nu;
 }

// Total N decay width

Double_t GammaTot(Double_t MN) {

  Double_t gammaTot = 0.;
  
  if (MN < 2*e) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.);
  else if (MN >= 2*e && MN < (e+mu)) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e);
  else if (MN >= (e+mu) && MN < (pi0)) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu);
  else if (MN >= (pi0) && MN < (e+pi))
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi);
  else if (MN >= (e+pi) && MN < 2*mu) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi);
  else if (MN >= 2*mu && MN < (mu+pi)) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu);
  else if (MN >= (mu+pi) && MN < (K+e)) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi);
  else if (MN >= (K+e) && MN < (eta))
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK);
  else if (MN >= (eta) && MN < (K+mu))
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta);
  else if (MN >= (K+mu) && MN < (rho0))
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta) + Gamma2(MN, K, mu, fK);
  else if (MN >= (rho0) && MN < (rho+e))
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta) + Gamma2(MN, K, mu, fK) + Gamma2(MN, rho0, 0., fRho);
  else if (MN >= (rho+e) && MN < (rho+mu)) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta) + Gamma2(MN, K, mu, fK) + Gamma2(MN, rho0, 0., fRho) + Gamma2(MN, rho, e, fRho);
  else if (MN >= (mu+rho) && MN < (etaprime))
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta) + Gamma2(MN, K, mu, fK) + Gamma2(MN, rho0, 0., fRho) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho);
  else if (MN >= (etaprime) && MN < (e+tau)) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta) + Gamma2(MN, K, mu, fK) + Gamma2(MN, rho0, 0., fRho) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + Gamma2(MN, etaprime, 0., fEtaprime);
  
  else if (MN >= (e+tau) && MN < (mu+tau)) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta) + Gamma2(MN, K, mu, fK) + Gamma2(MN, rho0, 0., fRho) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + Gamma2(MN, etaprime, 0., fEtaprime) + GammaLeptonNu3(MN, e, tau);
  else if (MN >= (mu+tau) && MN < (pi+tau)) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta) + Gamma2(MN, K, mu, fK) + Gamma2(MN, rho0, 0., fRho) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + Gamma2(MN, etaprime, 0., fEtaprime) + GammaLeptonNu3(MN, e, tau) + GammaLeptonNu3(MN, mu, tau);
  else if (MN >= (pi+tau) && MN < (K+tau)) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta) + Gamma2(MN, K, mu, fK) + Gamma2(MN, rho0, 0., fRho) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + Gamma2(MN, etaprime, 0., fEtaprime) + GammaLeptonNu3(MN, e, tau) + GammaLeptonNu3(MN, mu, tau) + Gamma2(MN, pi, tau, fPi);
  else if (MN >= (K+tau) && MN < (rho+tau)) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta) + Gamma2(MN, K, mu, fK) + Gamma2(MN, rho0, 0., fRho) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + Gamma2(MN, etaprime, 0., fEtaprime) + GammaLeptonNu3(MN, e, tau) + GammaLeptonNu3(MN, mu, tau) + Gamma2(MN, pi, tau, fPi) + Gamma2(MN, K, tau, fK);
  else if (MN >= (rho+tau) && MN < 2*tau)
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta) + Gamma2(MN, K, mu, fK) + Gamma2(MN, rho0, 0., fRho) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + Gamma2(MN, etaprime, 0., fEtaprime) + GammaLeptonNu3(MN, e, tau) + GammaLeptonNu3(MN, mu, tau) + Gamma2(MN, pi, tau, fPi) + Gamma2(MN, K, tau, fK) + Gamma2(MN, rho, tau, fRho);
  else if (MN >= 2*tau) 
    gammaTot = GammaLeptonNu3(MN, 0., 0.) + GammaLeptonNu3(MN, e, e) + GammaLeptonNu3(MN, e, mu) + Gamma2(MN, pi0, 0., fPi) + Gamma2(MN, pi, e, fPi) + GammaLeptonNu3(MN, mu, mu) + Gamma2(MN, pi, mu, fPi) + Gamma2(MN, K, e, fK) + Gamma2(MN, eta, 0., fEta) + Gamma2(MN, K, mu, fK) + Gamma2(MN, rho0, 0., fRho) + Gamma2(MN, rho, e, fRho) + Gamma2(MN, rho, mu, fRho) + Gamma2(MN, etaprime, 0., fEtaprime) + GammaLeptonNu3(MN, e, tau) + GammaLeptonNu3(MN, mu, tau) + Gamma2(MN, pi, tau, fPi) + Gamma2(MN, rho, tau, fRho) + GammaLeptonNu3(MN, tau, tau);

  return 2.*gammaTot;
}

// Function to distinguish between production/decay and 2-body/3-body

Double_t ComputeBR(Double_t Mass1, Double_t Mass2, Double_t Mass3, Double_t Mass4, Double_t Factor, Bool_t Prod, Bool_t TwoBody) {

  Double_t brt = 0.;
  Double_t gammat = GammaTot(Mass1);
  
  if (Prod == kTRUE) {
    if (TwoBody == kTRUE) {
      if (Mass1 != tau && Mass3 != pi) {
	if (Factor > 0. && Mass1 >= (Mass2+Mass3)) 
	  brt = TwoBodyBR(Mass1, Mass2, Mass3);
	else
	  brt = 0.;
      }
      else if (Mass1 == tau && (Mass3 == pi || Mass3 == rho || Mass3 == K)) {
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
      if (Mass2 == pi) {
	if (Gamma2(Mass1, Mass2, Mass3, fPi) > 0.)
	  brt = Gamma2(Mass1, Mass2, Mass3, fPi)/gammat;
	else {
	  brt = 0.;
	}
      }
      else if (Mass2 == rho) {
	if (Gamma2(Mass1, Mass2, Mass3, fRho) > 0.)
	  brt = Gamma2(Mass1, Mass2, Mass3, fRho)/gammat;
	else {
	  brt = 0.;
	}
      }
      else if (Mass2 == pi0) {
        if (Gamma2(Mass1, Mass2, Mass3, fPi) > 0.)
          brt = Gamma2(Mass1, Mass2, Mass3, fPi)/gammat;
        else {
          brt = 0.;
        }
      }
      else if (Mass2 == K) {
        if (Gamma2(Mass1, Mass2, Mass3, fK) > 0.)
          brt = Gamma2(Mass1, Mass2, Mass3, fK)/gammat;
        else {
          brt = 0.;
        }
      }
      else if (Mass2 == eta) {
        if (Gamma2(Mass1, Mass2, Mass3, fEta) > 0.)
          brt = Gamma2(Mass1, Mass2, Mass3, fEta)/gammat;
        else {
          brt = 0.;
        }
      }
      else if (Mass2 == etaprime) {
        if (Gamma2(Mass1, Mass2, Mass3, fEtaprime) > 0.)
          brt = Gamma2(Mass1, Mass2, Mass3, fEtaprime)/gammat;
        else {
          brt = 0.;
        }
      }
      else if (Mass2 == rho0) {
	if (Gamma2(Mass1, Mass2, Mass3, fRho) > 0.)
	  brt = Gamma2(Mass1, Mass2, Mass3, fRho)/gammat;
	else {
	  brt = 0.;
	}
      }
      else {
	cout<<"[ComputeBR] Unknown N 2-body decay mode"<<endl;
	exit(1);
      }
    }
    else if (TwoBody == kFALSE) {
      if (GammaLeptonNu3(Mass1, Mass2, Mass3) > 0. && gammat > 0.)
	brt = GammaLeptonNu3(Mass1, Mass2, Mass3)/gammat;
      else {
	brt = 0.;
      }
    }
  }

  return brt;
}

// Function for partial decay widths

Double_t DecayWidth(Double_t Mass1, Double_t Mass2, Double_t Mass3, Bool_t TwoBody) {

  Double_t gamma = 0.;
  
  if (TwoBody == kTRUE) {
    if (Mass2 == pi) {
      if (Gamma2(Mass1, Mass2, Mass3, fPi) > 0.)
	gamma = Gamma2(Mass1, Mass2, Mass3, fPi);
      else {
	gamma = 0.;
      }
    }
    else if (Mass2 == rho) {
      if (Gamma2(Mass1, Mass2, Mass3, fRho) > 0.)
	gamma = Gamma2(Mass1, Mass2, Mass3, fRho);
      else {
	gamma = 0.;
      }
    }
    else if (Mass2 == pi0) {
      if (Gamma2(Mass1, Mass2, Mass3, fPi) > 0.)
	gamma = Gamma2(Mass1, Mass2, Mass3, fPi);
      else {
	gamma = 0.;
      }
    }
    else if (Mass2 == K) {
      if (Gamma2(Mass1, Mass2, Mass3, fK) > 0.)
	gamma = Gamma2(Mass1, Mass2, Mass3, fK);
      else {
	gamma = 0.;
      }
    }
    else if (Mass2 == eta) {
      if (Gamma2(Mass1, Mass2, Mass3, fEta) > 0.)
	gamma = Gamma2(Mass1, Mass2, Mass3, fEta);
      else {
	gamma = 0.;
      }
    }
    else if (Mass2 == etaprime) {
      if (Gamma2(Mass1, Mass2, Mass3, fEtaprime) > 0.)
	gamma = Gamma2(Mass1, Mass2, Mass3, fEtaprime);
      else {
	gamma = 0.;
      }
    }
    else if (Mass2 == rho0) {
      if (Gamma2(Mass1, Mass2, Mass3, fRho) > 0.)
	gamma = Gamma2(Mass1, Mass2, Mass3, fRho);
      else {
	gamma = 0.;
      }
    }
    else {
      cout<<"[ComputeBR] Unknown N 2-body decay mode"<<endl;
      exit(1);
    }
  }
  else if (TwoBody == kFALSE) {
    if (GammaLeptonNu3(Mass1, Mass2, Mass3) > 0.)
      gamma = GammaLeptonNu3(Mass1, Mass2, Mass3);
    else {
      gamma = 0.;
    }
  }

  return gamma;
}

// Function for scan on N mass

void MassScan(Double_t Mass1, Double_t Mass2, Double_t Mass3, Bool_t Prod, Bool_t TwoBody, Bool_t Gamma, std::string Title, TMultiGraph* M, Double_t factor) {

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
	BR = factor*ComputeBR(Mass1, MN, Mass2, 0., PSF, Prod, TwoBody);
	xMN[MN] = MN/1000.;
	yBR[MN] = BR;
      }
    }
    else if (TwoBody == kFALSE) {
      for (Int_t MN = InitialMass; MN < Mass; MN += step) {
        BR = factor*ComputeBR(Mass1, MN, Mass2, Mass3, 0., Prod, TwoBody);
        xMN[MN] = MN/1000.;
        yBR[MN] = BR;
      }
    }
  }
  else if (Prod == kFALSE) {
    if (Gamma == kTRUE) {
      for (Int_t MN = InitialMass; MN < Mass; MN += step) {
	BR = DecayWidth(MN, Mass1, Mass2, TwoBody);
	xMN[MN] = MN/1000.;
	yBR[MN] = BR;
      }
    }
    else {
      if (TwoBody == kTRUE) {
	for (Int_t MN = InitialMass; MN < Mass; MN += step) {
	  BR = factor*ComputeBR(MN, Mass1, Mass2, 0., 0., Prod, TwoBody);
	  xMN[MN] = MN/1000.;
	  yBR[MN] = BR;
	}
      }
      else if (TwoBody == kFALSE) {
	for (Int_t MN = InitialMass; MN < Mass; MN += step) {
	  BR = factor*ComputeBR(MN, Mass1, Mass2, 0., 0., Prod, TwoBody);  // factor to be removed/changed when taking into account coupling
	  xMN[MN] = MN/1000.;
	  yBR[MN] = BR;
	}
      }
    }
  }
  
  TGraph* gr = new TGraph(Mass, xMN, yBR);
  
  gr->SetNameTitle(Title.c_str(), Title.c_str());
  gr->SetLineWidth(5);
  gr->SetFillColor(0);
  
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

// Function for Dalitz plots of 3-body N production

void PhaseSpace(Double_t Mass1, Double_t Mass3, Double_t Mass4, std::string Title, Int_t model) {

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

    Double_t wMax = 0.;
    
    for (Int_t n = 0; n < 1000000; n++) {
      Double_t W = Event.Generate();
      TLorentzVector *p1 = Event.GetDecay(0);
      TLorentzVector *p2 = Event.GetDecay(1);
      TLorentzVector *p3 = Event.GetDecay(2);
      TLorentzVector p12 = *p1 + *p2;
      TLorentzVector p13 = *p1 + *p3;

      if (W > wMax)
	wMax = W;
      
      h2->Fill(p12.M2(), p13.M2());
    }
    
    cout<<Title<<" "<<wMax<<endl;

    h2->Draw("colz");
    h2->SetTitle(Form("%s, m_{N} = %.1f GeV, model %i",Title.c_str(), Mass2, model));
    h2->GetXaxis()->SetTitle("Neutrino-meson invariant mass [GeV^{2}/c^{4}]");
    h2->GetYaxis()->SetTitle("Neutrino-lepton invariant mass [GeV^{2}/c^{4}]");
    h2->GetXaxis()->SetRange(h2->FindFirstBinAbove(0., 1)-2, h2->FindLastBinAbove(0., 1)+2);
    h2->GetYaxis()->SetRange(h2->FindFirstBinAbove(0., 2)-2, h2->FindLastBinAbove(0., 2)+2);
    gPad->SetLogz();
    
    if (model == 1)
      gPad->SaveAs(Form("/home/li/Desktop/HeavyNeutrino/DalitzPlots/Model I/%s_%.1f_model%i.png", Title.c_str(), Mass2, model));
    else if (model == 2)
      gPad->SaveAs(Form("/home/li/Desktop/HeavyNeutrino/DalitzPlots/Model II/%s_%.1f_model%i.png", Title.c_str(), Mass2, model));
    else if (model == 3)
      gPad->SaveAs(Form("/home/li/Desktop/HeavyNeutrino/DalitzPlots/Model III/%s_%.1f_model%i.png", Title.c_str(), Mass2, model));
  }
}

// Making multigraph

TGraph* SumAllGraphs(TMultiGraph* M) {

  TList* list = M->GetListOfGraphs();
  TGraph* g = new TGraph();

  for (Int_t i = 0; i < ((TGraph*)list->At(0))->GetN(); ++i) {
    Double_t val = 0.;
    Double_t x, y;
    for (Int_t j = 0; j < list->GetEntries(); ++j) {
      TGraph* gg = (TGraph*)list->At(j);
      gg->GetPoint(i, x, y);
      val += y;
    }
    g->SetPoint(i, x, val);
  }

  return g;
}

// Function to call al Dalitz plots

void AllDalitz(Int_t model) {

  PhaseSpace(D,  K0,     e,  "D #rightarrow K^{0}eN"  ,     model);  // 3-body HNL production: Dalitz plots
  PhaseSpace(D,  pi0,    e,  "D #rightarrow #pi^{0}eN" ,    model);
  PhaseSpace(D0, K,      e,  "D^{0} #rightarrow KeN"  ,     model);
  PhaseSpace(D0, pi,     e,  "D^{0} #rightarrow #pieN" ,    model);
  PhaseSpace(D,  K0,     mu, "D #rightarrow K^{0}#muN" ,    model);
  PhaseSpace(D,  pi0,    mu, "D #rightarrow #pi^{0}#muN",   model);
  PhaseSpace(D0, K,      mu, "D^{0} #rightarrow K#muN" ,    model);
  PhaseSpace(D0, pi,     mu, "D^{0} #rightarrow #pi#muN",   model);
  PhaseSpace(D,  K0Star, e,  "D #rightarrow K^{0*}eN" ,     model);
  PhaseSpace(D0, KStar,  e,  "D^{0} #rightarrow K^{*}eN" ,  model);
  PhaseSpace(D,  K0Star, mu, "D #rightarrow K^{0*}#muN",    model);
  PhaseSpace(D0, KStar,  mu, "D^{0} #rightarrow K^{*}#muN", model);
}

// Function to call all N production modes

void AllProd (Int_t model, TMultiGraph* M) {

  MassScan(D,   e,      0., 1, 1, 0, "D #rightarrow Ne",                M, 1.);//sigmacc*2.*ffD);  // HNL production via two-body decay
  MassScan(D,   mu,     0., 1, 1, 0, "D #rightarrow N#mu",              M, 1.);//sigmacc*2.*ffD);
  MassScan(DS,  e,      0., 1, 1, 0, "D_{S} #rightarrow Ne",            M, 1.);//sigmacc*2.*ffDS);
  MassScan(DS,  mu,     0., 1, 1, 0, "D_{S} #rightarrow N#mu",          M, 1.);//sigmacc*2.*ffDS);
  MassScan(DS,  tau,    0., 1, 1, 0, "D_{S} #rightarrow N#tau",         M, 1.);//sigmacc*2.*ffDS);
  MassScan(tau, pi,     0., 1, 1, 0, "#tau #rightarrow N#pi",           M, 1.);//sigmacc*2.*ffDS);
  MassScan(tau, K,      0., 1, 1, 0, "#tau #rightarrow NK",             M, 1.);//sigmacc*2.*ffDS);
  MassScan(tau, rho,    0., 1, 1, 0, "#tau #rightarrow N#rho",          M, 1.);//sigmacc*2.*ffDS);
  MassScan(D,   K0,     e,  1, 0, 0, "D #rightarrow K^{0}eN",           M, 1.);//sigmacc*2.*ffD);  // HNL production via three-body decay
  MassScan(D,   pi0,    e,  1, 0, 0, "D #rightarrow #pi^{0}eN",         M, 1.);//sigmacc*2.*ffD);
  MassScan(D0,  K,      e,  1, 0, 0, "D^{0} #rightarrow KeN",           M, 1.);//sigmacc*2.*ffD0);
  MassScan(D0,  pi,     e,  1, 0, 0, "D^{0} #rightarrow #pieN",         M, 1.);//sigmacc*2.*ffD0);
  MassScan(D,   K0,     mu, 1, 0, 0, "D #rightarrow K^{0}#muN",         M, 1.);//sigmacc*2.*ffD);
  MassScan(D,   pi0,    mu, 1, 0, 0, "D #rightarrow #pi^{0}#muN",       M, 1.);//sigmacc*2.*ffD);
  MassScan(D0,  K,      mu, 1, 0, 0, "D^{0} #rightarrow K#muN",         M, 1.);//sigmacc*2.*ffD0);
  MassScan(D0,  pi,     mu, 1, 0, 0, "D^{0} #rightarrow #pi#muN",       M, 1.);//sigmacc*2.*ffD0);
  MassScan(D,   K0Star, e,  1, 0, 0, "D #rightarrow K^{0*}eN",          M, 1.);//sigmacc*2.*ffD);
  MassScan(D0,  KStar,  e,  1, 0, 0, "D^{0} #rightarrow K^{*}eN",       M, 1.);//sigmacc*2.*ffD0);
  MassScan(D,   K0Star, mu, 1, 0, 0, "D #rightarrow K^{0*}#muN",        M, 1.);//sigmacc*2.*ffD);
  MassScan(D0,  KStar,  mu, 1, 0, 0, "D^{0} #rightarrow K^{*}#muN",     M, 1.);//sigmacc*2.*ffD0);
  MassScan(tau, 0.1,    e,  1, 0, 0, "#tau #rightarrow Ne#nu_{#tau}",   M, 1.);//sigmacc*2.*ffDS);
  MassScan(tau, 0.01,   e,  1, 0, 0, "#tau #rightarrow Ne#nu_{e}",      M, 1.);//sigmacc*2.*ffDS);
  MassScan(tau, 0.1,    mu, 1, 0, 0, "#tau #rightarrow N#mu#nu_{#tau}", M, 1.);//sigmacc*2.*ffDS);
  MassScan(tau, 0.01,   mu, 1, 0, 0, "#tau #rightarrow N#mu#nu_{#mu}",  M, 1.);//sigmacc*2.*ffDS);

  M->Draw("AC");
  M->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
  M->GetYaxis()->SetTitle("BR");//"f*BR");
  M->GetXaxis()->SetTitleOffset(1.2);
  M->GetYaxis()->SetTitleOffset(1.2);
  M->GetXaxis()->SetTitleSize(labelSize);
  M->GetYaxis()->SetTitleSize(labelSize);
  M->GetXaxis()->SetLabelSize(labelSize);
  M->GetYaxis()->SetLabelSize(labelSize);
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  M->SetMinimum(1.E-8);
  M->SetMaximum(1.1);
  M->GetXaxis()->SetLimits(0., 2.2);
  gPad->BuildLegend(0.818, 0.223, 0.984, 0.881);
  gPad->Update();
  gPad->Modified();
  gPad->Write();
  
  gPad->SaveAs(Form("/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/NProdGraph_%i.pdf", model));
  gPad->SaveAs(Form("/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/NProdGraph_%i.png", model));

  new TCanvas;
  TGraph* sum = SumAllGraphs(M);
  sum->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
  sum->GetYaxis()->SetTitle("Sum of BRs");    
  sum->SetTitle("Sum of N production modes vs N mass, model II (1:16:3.8)");
  sum->GetXaxis()->SetTitleOffset(0.9);
  sum->GetYaxis()->SetTitleOffset(1.);
  sum->GetXaxis()->SetTitleSize(labelSize);
  sum->GetYaxis()->SetTitleSize(labelSize);
  sum->GetXaxis()->SetLabelSize(labelSize);
  sum->GetYaxis()->SetLabelSize(labelSize);
  gPad->SetGridx();
  gPad->SetGridy();
  sum->GetXaxis()->SetLimits(0., 2.2);
  sum->SetLineColor(2);
  sum->SetLineWidth(3);
  sum->Draw();
  gPad->Update();
  gPad->Write();
  gPad->SaveAs(Form("/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/NSumProdGraph_%i.pdf", model));
  gPad->SaveAs(Form("/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/NSumProdGraph_%i.png", model));
}

// Function to call all N decay modes

void AllDecay(Int_t model, TMultiGraph* M) {

  MassScan(0.,        0.,  0., 0, 0, 0, "N #rightarrow #nu#nu#nu",   M, 1.);  // HNL decay via two- and three-body decay
  MassScan(e,         e,   0., 0, 0, 0, "N #rightarrow ee#nu",       M, 1.);
  MassScan(e,         mu,  0., 0, 0, 0, "N #rightarrow e#mu#nu",     M, 1.);
  MassScan(pi0,       0.,  0., 0, 1, 0, "N #rightarrow #pi^{0}#nu",  M, 1.);
  MassScan(pi,        e,   0., 0, 1, 0, "N #rightarrow #pie",        M, 1.);
  MassScan(mu,        mu,  0., 0, 0, 0, "N #rightarrow #mu#mu#nu",   M, 1.);
  MassScan(pi,        mu,  0., 0, 1, 0, "N #rightarrow #pi#mu",      M, 1.);
  MassScan(K,         e,   0., 0, 1, 0, "N #rightarrow Ke",          M, 1.);
  MassScan(eta,       0.,  0., 0, 1, 0, "N #rightarrow #eta#nu",     M, 1.);
  MassScan(K,         mu,  0., 0, 1, 0, "N #rightarrow K#mu",        M, 1.);
  MassScan(rho0,      0.,  0., 0, 1, 0, "N #rightarrow #rho^{0}#nu", M, 1.);
  MassScan(rho,       e,   0., 0, 1, 0, "N #rightarrow #rhoe",       M, 1.);
  MassScan(rho,       mu,  0., 0, 1, 0, "N #rightarrow #rho#mu",     M, 1.);
  MassScan(etaprime,  0.,  0., 0, 1, 0, "N #rightarrow #eta'#nu",    M, 1.);
  MassScan(e,         tau, 0., 0, 0, 0, "N #rightarrow e#tau#nu",    M, 1.);
  MassScan(mu,        tau, 0., 0, 0, 0, "N #rightarrow #mu#tau#nu",  M, 1.);
  MassScan(pi,        tau, 0., 0, 1, 0, "N #rightarrow #pi#tau",     M, 1.);
  MassScan(K,         tau, 0., 0, 1, 0, "N #rightarrow K#tau",       M, 1.);
  MassScan(rho,       tau, 0., 0, 1, 0, "N #rightarrow #rho#tau",    M, 1.);
  MassScan(tau,       tau, 0., 0, 0, 0, "N #rightarrow #tau#tau#nu", M, 1.);
  
  M->Draw("AC");
  M->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
  M->GetYaxis()->SetTitle("BR");
  M->GetXaxis()->SetTitleOffset(1.2);
  M->GetYaxis()->SetTitleOffset(1.2);
  M->GetXaxis()->SetTitleSize(labelSize);
  M->GetYaxis()->SetTitleSize(labelSize);
  M->GetXaxis()->SetLabelSize(labelSize);
  M->GetYaxis()->SetLabelSize(labelSize);
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  M->SetMinimum(1.E-7);
  M->SetMaximum(1.5);
  M->GetXaxis()->SetLimits(0., 5.);
  gPad->BuildLegend(0.841, 0.218, 0.985, 0.862);
  gPad->Update();
  gPad->Modified();
  gPad->Write();
 
  gPad->SaveAs(Form("/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/NDecayGraph_%i.pdf", model));
  gPad->SaveAs(Form("/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/NDecayGraph_%i.png", model));
}

// Function to call all N partial decay widths

void AllGamma(Int_t model, TMultiGraph* M) {

  MassScan(0.,        0.,  0., 0, 0, 1, "N #rightarrow #nu#nu#nu",   M, 1.);  // HNL decay via two- and three-body decay
  MassScan(e,         e,   0., 0, 0, 1, "N #rightarrow ee#nu",       M, 1.);
  MassScan(e,         mu,  0., 0, 0, 1, "N #rightarrow e#mu#nu",     M, 1.);
  MassScan(pi0,       0.,  0., 0, 1, 1, "N #rightarrow #pi^{0}#nu",  M, 1.);
  MassScan(pi,        e,   0., 0, 1, 1, "N #rightarrow #pie",        M, 1.);
  MassScan(mu,        mu,  0., 0, 0, 1, "N #rightarrow #mu#mu#nu",   M, 1.);
  MassScan(pi,        mu,  0., 0, 1, 1, "N #rightarrow #pi#mu",      M, 1.);
  MassScan(K,         e,   0., 0, 1, 1, "N #rightarrow Ke",          M, 1.);
  MassScan(eta,       0.,  0., 0, 1, 1, "N #rightarrow #eta#nu",     M, 1.);
  MassScan(K,         mu,  0., 0, 1, 1, "N #rightarrow K#mu",        M, 1.);
  MassScan(rho0,      0.,  0., 0, 1, 1, "N #rightarrow #rho^{0}#nu", M, 1.);
  MassScan(rho,       e,   0., 0, 1, 1, "N #rightarrow #rhoe",       M, 1.);
  MassScan(rho,       mu,  0., 0, 1, 1, "N #rightarrow #rho#mu",     M, 1.);
  MassScan(etaprime,  0.,  0., 0, 1, 1, "N #rightarrow #eta'#nu",    M, 1.);
  MassScan(e,         tau, 0., 0, 0, 1, "N #rightarrow e#tau#nu",    M, 1.);

  MassScan(mu,        tau, 0., 0, 0, 1, "N #rightarrow #mu#tau#nu",  M, 1.);
  MassScan(pi,        tau, 0., 0, 1, 1, "N #rightarrow #pi#tau",     M, 1.);
  MassScan(K,         tau, 0., 0, 1, 1, "N #rightarrow K#tau",       M, 1.);
  MassScan(rho,       tau, 0., 0, 1, 1, "N #rightarrow #rho#tau",    M, 1.);
  MassScan(tau,       tau, 0., 0, 0, 1, "N #rightarrow #tau#tau#nu", M, 1.);

  M->Draw("AC");
  M->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
  M->GetYaxis()->SetTitle("Partial decay width [MeV]");
  M->GetXaxis()->SetTitleOffset(1.2);
  M->GetYaxis()->SetTitleOffset(1.2);
  M->GetXaxis()->SetTitleSize(labelSize);
  M->GetYaxis()->SetTitleSize(labelSize);
  M->GetXaxis()->SetLabelSize(labelSize);
  M->GetYaxis()->SetLabelSize(labelSize);
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  M->SetMinimum(1.E-18);
  M->SetMaximum(1.E-4);
  M->GetXaxis()->SetLimits(0., 10.);
  gPad->BuildLegend(0.841, 0.218, 0.985, 0.862);
  gPad->Update();
  gPad->Modified();
  gPad->Write();
  
  gPad->SaveAs(Form("/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/NWidthGraph_%i.pdf", model));
  gPad->SaveAs(Form("/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/NWidthGraph_%i.png", model));
}

// Main

Int_t AllInOnePlot(Int_t mode, Int_t model) {

  // mode = 0 (Dalitz), 1 (prod), 2 (decay), 3 (width); model = 1, 2, 3 (Shaposhnikov)
  
  TCanvas *c = new TCanvas();

  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.2);
  //c->SetWindowSize(20000., 12000.);

  TGaxis::SetMaxDigits(2);

  if (model == 1 || model == 2 || model == 3)
    SetModel(model);
  else {
    cout<<"[GeneralPlots] Unknown model:"<<endl;
    cout<<"1: 52:1:1"<<endl;
    cout<<"2: 1:16:3.8"<<endl;
    cout<<"3: 0.061:1:4.3"<<endl;
    exit(1);
  }
    
  TMultiGraph* Mprod = new TMultiGraph("Mprod", "");
  Mprod->SetName("Mprod");
  TMultiGraph* Mdecay = new TMultiGraph("Mdecay", "");
  Mdecay->SetName("Mdecay");
  TMultiGraph* Mgamma = new TMultiGraph("Mgamma", "");
  Mdecay->SetName("Mgamma");
  std::string HistoTitle = "";
  
  if (model == 1) 
    HistoTitle = "I (52:1:1)";
  else if (model == 2)
    HistoTitle = "II (1:16:3.8)";
  else if (model == 3)
    HistoTitle = "III (0.061:1:4.3)";

  Mprod->SetTitle(("N production modes vs N mass, model " + HistoTitle).c_str());
  Mdecay->SetTitle(("N decay modes vs N mass, model " + HistoTitle).c_str());
  Mgamma->SetTitle(("N partial decay widths vs N mass, model " + HistoTitle).c_str());
    
  TGaxis::SetMaxDigits(2);

  if (mode == 0) {
    AllDalitz(model);
  }
  else if (mode == 1) {
    AllProd(model, Mprod);
  }
  else if (mode == 2) {
    AllDecay(model, Mdecay);
  }
  else if (mode == 3) {
    AllGamma(model, Mgamma);
  }
  else {
    cout<<"[GeneralPlots] Unknown mode:"<<endl;
    cout<<"0: Dalitz plots for 3-body N production"<<endl;
    cout<<"1: BR plots for N production modes"<<endl;
    cout<<"2: BR plots for N decay modes"<<endl;
    cout<<"3: Partial decay width plots for N decay modes"<<endl;
    exit(1);
  }
  
  exit(0);
}
