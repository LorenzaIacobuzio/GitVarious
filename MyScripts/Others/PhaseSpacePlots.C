Double_t PhaseSpace(Double_t mesonMass, Double_t MN, Double_t leptonMass) {
  
  Double_t phaseSpace = 0.;

  phaseSpace = TMath::Power(mesonMass*mesonMass-MN*MN-leptonMass*leptonMass,2) - 4.*MN*MN*leptonMass*leptonMass;

  return phaseSpace;
}

Double_t PhaseSpaceFactor(Double_t mesonMass, Double_t MN, Double_t leptonMass, Double_t phaseSpace) {

  Double_t factor = 0.;

  if(phaseSpace > 0.)
    factor = (mesonMass*mesonMass*(MN*MN + leptonMass*leptonMass)-TMath::Power(MN*MN-leptonMass*leptonMass,2))*TMath::Power(phaseSpace,0.5)/(leptonMass*leptonMass*TMath::Power(mesonMass*mesonMass-leptonMass*leptonMass,2));

  return factor;
}

Double_t GammaPiE(Double_t MN) {

  Double_t Mpi = 139.6;
  Double_t Me = 0.511;
  Double_t GF = 1.1E-11; // MeV^-2                             
  Double_t fPi = 130.41; // MeV                                   
  Double_t gamma_pi_e = 0.;

  gamma_pi_e = GF*GF*fPi*fPi/(8*TMath::Pi()*MN*MN)*(Me*Me + MN*MN)*(MN*MN - Mpi*Mpi - Me*Me)*TMath::Power(MN*MN*MN*MN + Mpi*Mpi*Mpi*Mpi + Me*Me*Me*Me - 2.*MN*Me*MN*Me - 2.*MN*Mpi*MN*Mpi - 2.*Mpi*Me*Mpi*Me,0.5)/(2*MN);

  return gamma_pi_e;
}

Double_t GammaPiMu(Double_t MN) {

  Double_t Mpi = 139.6;
  Double_t Mmu = 105.;
  Double_t GF = 1.1E-11;
  Double_t fPi = 130.41;
  Double_t gamma_pi_mu = 0.;

  gamma_pi_mu = GF*GF*fPi*fPi/(8*TMath::Pi()*MN*MN)*(Mmu*Mmu + MN*MN)*(MN*MN - Mpi*Mpi - Mmu*Mmu)*TMath::Power(MN*MN*MN*MN + Mpi*Mpi*Mpi*Mpi + Mmu*Mmu*Mmu*Mmu - 2.*MN*Mmu*MN*Mmu - 2.*MN*Mpi*MN*Mpi - 2.*Mpi*Mmu*Mpi*Mmu,0.5)/(2*MN);

  return gamma_pi_mu;
}

void DrawGraph(TGraph* g) {

  g->GetXaxis()->SetTitle("N mass [GeV]");
  g->GetYaxis()->SetTitle("Fl");
  g->SetLineWidth(3);
  g->Draw("AC");
  g->Write();
}

void HNLPlots() {

  Double_t Me = 0.511;
  Double_t Mmu = 105.;
  Double_t Mpi = 139.6;
  Double_t DMass = 1870.;
  Double_t DSMass = 1969.;
  Double_t PS = 0.;
  Double_t PSF = 0.;
  Double_t PiE = 0.;
  Double_t PiMu = 0.;
  Int_t Mass = 300;
  Double_t xMN[Mass];
  Double_t yPSF[Mass];
  Double_t yGammaPiE[Mass];
  Double_t yGammaPiMu[Mass];

  TFile *f = new TFile("Plots.root","RECREATE");

  for (Int_t MN = 0; MN < Mass; MN++) {
    PS = PhaseSpace(DMass, MN, Me);
    PSF = PhaseSpaceFactor(DMass, MN, Me, PS);
    xMN[MN] = MN;
    yPSF[MN] = PSF;
  }

  TGraph* gr1 = new TGraph(Mass, xMN, yPSF);
  gr1->SetTitle("Phase space factor vs N mass, for electrons");
  DrawGraph(gr1);

  for (Int_t MN = 0; MN < Mass; MN++) {
    PS = PhaseSpace(DMass, MN, Mmu);
    PSF = PhaseSpaceFactor(DMass, MN, Mmu, PS);
    xMN[MN] = MN;
    yPSF[MN] = PSF;
  }

  TGraph* gr2 = new TGraph(Mass, xMN, yPSF);
  gr2->SetTitle("Phase space factor vs N mass, for muons");
  DrawGraph(gr2);

  for (Int_t MN = 0; MN < Mass; MN++) {
    PS = PhaseSpace(DSMass, MN, Me);
    PSF = PhaseSpaceFactor(DSMass, MN, Me, PS);
    xMN[MN] = MN;
    yPSF[MN] = PSF;
  }

  TGraph* gr3 = new TGraph(Mass, xMN, yPSF);
  gr3->SetTitle("Phase space factor vs N mass, for electrons");
  DrawGraph(gr3);

  for (Int_t MN = 0; MN < Mass; MN++) {
    PS = PhaseSpace(DSMass, MN, Mmu);
    PSF = PhaseSpaceFactor(DSMass, MN, Mmu, PS);
    xMN[MN] = MN;
    yPSF[MN] = PSF;
  }

  TGraph* gr4 = new TGraph(Mass, xMN, yPSF);
  gr4->SetTitle("Phase space factor vs N mass, for muons");
  DrawGraph(gr4);

  for (Int_t MN = 0; MN < Mass; MN++) {
    if (MN > Mpi + Me)
      PiE = GammaPiE(MN)*1.E13;
    else
      PiE = 0.;
    if (MN > Mpi + Mmu)
      PiMu = GammaPiMu(MN)*1.E13;
    else
      PiMu = 0.;

    xMN[MN] = MN;
    yGammaPiE[MN] = PiE;
    yGammaPiMu[MN] = PiMu;
  }

  TGraph* gr5 = new TGraph(Mass, xMN, yGammaPiE);

  gr5->GetXaxis()->SetTitle("N mass [GeV]");
  gr5->GetYaxis()->SetTitle("GammaPiE");
  gr5->SetTitle("Decay amplitude of N into pion-electron pairs vs N mass");
  gr5->SetLineWidth(3);
  gr5->Draw("AC");
  gr5->Write();

  TGraph* gr6 = new TGraph(Mass, xMN, yGammaPiMu);

  gr6->GetXaxis()->SetTitle("N mass [GeV]");
  gr6->GetYaxis()->SetTitle("GammaPiMu");
  gr6->SetTitle("Decay amplitude of N into pion-muon pairs vs N mass");
  gr6->SetLineWidth(3);
  gr6->Draw("AC");
  gr6->Write();

  f->Write();
}
