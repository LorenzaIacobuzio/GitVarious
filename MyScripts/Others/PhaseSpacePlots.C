Double_t PhaseSpace(Double_t mesonMass, Double_t MN, Double_t leptonMass) {

  Double_t phaseSpace = TMath::Power(mesonMass*mesonMass-MN*MN-leptonMass*leptonMass,2) - 4.*MN*MN*leptonMass*leptonMass;

  return phaseSpace;
}

Double_t PhaseSpaceFactor(Double_t mesonMass, Double_t MN, Double_t leptonMass, Double_t phaseSpace) {

  Double_t factor = 0.;

  if(phaseSpace > 0.) {
    factor = ((mesonMass*mesonMass*(MN*MN + leptonMass*leptonMass)-TMath::Power(MN*MN-leptonMass*leptonMass,2))*TMath::Power(phaseSpace,0.5))/(leptonMass*leptonMass*TMath::Power(mesonMass*mesonMass-leptonMass*leptonMass,2));
  }

  return factor;
}

Double_t GammaPiE(Double_t MN) {

  Double_t Mpi = 139.6;
  Double_t Me = 0.511;
  Double_t GF = 1.1E-11; // MeV^-2                             
  Double_t fPi = 130.41; // MeV                                   
  Double_t gamma_pi_e = 0.;
  
  gamma_pi_e = (GF*GF*fPi*fPi/(8*TMath::Pi()*MN*MN)*(Me*Me + MN*MN)*(MN*MN - Mpi*Mpi - Me*Me)*TMath::Power(MN*MN*MN*MN + Mpi*Mpi*Mpi*Mpi + Me*Me*Me*Me - 2.*MN*Me*MN*Me - 2.*MN*Mpi*MN*Mpi - 2.*Mpi*Me*Mpi*Me,0.5))/(2*MN);

  return gamma_pi_e;
}

Double_t GammaPiMu(Double_t MN) {

  Double_t Mpi = 139.6;
  Double_t Mmu = 105.;
  Double_t GF = 1.1E-11;
  Double_t fPi = 130.41;
  Double_t gamma_pi_mu = 0.;

  gamma_pi_mu = (GF*GF*fPi*fPi/(8*TMath::Pi()*MN*MN)*(Mmu*Mmu + MN*MN)*(MN*MN - Mpi*Mpi - Mmu*Mmu)*TMath::Power(MN*MN*MN*MN + Mpi*Mpi*Mpi*Mpi + Mmu*Mmu*Mmu*Mmu - 2.*MN*Mmu*MN*Mmu - 2.*MN*Mpi*MN*Mpi - 2.*Mpi*Mmu*Mpi*Mmu,0.5))/(2*MN);

  return gamma_pi_mu;
}

Double_t ComputeBR(Double_t mesonMass, Double_t leptonMass, Double_t Factor) {

  Double_t Me = 0.511;
  Double_t Mmu = 105.;
  Double_t DMass = 1870.;
  Double_t DSMass = 1969.;
  Double_t De2BR = 9.2E-9;
  Double_t Dm2BR = 3.74E-4;
  Double_t DSe2BR = 1.4E-7;
  Double_t DSm2BR = 5.56E-3;
  Double_t brt = 0.;

  if (Factor > 0.) {
    if (mesonMass == DMass && leptonMass == Me)
      brt = De2BR*Factor;
    
    else if (mesonMass == DMass && leptonMass == Mmu)
      brt = Dm2BR*Factor;
    
    else if (mesonMass == DSMass && leptonMass == Me)
      brt = DSe2BR*Factor;
    
    else if (mesonMass == DSMass && leptonMass == Mmu)
      brt = DSm2BR*Factor;
  }
  
  return brt;
}

Double_t ComputeTotalBR(Double_t mesonMass, Double_t MN) {

  Double_t Me = 0.511;
  Double_t Mmu = 105.;
  Double_t DMass = 1870.;
  Double_t DSMass = 1969.;
  Double_t De2BR = 9.2E-9;
  Double_t Dm2BR = 3.74E-4;
  Double_t DSe2BR = 1.4E-7;
  Double_t DSm2BR = 5.56E-3;
  Double_t brt = 0.;

  Double_t phaseSpaceE = PhaseSpace(mesonMass, MN, Me);
  Double_t phaseSpaceMu = PhaseSpace(mesonMass, MN, Mmu);
  Double_t FactorE = PhaseSpaceFactor(mesonMass, MN, Me, phaseSpaceE);
  Double_t FactorMu = PhaseSpaceFactor(mesonMass, MN, Mmu, phaseSpaceMu);
  
  if (FactorE > 0. && FactorMu > 0.) {
    if (mesonMass == DMass)
      brt = De2BR*FactorE/(De2BR*FactorE+Dm2BR*FactorMu);
    
    else if (mesonMass == DSMass)
      brt = DSe2BR*FactorE/(DSe2BR*FactorE+DSm2BR*FactorMu);
  }
  
  return brt;
}

void PhaseSpacePlots() {

  Double_t Me = 0.511;
  Double_t Mmu = 105.;
  Double_t Mpi = 139.6;
  Double_t DMass = 1870.;
  Double_t DSMass = 1969.;
  Double_t PS = 0.;
  Double_t PSF = 0.;
  Double_t PiE = 0.;
  Double_t PiMu = 0.;
  Double_t BR = 0.;
  Double_t BRratio = 0.;
  const int Mass = 10000;
  const int step = 1;
  const int Masses = Mass/step;
  Double_t xMN[Masses];
  Double_t yPSF[Masses];
  Double_t yGammaPiE[Masses];
  Double_t yGammaPiMu[Masses];
  Double_t yBR[Masses];
  Double_t yBRratio[Masses];
  TString name = "";
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;
  
  TFile *f = new TFile("/home/li/Desktop/Plots.root","RECREATE");

  for (Int_t MN = 0; MN < Mass; MN += step) {
    PS = PhaseSpace(DMass, MN, Me);
    PSF = PhaseSpaceFactor(DMass, MN, Me, PS);
    BR = ComputeBR(DMass, Me, PSF);
    BRratio = ComputeTotalBR(DMass, MN);
    xMN[MN] = MN;
    if (PSF >= 0.)
      yPSF[MN] = PSF;
    yBR[MN] = BR;
    yBRratio[MN] = BRratio;
  }

  TGraph* gr1 = new TGraph(Mass, xMN, yPSF);
  gr1->SetName("gr1");

  TGraph* gr7 = new TGraph(Mass, xMN, yBR);
  gr7->SetName("gr7");

  TGraph* gr11 = new TGraph(Mass, xMN, yBRratio);
  gr11->SetName("gr11");
  
  for (Int_t MN = 0; MN < Mass; MN += step) {
    PS = PhaseSpace(DMass, MN, Mmu);
    PSF = PhaseSpaceFactor(DMass, MN, Mmu, PS);
    BR = ComputeBR(DMass, Mmu, PSF);
    xMN[MN] = MN;
    if (PSF >= 0.)
      yPSF[MN] = PSF;
    yBR[MN] = BR;
  }

  TGraph* gr2 = new TGraph(Mass, xMN, yPSF);
  gr2->SetName("gr2");

  TGraph* gr8 = new TGraph(Mass, xMN, yBR);
  gr8->SetName("gr8");

  for (Int_t MN = 0; MN < Mass; MN += step) {
    PS = PhaseSpace(DSMass, MN, Me);
    PSF = PhaseSpaceFactor(DSMass, MN, Me, PS);
    BR = ComputeBR(DSMass, Me, PSF);
    BRratio = ComputeTotalBR(DSMass, MN);
    xMN[MN] = MN;
    if (PSF >= 0.)
      yPSF[MN] = PSF;
    yBR[MN] = BR;
    yBRratio[MN] = BRratio;
  }

  TGraph* gr3 = new TGraph(Mass, xMN, yPSF);
  gr3->SetName("gr3");

  TGraph* gr9 = new TGraph(Mass, xMN, yBR);
  gr9->SetName("gr9");

  TGraph* gr12 = new TGraph(Mass, xMN, yBRratio);
  gr12->SetName("gr12");

  for (Int_t MN = 0; MN < Mass; MN += step) {
    PS = PhaseSpace(DSMass, MN, Mmu);
    PSF = PhaseSpaceFactor(DSMass, MN, Mmu, PS);
    BR = ComputeBR(DSMass, Mmu, PSF);
    xMN[MN] = MN;
    if (PSF >= 0.)
      yPSF[MN] = PSF;
    yBR[MN] = BR;
  }

  TGraph* gr4 = new TGraph(Mass, xMN, yPSF);
  gr4->SetName("gr4");

  TGraph* gr10 = new TGraph(Mass, xMN, yBR);
  gr10->SetName("gr10");

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

  TMultiGraph* M = new TMultiGraph("M", "BR times phase space vs N mass");
  M->SetName("M");
  TMultiGraph* M2 = new TMultiGraph("M2", "BR vs N mass");
  M2->SetName("M2");
  TMultiGraph* M3 = new TMultiGraph("M3", "N semi-leptonic partial decay widths vs N mass");
  M3->SetName("M3");

  gr1->SetLineColor(2);
  gr1->SetLineWidth(3);
  gr2->SetLineColor(4);
  gr2->SetLineWidth(3);
  gr3->SetLineColor(6);
  gr3->SetLineWidth(3);
  gr4->SetLineColor(9);
  gr4->SetLineWidth(3);

  gr11->SetLineColor(2);
  gr11->SetLineWidth(3);
  gr12->SetLineColor(4);
  gr12->SetLineWidth(3);

  gr5->SetLineColor(2);
  gr5->SetLineWidth(5);
  gr6->SetLineColor(kGreen+1);
  gr6->SetLineWidth(2);

  gr7->SetLineColor(2);
  gr7->SetLineWidth(5);
  gr8->SetLineColor(kGreen+1);
  gr8->SetLineWidth(2);
  gr9->SetLineColor(6);
  gr9->SetLineWidth(3);
  gr10->SetLineColor(4);
  gr10->SetLineWidth(3);

  auto legend = new TLegend(0.1,0.65,0.4,0.9);
  legend->AddEntry(gr7,  "BR(D->nue)*Fe", "l");
  legend->AddEntry(gr8,  "BR(D->numu)*Fmu", "l");
  legend->AddEntry(gr9,  "BR(DS->nue)*Fe", "l");
  legend->AddEntry(gr10, "BR(DS->numu)*Fmu", "l");

  auto legend2 = new TLegend(0.1,0.65,0.45,0.9);
  legend2->AddEntry(gr11,  "BRe/tot for D", "l");
  legend2->AddEntry(gr12,  "BRe/tot for DS", "l");

  auto legend3 = new TLegend(0.1,0.75,0.35,0.9);
  legend3->AddEntry(gr5,  "N->pie", "l");
  legend3->AddEntry(gr6,  "N->pimu", "l");

  M->Add(gr7);
  M->Add(gr8);
  M->Add(gr9);
  M->Add(gr10);
  M->Draw("AC");
  M->GetXaxis()->SetTitle("N mass [MeV]");
  M->GetXaxis()->SetRangeUser(-50, 2100);
  M->GetYaxis()->SetTitle("BR*F");
  M->GetYaxis()->SetTitleOffset(1.);
  M->GetXaxis()->SetTitleOffset(0.9);
  M->GetXaxis()->SetTitleSize(labelSize);
  M->GetYaxis()->SetTitleSize(labelSize);
  M->GetXaxis()->SetLabelSize(labelSize);
  M->GetYaxis()->SetLabelSize(labelSize);
  TGaxis::SetMaxDigits(2);
  gPad->Update();
  M->Write();

  gr1->GetYaxis()->SetTitle("F(D,DS->Ne)");
  gr1->SetTitle("Phase space factor vs N mass");
  gr1->GetXaxis()->SetTitle("N mass [MeV]");
  gr1->GetXaxis()->SetRangeUser(0, 2000);
  gr1->GetXaxis()->SetTitleOffset(0.9);
  gr1->GetXaxis()->SetTitleSize(labelSize);
  gr1->GetYaxis()->SetTitleSize(labelSize);
  gr1->GetXaxis()->SetLabelSize(labelSize);
  gr1->GetYaxis()->SetLabelSize(labelSize); 
  gr1->Draw("AC");
  gPad->Update();
  gr1->Write();

  gr2->GetYaxis()->SetTitle("F(D,DS->Nmu)");
  gr2->GetYaxis()->SetTitleOffset(1.);
  gr2->SetTitle("Phase space factor vs N mass");
  gr2->GetXaxis()->SetTitle("N mass [MeV]");
  gr2->GetXaxis()->SetRangeUser(0, 2000);
  gr2->GetXaxis()->SetTitleOffset(0.9);
  gr2->GetXaxis()->SetTitleSize(labelSize);
  gr2->GetYaxis()->SetTitleSize(labelSize);
  gr2->GetXaxis()->SetLabelSize(labelSize);
  gr2->GetYaxis()->SetLabelSize(labelSize); 
  gr2->Draw("AC");
  gPad->Update();
  gr2->Write();

  M2->Add(gr11);
  M2->Add(gr12);
  M2->Draw("AC");
  M2->GetXaxis()->SetTitle("N mass [MeV]");
  M2->GetXaxis()->SetRangeUser(-500, 2500);
  M2->GetYaxis()->SetTitle("BRe/tot");
  M2->GetXaxis()->SetTitleOffset(0.9);
  M2->GetXaxis()->SetTitleSize(labelSize);
  M2->GetYaxis()->SetTitleSize(labelSize);
  M2->GetXaxis()->SetLabelSize(labelSize);
  M2->GetYaxis()->SetLabelSize(labelSize);
  gPad->Update();
  M2->Write();

  M3->Add(gr5);
  M3->Add(gr6);
  M3->Draw("AC");
  M3->GetXaxis()->SetTitle("N mass [MeV]");
  M3->GetYaxis()->SetTitle("Decay width [MeV]");
  M3->GetYaxis()->SetTitleOffset(1.);
  M3->GetXaxis()->SetTitleOffset(0.9);
  M3->GetXaxis()->SetTitleSize(labelSize);
  M3->GetYaxis()->SetTitleSize(labelSize);
  M3->GetXaxis()->SetLabelSize(labelSize);
  M3->GetYaxis()->SetLabelSize(labelSize);
  gPad->Update();
  M3->Write();

  f->Write();
}
