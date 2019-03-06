void Save(TString path, TCanvas *c, TH1D* h, TString x, Double_t labelSize, Double_t titleSize) {

  h->Draw();
  h->GetXaxis()->SetTitle(x);
  h->SetFillColor(38);
  h->SetTitleSize(titleSize, "t");
  h->GetXaxis()->SetTitleSize(labelSize);
  h->GetXaxis()->SetLabelSize(labelSize);
  h->GetYaxis()->SetTitleSize(labelSize);
  h->GetYaxis()->SetLabelSize(labelSize);
  h->GetXaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat(0);
  gPad->Update();
  c->SaveAs(path + h->GetName() + ".pdf");
  c->SaveAs(path + h->GetName() + ".png");
}

void Save(TString path, TCanvas *c, TH2D* h, TString x, TString y, Double_t labelSize, Double_t titleSize) {

  TString name = h->GetName();

  if (name.Contains("Mass"))
    h->Draw("colz");
  else
    h->Draw("text");

  h->GetXaxis()->SetTitle(x);
  h->GetYaxis()->SetTitle(y);
  h->SetTitleSize(titleSize, "t");
  h->GetXaxis()->SetTitleSize(labelSize);
  h->GetYaxis()->SetTitleSize(labelSize);
  h->GetXaxis()->SetLabelSize(labelSize);
  h->GetYaxis()->SetLabelSize(labelSize);
  h->GetXaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat(0);
  gPad->Update();
  c->SaveAs(path + h->GetName() + ".pdf");
  c->SaveAs(path + h->GetName() + ".png");
}

void TreePlots(TString dir, TString histo1) {

  // dir = output dir, histo1 = histo to do cosmetics on

  TCanvas *c = new TCanvas();  
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;  
  TString path = "";
  Int_t counterComb = 0;
  Int_t counterPrompt = 0;
  Int_t counterPar = 0;

  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.25);
  c->SetTopMargin(0.15);
  c->SetGrid();
  c->RedrawAxis();

  if (dir != "")
    path = dir;
  else {
    if (histo1.Contains("2016"))
      path = "/Users/lorenza/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/2/2016/";
    else if (histo1.Contains("2017"))
      path = "/Users/lorenza/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/2/2017/";
    else if (histo1.Contains("Capped"))
      path = "/Users/lorenza/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/2/K3piCapped/";
    else if (!histo1.Contains("Capped") && histo1.Contains("K3pi"))
      path = "/Users/lorenza/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/2/K3pi/";

  //path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/2/HeavyNeutrino/";
  }
  
  TFile *f = TFile::Open(histo1);

  if (f == 0) {
    cout << "Error: cannot open " << histo1 << endl;
    return;
  }

  TTree* tree = (TTree*)f->Get("Passed");
  Double_t Weight;
  Double_t CHODTime1;
  Double_t CHODTime2;
  Double_t CDA;
  Double_t Zvertex;
  Double_t CDALine;
  Double_t ZCDALine;
  Double_t BeamlineDist;
  Double_t xSR;
  Double_t ySR;
  Double_t MuEoP;
  Double_t PiEoP;
  Double_t R;
  Double_t energyPi;
  Double_t energyMu;
  Double_t invMass;
  Double_t L0TPTime;
  TVector3 *Mom1 = new TVector3();
  TVector3 *Mom2 = new TVector3();
  TVector3 *TotMom = new TVector3();
  TVector3 *Vertex = new TVector3();
  TVector3 *threeMomPi = new TVector3();
  TVector3 *threeMomMu = new TVector3();
  Bool_t Target;
  Bool_t K3pi;
  Bool_t autoPass;
  Int_t Assoc;
  //TRecoCedarCandidate *KTAGcand;
  
  tree->SetBranchAddress("Weight", &Weight);
  tree->SetBranchAddress("CHODTime1", &CHODTime1);
  tree->SetBranchAddress("CHODTime2", &CHODTime2);
  tree->SetBranchAddress("CDA", &CDA);
  tree->SetBranchAddress("Zvertex", &Zvertex);
  tree->SetBranchAddress("CDALine", &CDALine);
  tree->SetBranchAddress("ZCDALine", &ZCDALine);
  tree->SetBranchAddress("BeamlineDist", &BeamlineDist);
  tree->SetBranchAddress("xSR", &xSR);
  tree->SetBranchAddress("ySR", &ySR);
  tree->SetBranchAddress("MuEoP", &MuEoP);
  tree->SetBranchAddress("PiEoP", &PiEoP);
  tree->SetBranchAddress("R", &R);
  tree->SetBranchAddress("energyPi", &energyPi);
  tree->SetBranchAddress("energyMu", &energyMu);
  tree->SetBranchAddress("invMass", &invMass);
  tree->SetBranchAddress("L0TPTime", &L0TPTime);
  tree->SetBranchAddress("Mom1", &Mom1);
  tree->SetBranchAddress("Mom2", &Mom2);
  tree->SetBranchAddress("TotMom", &TotMom);
  tree->SetBranchAddress("Vertex", &Vertex);
  tree->SetBranchAddress("threeMomPi", &threeMomPi);
  tree->SetBranchAddress("threeMomMu", &threeMomMu);
  tree->SetBranchAddress("Target", &Target);
  tree->SetBranchAddress("K3pi", &K3pi);
  tree->SetBranchAddress("autoPass", &autoPass);
  tree->SetBranchAddress("Assoc", &Assoc);
  //tree->SetBranchAddress("KTAGcand", KTAGcand);

  // Combinatorial
    
  TH1D *hDistComb = new TH1D("hDistComb", "Combinatorial background studies", 100, 0., 1000.);
  TH1D *hTimeComb = new TH1D("hTimeComb", "Combinatorial background studies", 100, -15., 15.);
  TH1D *hZComb = new TH1D("hZComb", "Combinatorial background studies", 500, 100., 190.);
  TH1D *hCDAComb = new TH1D("hCDAComb", "Combinatorial background studies", 100, 0., 0.5);
  TH1D *hInvMassComb = new TH1D("hInvMassComb", "Combinatorial background studies", 300, 0.2, 2.);
  TH1D *hInvMassCombSR = new TH1D("hInvMassCombSR", "Combinatorial background studies", 300, 0.2, 2.);
  TH1D *hMomPiComb = new TH1D("hMomPiComb", "Combinatorial background studies", 100, 0., 200.);
  TH1D *hMomMuComb = new TH1D("hMomMuComb", "Combinatorial background studies", 100, 0., 200.);
  TH2D *hDistvsMassComb = new TH2D("hDistvsMassComb", "Combinatorial background studies", 300, 0.2, 2., 50, 0., 1000.);
  TH2D *hSRComb = new TH2D("hSRComb", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalComb = new TH2D("hSRFinalComb", "Signal region", 500, -50., 50., 50, 0., 0.1);

  // Prompt

  TH1D *hDistPrompt = new TH1D("hDistPrompt", "Prompt background studies", 100, 0., 1000.);
  TH1D *hTimePrompt = new TH1D("hTimePrompt", "Prompt background studies", 100, -15., 15.);
  TH1D *hZPrompt = new TH1D("hZPrompt", "Prompt background studies", 500, 100., 190.);
  TH1D *hCDAPrompt = new TH1D("hCDAPrompt", "Prompt background studies", 100, 0., 0.5);
  TH1D *hInvMassPrompt = new TH1D("hInvMassPrompt", "Prompt background studies", 300, 0.2, 2.);
  TH1D *hInvMassPromptSR = new TH1D("hInvMassPromptSR", "Combinatorial background studies", 300, 0.2, 2.);
  TH1D *hMomPiPrompt = new TH1D("hMomPiPrompt", "Prompt background studies", 100, 0., 200.);
  TH1D *hMomMuPrompt = new TH1D("hMomMuPrompt", "Prompt background studies", 100, 0., 200.);
  TH2D *hDistvsMassPrompt = new TH2D("hDistvsMassPrompt", "Prompt background studies", 300, 0.2, 2., 50, 0., 1000.);
  TH2D *hSRPrompt = new TH2D("hSRPrompt", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalPrompt = new TH2D("hSRFinalPrompt", "Signal region", 500, -50., 50., 50, 0., 0.1);

  // Parasitic

  TH1D *hDistPar = new TH1D("hDistPar", "Parasitic background studies", 100, 0., 1000.);
  TH1D *hTimePar = new TH1D("hTimePar", "Parasitic background studies", 100, -15., 15.);
  TH1D *hZPar = new TH1D("hZPar", "Parasitic background studies", 500, 100., 190.);
  TH1D *hCDAPar = new TH1D("hCDAPar", "Parasitic background studies", 100, 0., 0.5);
  TH1D *hInvMassPar = new TH1D("hInvMassPar", "Parasitic background studies", 300, 0.2, 2.);
  TH1D *hInvMassParSR = new TH1D("hInvMassParSR", "Combinatorial background studies", 300, 0.2, 2.);
  TH1D *hMomPiPar = new TH1D("hMomPiPar", "Parasitic background studies", 100, 0., 200.);
  TH1D *hMomMuPar = new TH1D("hMomMuPar", "Parasitic background studies", 100, 0., 200.);
  TH2D *hDistvsMassPar = new TH2D("hDistvsMassPar", "Parasitic background studies", 300, 0.2, 2., 50, 0., 1000.);
  TH2D *hSRPar = new TH2D("hSRPar", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalPar = new TH2D("hSRFinalPar", "Signal region", 500, -50., 50., 50, 0., 0.1);

  // SR

  TH1D *hDistSR = new TH1D("hDistSR", "Signal region background studies", 100, 0., 1000.);
  TH1D *hTimeSR = new TH1D("hTimeSR", "Signal region background studies", 100, -15., 15.);
  TH1D *hZSR = new TH1D("hZSR", "Signal region background studies", 500, 100., 190.);
  TH1D *hCDASR = new TH1D("hCDASR", "Signal region background studies", 100, 0., 0.5);
  TH1D *hInvMassSR = new TH1D("hInvMassSR", "Signal region background studies", 300, 0.2, 2.);
  TH1D *hMomPiSR = new TH1D("hMomPiSR", "Signal region background studies", 100, 0., 200.);
  TH1D *hMomMuSR = new TH1D("hMomMuSR", "Signal region background studies", 100, 0., 200.);
  TH2D *hDistvsMassSR = new TH2D("hDistvsMassSR", "Signal region background studies", 300, 0.2, 2., 50, 0., 1000.);
  TH2D *hSRSR = new TH2D("hSRSR", "Signal region", 500, -50., 50., 50, 0., 0.1);
  
  for(Int_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    if ((ZCDALine < -10000.) || (ZCDALine > 35000.) || ((ZCDALine >= -10000. && ZCDALine <= 35000.) && CDALine > 40.)) {
      hDistSR->Fill(BeamlineDist, Weight);
      hTimeSR->Fill(CHODTime1-CHODTime2, Weight);
      hZSR->Fill(Zvertex/1000., Weight);
      hCDASR->Fill(CDALine/1000., Weight);
      hInvMassSR->Fill(invMass/1000., Weight);
      hMomPiSR->Fill(threeMomPi->Mag()/1000., Weight);
      hMomMuSR->Fill(threeMomMu->Mag()/1000., Weight);
      hDistvsMassSR->Fill(invMass/1000., BeamlineDist, Weight);
      hSRSR->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }
    
    if ((CHODTime1-CHODTime2 < -3. && CHODTime1-CHODTime2 > -5.) || (CHODTime1-CHODTime2 > 3. && CHODTime1-CHODTime2 < 5.)) {
      hDistComb->Fill(BeamlineDist, Weight);
      hTimeComb->Fill(CHODTime1-CHODTime2, Weight);
      hZComb->Fill(Zvertex/1000., Weight);
      hCDAComb->Fill(CDALine/1000., Weight);
      hInvMassComb->Fill(invMass/1000., Weight);
      hMomPiComb->Fill(threeMomPi->Mag()/1000., Weight);
      hMomMuComb->Fill(threeMomMu->Mag()/1000., Weight);
      hDistvsMassComb->Fill(invMass/1000., BeamlineDist, Weight);
      hSRComb->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }

    if (Zvertex >= 100000. && Zvertex <= 120000.) {
      hDistPrompt->Fill(BeamlineDist, Weight);
      hTimePrompt->Fill(CHODTime1-CHODTime2, Weight);
      hZPrompt->Fill(Zvertex/1000., Weight);
      hCDAPrompt->Fill(CDALine/1000., Weight);
      hInvMassPrompt->Fill(invMass/1000., Weight);
      hMomPiPrompt->Fill(threeMomPi->Mag()/1000., Weight);
      hMomMuPrompt->Fill(threeMomMu->Mag()/1000., Weight);
      hDistvsMassPrompt->Fill(invMass/1000., BeamlineDist, Weight);
      hSRPrompt->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }

    if (BeamlineDist >= 100. && BeamlineDist <= 150.) {
      hDistPar->Fill(BeamlineDist, Weight);
      hTimePar->Fill(CHODTime1-CHODTime2, Weight);
      hZPar->Fill(Zvertex/1000., Weight);
      hCDAPar->Fill(CDALine/1000., Weight);
      hInvMassPar->Fill(invMass/1000., Weight);
      hMomPiPar->Fill(threeMomPi->Mag()/1000., Weight);
      hMomMuPar->Fill(threeMomMu->Mag()/1000., Weight);
      hDistvsMassPar->Fill(invMass/1000., BeamlineDist, Weight);
      hSRPar->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }

    if (ZCDALine >= -10000. && ZCDALine <= 35000. && CDALine <= 40.) {
      if ((CHODTime1-CHODTime2 < -3. && CHODTime1-CHODTime2 > -5.) || (CHODTime1-CHODTime2 > 3. && CHODTime1-CHODTime2 < 5.)) {
	if (threeMomPi->Mag()/1000. < 70. || threeMomPi->Mag()/1000. > 80.) {
	  hSRFinalComb->Fill(ZCDALine/1000., CDALine/1000., Weight);
	  counterComb++;
	  hInvMassCombSR->Fill(invMass/1000., Weight);
	}
      }
      if (Zvertex >= 100000. && Zvertex <= 120000.) {
	if (threeMomPi->Mag()/1000. < 70. || threeMomPi->Mag()/1000. > 80.) {
	  hSRFinalPrompt->Fill(ZCDALine/1000., CDALine/1000., Weight);
	  counterPrompt++;
	  hInvMassPromptSR->Fill(invMass/1000., Weight);
	}
      }
      if (BeamlineDist >= 100. && BeamlineDist <= 150.) {
	if (threeMomPi->Mag()/1000. < 70. || threeMomPi->Mag()/1000. > 80.) {
	  hSRFinalPar->Fill(ZCDALine/1000., CDALine/1000., Weight);
	  counterPar++;
	  hInvMassParSR->Fill(invMass/1000., Weight);
	}
      }
    }
  }

  cout<<"Number of events in SR and time sidebands: "<<counterComb<<", and Z sidebands: "<<counterPrompt<<", and beamdist sidebands: "<<counterPar<<endl;

  Save(path, c, hDistSR, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hTimeSR, "Track time difference [ns]", labelSize, titleSize);
  Save(path, c, hZSR, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path, c, hCDASR, "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path, c, hInvMassSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path, c, hMomPiSR, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path, c, hMomMuSR, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path, c, hDistvsMassSR, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hSRSR, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);

  Save(path, c, hDistComb, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hTimeComb, "Track time difference [ns]", labelSize, titleSize);
  Save(path, c, hZComb, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path, c, hCDAComb, "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path, c, hInvMassComb, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path, c, hInvMassCombSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path, c, hMomPiComb, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path, c, hMomMuComb, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path, c, hDistvsMassComb, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hSRComb, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path, c, hSRFinalComb, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);

  Save(path, c, hDistPrompt, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hTimePrompt, "Track time difference [ns]", labelSize, titleSize);
  Save(path, c, hZPrompt, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path, c, hCDAPrompt, "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path, c, hInvMassPrompt, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path, c, hInvMassPromptSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path, c, hMomPiPrompt, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path, c, hMomMuPrompt, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path, c, hDistvsMassPrompt, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hSRPrompt, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path, c, hSRFinalPrompt, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);

  Save(path, c, hDistPar, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hTimePar, "Track time difference [ns]", labelSize, titleSize);
  Save(path, c, hZPar, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path, c, hCDAPar, "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path, c, hInvMassPar, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path, c, hInvMassParSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path, c, hMomPiPar, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path, c, hMomMuPar, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path, c, hDistvsMassPar, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hSRPar, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path, c, hSRFinalPar, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);

  tree->ResetBranchAddresses();
}
