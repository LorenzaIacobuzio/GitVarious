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


void Analyzer(TString dir, TString histo1, TString an, TCanvas* c, Int_t counterComb, Int_t counterPrompt, Int_t counterPar, Int_t counterParPosNeg) {

  TString path = "";
  TString analyzer = "";
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;  
  Double_t ZCDALineMax = 35000.;
  Double_t ZCDALineMin = -10000.;
  Double_t CDALineMax = 40.;
  Double_t Time1Max = -3.;
  Double_t Time1Min = -6.;
  Double_t Time2Max = 6.;
  Double_t Time2Min = 3.;
  Double_t CDAMax = 10.;
  Double_t ZVertexMax = 120000.;
  Double_t ZVertexMin = 100000.;
  Double_t BeamdistMax = 150.;
  Double_t BeamdistMin = 100.;
  
  if (an.Contains("Pos"))
    analyzer = "Pos";
  else if (an.Contains("Neg"))
    analyzer = "Neg";
  else if (!an.Contains("Pos") && !an.Contains("Neg"))
    analyzer = "Zero";
  
  if (dir != "")
    path = dir;
  else {
    if (histo1.Contains("2016"))
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/2016/" + analyzer;
    else if (histo1.Contains("2017"))
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/2017/" + analyzer;
    else if (histo1.Contains("2018"))
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/2018/" + analyzer;    
    else
      cout<<"I don't know which directory to put the plots into. Histo name is weird!"<<endl;
    _exit(1);
  }
  
  TFile *f = TFile::Open(histo1);
  
  if (f == 0) {
    cout << "Error: cannot open " << histo1 << endl;
    return;
  }

  TString treeName = an + "Passed";
  TTree* tree = (TTree*)f->Get(treeName);
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
  Double_t xGTK31;
  Double_t yGTK31;
  Double_t xGTK32;
  Double_t yGTK32;
  Double_t BeamCDA1;
  Double_t BeamCDA2;
  Double_t EtotLKr;
  TVector3 *Mom1 = new TVector3();
  TVector3 *Mom2 = new TVector3();
  TVector3 *TotMom = new TVector3();
  TVector3 *Vertex = new TVector3();
  TVector3 *Pos1 = new TVector3();
  TVector3 *Pos2 = new TVector3();
  TVector3 *BeamVtx1 = new TVector3();
  TVector3 *BeamVtx2 = new TVector3();
  Bool_t Target;
  Bool_t K3pi;
  Bool_t autoPass;
  Int_t Assoc;
  Bool_t CHANTIAssoc1;
  Bool_t CHANTIAssoc2;
  Int_t Charge1;
  Int_t Charge2;
  Int_t nSec;
  Int_t nCHOD;
  Int_t nNewCHOD;
  Int_t nLKr;
  std::vector<Double_t> *RICHInfo1 = new std::vector<Double_t>;
  std::vector<Double_t> *RICHInfo2 = new std::vector<Double_t>;
  
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
  tree->SetBranchAddress("xGTK31", &xGTK31);
  tree->SetBranchAddress("yGTK31", &yGTK31);
  tree->SetBranchAddress("xGTK32", &xGTK32);
  tree->SetBranchAddress("yGTK32", &yGTK32);
  tree->SetBranchAddress("BeamCDA1", &BeamCDA1);
  tree->SetBranchAddress("BeamCDA2", &BeamCDA2);
  tree->SetBranchAddress("EtotLKr", &EtotLKr);
  tree->SetBranchAddress("Mom1", &Mom1);
  tree->SetBranchAddress("Mom2", &Mom2);
  tree->SetBranchAddress("TotMom", &TotMom);
  tree->SetBranchAddress("Vertex", &Vertex);
  tree->SetBranchAddress("Pos1", &Pos1);
  tree->SetBranchAddress("Pos2", &Pos2);
  tree->SetBranchAddress("BeamVtx1", &BeamVtx1);
  tree->SetBranchAddress("BeamVtx2", &BeamVtx2);
  tree->SetBranchAddress("Target", &Target);
  tree->SetBranchAddress("K3pi", &K3pi);
  tree->SetBranchAddress("autoPass", &autoPass);
  tree->SetBranchAddress("Assoc", &Assoc);
  tree->SetBranchAddress("CHANTIAssoc1", &CHANTIAssoc1);
  tree->SetBranchAddress("CHANTIAssoc2", &CHANTIAssoc2);
  tree->SetBranchAddress("Charge1", &Charge1);
  tree->SetBranchAddress("Charge2", &Charge2);
  tree->SetBranchAddress("nSec", &nSec);
  tree->SetBranchAddress("nCHOD", &nCHOD);
  tree->SetBranchAddress("nNewCHOD", &nNewCHOD);
  tree->SetBranchAddress("nLKr", &nLKr);
  tree->SetBranchAddress("RICHInfo1", &RICHInfo1);
  tree->SetBranchAddress("RICHInfo2", &RICHInfo2);

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
  TH2D *hGTK3Comb = new TH2D("hGTK3Comb", "Combinatorial background studies", 100, 0., 200., 100, 0., 1.);
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
    
    if ((ZCDALine < ZCDALineMin) || (ZCDALine > ZCDALineMax) || ((ZCDALine >= ZCDALineMin && ZCDALine <= ZCDALineMax) && CDALine > CDALineMax) && CDA < CDAMax) { // all events outside blinded region
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
    
    if ((CHODTime1-CHODTime2 < Time1Max && CHODTime1-CHODTime2 > Time1Min) || (CHODTime1-CHODTime2 > Time2Min && CHODTime1-CHODTime2 < Time2Max) && CDA < CDAMax) { // all events inside time sidebands
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

    if (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax && CDA < CDAMax) { // all events inside vertex sidebands
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

    if (BeamlineDist >= BeamdistMin && BeamlineDist <= BeamdistMax && CDA < CDAMax) { // all events inside beam distance sidebands
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
    
    if (ZCDALine >= ZCDALineMin && ZCDALine <= ZCDALineMax && CDALine <= CDALineMax) { // all events inside blinded region
      if ((CHODTime1-CHODTime2 < Time1Max && CHODTime1-CHODTime2 > Time1Min) || (CHODTime1-CHODTime2 > Time2Min && CHODTime1-CHODTime2 < Time2Max) && CDA < 100.) { // and time sidebands (enriched sample)
	hSRFinalComb->Fill(ZCDALine/1000., CDALine/1000., Weight);
	counterComb++;
	hInvMassCombSR->Fill(invMass/1000., Weight);
	hGTK3Comb->Fill(Mom1->Mag()/1000., TMath::Sqrt(xGTK31*xGTK31 + yGTK31*yGTK31)/1000.);
	hGTK3Comb->Fill(Mom2->Mag()/1000., TMath::Sqrt(xGTK32*xGTK32 + yGTK32*yGTK32)/1000.);
      }
      if (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax && CDA < CDAMax) { // and vertex sidebands
	hSRFinalPrompt->Fill(ZCDALine/1000., CDALine/1000., Weight);
	counterPrompt++;
	hInvMassPromptSR->Fill(invMass/1000., Weight);
      }
      if (BeamlineDist >= BeamdistMin && BeamlineDist <= Beamdistmax && CDA < CDAMax) { // and beam distance sidebands
	hSRFinalPar->Fill(ZCDALine/1000., CDALine/1000., Weight);
	counterPar++;
	hInvMassParSR->Fill(invMass/1000., Weight);
      }
      if ((histo1.Contains("Pos") || histo1.Contains("Neg")) && CDA < 10.) { // and total charge +-2
	hSRFinalPar->Fill(ZCDALine/1000., CDALine/1000., Weight);
	counterParPosNeg++;
	hInvMassParSR->Fill(invMass/1000., Weight);	
      }
    }
  }
  
  cout<<histo1<<"; number of events in SR and time sidebands: "<<counterComb<<", and Z sidebands: "<<counterPrompt<<", and beamdist sidebands: "<<counterPar<<", and in time with charge +-2: "<<counterParPosNeg<<endl;
  
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

void TreePlots(TString dir, TString histo1) {

  // dir = output dir, histo1 = histo to do cosmetics on

  TCanvas *c = new TCanvas();  
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;  
  Int_t counterComb = 0;
  Int_t counterPrompt = 0;
  Int_t counterPar = 0;
  Int_t counterParPosNeg = 0;

  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.25);
  c->SetTopMargin(0.15);
  c->SetGrid();
  c->RedrawAxis();

  Analyzer(dir, histo1, "HeavyNeutrino", c, counterComb, counterPrompt, counterPar, counterParPosNeg);
  Analyzer(dir, histo1, "HeavyNeutrinoPos", c, counterComb, counterPrompt, counterPar, counterParPosNeg);
  Analyzer(dir, histo1, "HeavyNeutrinoNeg", c, counterComb, counterPrompt, counterPar, counterParPosNeg);
}
