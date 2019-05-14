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

void Analyzer(TString dir, TString histo1, TString an, TCanvas* c, Double_t &counterComb, Double_t &counterPar, Double_t &counterPrompt, Double_t &counterPromptPosNeg, Double_t &counterPromptPosNegFV) {
  
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
  Double_t ZEndFV = 180000.;
  Double_t xGTKMin = -30.;
  Double_t xGTKMax = 30.;
  Double_t yGTKMin = -13.5;
  Double_t yGTKMax = 13.5;
  
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
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/2016/" + analyzer + "/";
    else if (histo1.Contains("2017"))
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/2017/" + analyzer + "/";
    else if (histo1.Contains("2018"))
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/2018/" + analyzer + "/";    
    else if (histo1.Contains("Data"))
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/All/" + analyzer + "/";
    else
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/MC/" + analyzer + "/";
  }
  
  TFile *f = TFile::Open(histo1);
  
  if (f == 0) {
    cout << "Error: cannot open " << histo1 << endl;
    _exit(1);
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
  TH1D *hInvMassCombSREnriched = new TH1D("hInvMassCombSREnriched", "Enriched combinatorial background studies", 300, 0.2, 2.);
  TH1D *hMomPiComb = new TH1D("hMomPiComb", "Combinatorial background studies", 100, 0., 200.);
  TH1D *hMomMuComb = new TH1D("hMomMuComb", "Combinatorial background studies", 100, 0., 200.);
  TH1D *hGTK3IDHypComb = new TH1D("hGTK3IDHypComb", "Combinatorial background studies", 6, -0.5, 5.5);
  TH1D *hGTK3IDRingComb = new TH1D("hGTK3IDRingComb", "Combinatorial background studies", 100, 0., 250.);
  TH2D *hGTK3IDRingvsHypComb = new TH2D("hGTK3IDRingvsHypComb", "Combinatorial background studies", 6, -0.5, 5.5, 100, 0., 250.);
  TH2D *hDistvsMassComb = new TH2D("hDistvsMassComb", "Combinatorial background studies", 300, 0.2, 2., 50, 0., 1000.);
  TH2D *hGTK3XYComb = new TH2D("hGTK3XYComb", "Combinatorial background studies", 100, -50., 50., 100, -50., 50.);
  TH2D *hSRComb = new TH2D("hSRComb", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalComb = new TH2D("hSRFinalComb", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalCombEnriched = new TH2D("hSRFinalCombEnriched", "Enriched signal region", 500, -50., 50., 50, 0., 0.1);

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
    
    if (i%1000 == 0)
      cout<<"Processing event n."<<i<<endl;
    
    // 0-charge

    if ((ZCDALine < ZCDALineMin || ZCDALine > ZCDALineMax || (ZCDALine >= ZCDALineMin && ZCDALine <= ZCDALineMax && CDALine > CDALineMax)) && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // all events outside blinded region
      hDistSR->Fill(BeamlineDist, Weight);
      hTimeSR->Fill(CHODTime1-CHODTime2, Weight);
      hZSR->Fill(Zvertex/1000., Weight);
      hCDASR->Fill(CDALine/1000., Weight);
      hInvMassSR->Fill(invMass/1000., Weight);
      hMomPiSR->Fill(Mom1->Mag()/1000., Weight);
      hMomMuSR->Fill(Mom2->Mag()/1000., Weight);
      hDistvsMassSR->Fill(invMass/1000., BeamlineDist, Weight);
      hSRSR->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }
    
    if (((CHODTime1-CHODTime2 < Time1Max && CHODTime1-CHODTime2 > Time1Min) || (CHODTime1-CHODTime2 > Time2Min && CHODTime1-CHODTime2 < Time2Max)) && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // all events inside time sidebands
      hDistComb->Fill(BeamlineDist, Weight);
      hTimeComb->Fill(CHODTime1-CHODTime2, Weight);
      hZComb->Fill(Zvertex/1000., Weight);
      hCDAComb->Fill(CDALine/1000., Weight);
      hInvMassComb->Fill(invMass/1000., Weight);
      hMomPiComb->Fill(Mom1->Mag()/1000., Weight);
      hMomMuComb->Fill(Mom2->Mag()/1000., Weight);
      hDistvsMassComb->Fill(invMass/1000., BeamlineDist, Weight);
      hSRComb->Fill(ZCDALine/1000., CDALine/1000., Weight);

      // Spike at 75 GeV
      
      if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) {
	hGTK3XYComb->Fill(xGTK31, yGTK31);
	if (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax) {
	  if (RICHInfo1->at(0) == 99)
	    RICHInfo1->at(0) = 5;
	  hGTK3IDHypComb->Fill(RICHInfo1->at(0));
	  hGTK3IDRingComb->Fill(RICHInfo1->at(5));
	  hGTK3IDRingvsHypComb->Fill(RICHInfo1->at(0), RICHInfo1->at(5));
	}
      }
      if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) {
	hGTK3XYComb->Fill(xGTK32, yGTK32);
	if (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax) {
	  if (RICHInfo2->at(0) == 99)
	    RICHInfo2->at(0) = 5;
	  hGTK3IDHypComb->Fill(RICHInfo2->at(0));
	  hGTK3IDRingComb->Fill(RICHInfo2->at(5));
	  hGTK3IDRingvsHypComb->Fill(RICHInfo2->at(0), RICHInfo2->at(5));
	}
      }
    }
    
    if (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // all events inside vertex sidebands
      hDistPrompt->Fill(BeamlineDist, Weight);
      hTimePrompt->Fill(CHODTime1-CHODTime2, Weight);
      hZPrompt->Fill(Zvertex/1000., Weight);
      hCDAPrompt->Fill(CDALine/1000., Weight);
      hInvMassPrompt->Fill(invMass/1000., Weight);
      hMomPiPrompt->Fill(Mom1->Mag()/1000., Weight);
      hMomMuPrompt->Fill(Mom2->Mag()/1000., Weight);
      hDistvsMassPrompt->Fill(invMass/1000., BeamlineDist, Weight);
      hSRPrompt->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }
    
    if (BeamlineDist >= BeamdistMin && BeamlineDist <= BeamdistMax && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // all events inside beam distance sidebands
      hDistPar->Fill(BeamlineDist, Weight);
      hTimePar->Fill(CHODTime1-CHODTime2, Weight);
      hZPar->Fill(Zvertex/1000., Weight);
      hCDAPar->Fill(CDALine/1000., Weight);
      hInvMassPar->Fill(invMass/1000., Weight);
      hMomPiPar->Fill(Mom1->Mag()/1000., Weight);
      hMomMuPar->Fill(Mom2->Mag()/1000., Weight);
      hDistvsMassPar->Fill(invMass/1000., BeamlineDist, Weight);
      hSRPar->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }
    
    if (ZCDALine >= ZCDALineMin && ZCDALine <= ZCDALineMax && CDALine <= CDALineMax) { // all events inside blinded region...
      if (((CHODTime1-CHODTime2 < Time1Max && CHODTime1-CHODTime2 > Time1Min) || (CHODTime1-CHODTime2 > Time2Min && CHODTime1-CHODTime2 < Time2Max)) && CDA < 100. && !an.Contains("Pos") && !an.Contains("Neg")) { // ...and time sidebands (enriched sample)
	hSRFinalCombEnriched->Fill(ZCDALine/1000., CDALine/1000., Weight);
	hInvMassCombSREnriched->Fill(invMass/1000., Weight);
      }
      if (((CHODTime1-CHODTime2 < Time1Max && CHODTime1-CHODTime2 > Time1Min) || (CHODTime1-CHODTime2 > Time2Min && CHODTime1-CHODTime2 < Time2Max)) && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // ...and time sidebands (not enriched sample to study spike at 75 GeV)
	counterComb++;
	hSRFinalComb->Fill(ZCDALine/1000., CDALine/1000., Weight);
	hInvMassCombSR->Fill(invMass/1000., Weight);
      }
      if (BeamlineDist >= BeamdistMin && BeamlineDist <= BeamdistMax && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // ...and beam distance sidebands
	hSRFinalPar->Fill(ZCDALine/1000., CDALine/1000., Weight);
	counterPar++;
	hInvMassParSR->Fill(invMass/1000., Weight);
      }
      if (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // ...and vertex sidebands (0-charge)
	hSRFinalPrompt->Fill(ZCDALine/1000., CDALine/1000., Weight);
	counterPrompt++;
	hInvMassPromptSR->Fill(invMass/1000., Weight);	
      }
      if ((an.Contains("Pos") || an.Contains("Neg")) && Zvertex >= ZVertexMin && Zvertex <= ZVertexMax && CDA < CDAMax) { // ...and vertex sidebands (2-charge)
	counterPromptPosNeg++;
      }
      if ((an.Contains("Pos") || an.Contains("Neg")) && Zvertex >= ZVertexMax && Zvertex <= ZEndFV && CDA < CDAMax) { // ...and outside vertex sidebands (2-charge)
	counterPromptPosNegFV++;
      }
    }
  }
  
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
  Save(path, c, hInvMassCombSREnriched, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path, c, hMomPiComb, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path, c, hMomMuComb, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path, c, hGTK3IDHypComb, "RICH hypothesis", labelSize, titleSize);
  Save(path, c, hGTK3IDRingComb, "RICH radius [mm]", labelSize, titleSize);
  Save(path, c, hGTK3IDRingvsHypComb, "RICH hypothesis", "RICH radius [mm]", labelSize, titleSize);
  Save(path, c, hDistvsMassComb, "Reconstructed HNL mass [GeV/c^{2}]", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hGTK3XYComb, "X at GTK3 [m]", "Y at GTK3 [m]", labelSize, titleSize);
  Save(path, c, hSRComb, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path, c, hSRFinalComb, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path, c, hSRFinalCombEnriched, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);

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
  Double_t counterComb = 0;
  Double_t counterPar = 0;
  Double_t counterPrompt = 0;
  Double_t counterPromptFV = 0;
  Double_t counterPromptPosNeg = 0;
  Double_t counterPromptPosNegFV = 0;
  
  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.25);
  c->SetTopMargin(0.15);
  c->SetGrid();
  c->RedrawAxis();

  Analyzer(dir, histo1, "HeavyNeutrino", c, counterComb, counterPar, counterPrompt, counterPromptPosNeg, counterPromptPosNegFV);
  Analyzer(dir, histo1, "HeavyNeutrinoPos", c, counterComb, counterPar, counterPrompt, counterPromptPosNeg, counterPromptPosNegFV);
  Analyzer(dir, histo1, "HeavyNeutrinoNeg", c, counterComb, counterPar, counterPrompt, counterPromptPosNeg, counterPromptPosNegFV);

  counterPromptFV = counterPrompt*counterPromptPosNegFV/counterPromptPosNeg;
  
  cout<<"Number of events in SR and time sidebands: "<<counterComb<<", and beamdist sidebands: "<<counterPar<<", and Z sidebands (0-charge): "<<counterPrompt<<", and Z sidebands (2-charge): "<<counterPromptPosNeg<<", and FV (2-charge): "<<counterPromptPosNegFV<<", and FV (0-charge): "<<counterPromptFV<<endl;

  _exit(0);
}
