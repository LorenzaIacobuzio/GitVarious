void TH1Cosmetics(TH1* h1, Double_t labelSize, Double_t titleSize) {

  h1->SetFillColor(38);
  h1->SetTitleSize(titleSize, "t");
  h1->GetXaxis()->SetTitleSize(labelSize);
  h1->GetXaxis()->SetLabelSize(labelSize);
  h1->GetYaxis()->SetTitleSize(labelSize);
  h1->GetYaxis()->SetLabelSize(labelSize);
  h1->GetXaxis()->SetTitleOffset(1.4);
  h1->GetYaxis()->SetTitleOffset(1.4);
  gPad->Update();
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);
  gStyle->SetOptStat(0);
  gPad->SetLogy(0);
  gPad->Update();
}

void TH2Cosmetics(TH2* h2, Bool_t logScale, Double_t labelSize, Double_t titleSize, Bool_t data) {

  h2->SetTitleSize(titleSize, "t");
  h2->GetXaxis()->SetTitleSize(labelSize);
  h2->GetYaxis()->SetTitleSize(labelSize);
  h2->GetXaxis()->SetLabelSize(labelSize);
  h2->GetYaxis()->SetLabelSize(labelSize);
  h2->GetXaxis()->SetTitleOffset(1.4);
  h2->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat(0);
  gPad->Update();
  
  if (logScale == true) {
    gPad->SetLogy();
    gStyle->SetLineStyleString(9, "80 20");
    gStyle->SetLineStyleString(9, "80 20");
    gPad->Update();
  }
  else
    gPad->SetLogy(0);
  
  if ((TString)h2->GetName() == "XYSpec0Mu" || (TString)h2->GetName() == "XYSpec0Pi") {
    h2->GetXaxis()->SetRangeUser(-1., 1.);
    h2->GetYaxis()->SetRangeUser(-1., 1.);
  }
  else if ((TString)h2->GetName() == "BeamvsTar_Fin") {
    h2->GetXaxis()->SetRangeUser(0., 0.08);
    h2->GetYaxis()->SetRangeUser(0., 0.8);
  }
  else if ((TString)h2->GetName() == "BeamvsTAX_Fin") {
    h2->GetXaxis()->SetRangeUser(0., 0.1);
    h2->GetYaxis()->SetRangeUser(0., 0.8);
  }
  /*
  if (!data) {
    if ((TString)h2->GetName() == "SignalRegionTar_Fin") {
      h2->GetXaxis()->SetRangeUser(0., 20.);
      h2->GetYaxis()->SetRangeUser(0., 0.05);
    }
    else if ((TString)h2->GetName() == "SignalRegionTAX_Fin") {
      h2->GetXaxis()->SetRangeUser(0., 50.);
      h2->GetYaxis()->SetRangeUser(0., 0.1);
    }
  }
  */
}

void TGraphCosmetics(TGraph* g, Double_t labelSize, Double_t titleSize) {

  TString title = g->GetTitle();

  if (title.Contains("width")) {
    g->GetYaxis()->SetTitle("Decay width [MeV]");
    if (title.Contains("coupling"))
      g->GetXaxis()->SetTitle("Log(U^{2})");
    else if (title.Contains("mass"))
      g->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
  }
  else if (title.Contains("lifetime")) {
    g->GetYaxis()->SetTitle("Lifetime [ns]");
    if (title.Contains("coupling"))
      g->GetXaxis()->SetTitle("Log(U^{2})");
    else if (title.Contains("mass"))
      g->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
  }
  else if (!strcmp(g->GetName(), "MassScan/MeanMass")) {
    g->SetTitle("Mean probability of N reaching and decaying in the FV vs N mass; N mass [GeV/c^{2}]; Mean probability");
  }
  else if (!strcmp(g->GetName(), "CouplingScan/MeanCoupling")) {
    g->SetTitle("Mean probability of N reaching and decaying in the FV vs coupling; Log(U^{2}); Mean probability");
  }

  gPad->Update();
  g->GetXaxis()->SetTitleOffset(1.4);
  g->GetYaxis()->SetTitleOffset(1.4);
  g->GetXaxis()->SetTitleSize(labelSize);
  g->GetYaxis()->SetTitleSize(labelSize);
  g->GetXaxis()->SetLabelSize(labelSize);
  g->GetYaxis()->SetLabelSize(labelSize);
  gStyle->SetOptStat(0);

  if (title.Contains("N POT vs N K decays"))
    gPad->SetLogy(0);
  else
    gPad->SetLogy();
  
  gPad->Update();
}

void TGraphCosmetics(TGraphAsymmErrors* g, Double_t labelSize, Double_t titleSize) {

  TString title = g->GetTitle();

  if (title.Contains("Yield"))
    g->GetYaxis()->SetTitle("Yield per POT");
  else if (title.Contains("Acceptance"))
    g->GetYaxis()->SetTitle("Acceptance");
  
  if (title.Contains("Sensitivity")) {
    g->SetPoint(10, 1.2, -5.6);
    gPad->SetLogy(0);
    gPad->Update();
    g->SetLineColor(kRed);
    g->SetLineWidth(3);
    g->SetFillColor(kRed);
    g->SetFillStyle(3001);
  }

  if (!title.Contains("Sensitivity")) {
    if (title.Contains("coupling"))
      g->GetXaxis()->SetTitle("Log(U^{2})");
    else if (title.Contains("mass"))
      g->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    else if (title.Contains("momentum")) {
      g->GetXaxis()->SetTitle("N momentum [GeV/c]");
      if (title.Contains("Acceptance") && title.Contains("momentum"))
	gPad->SetLogy();
    }
  }
  
  g->GetXaxis()->SetTitleOffset(1.4);
  g->GetYaxis()->SetTitleOffset(1.4);
  gPad->Update();
  g->GetXaxis()->SetTitleSize(labelSize);
  g->GetYaxis()->SetTitleSize(labelSize);
  g->GetXaxis()->SetLabelSize(labelSize);
  g->GetYaxis()->SetLabelSize(labelSize);
  gStyle->SetOptStat(0);
  gPad->Update();
}

void TMultiGraphCosmetics(TMultiGraph *m, const char* x, const char* y, TCanvas* c, TString path, Double_t labelSize, Double_t titleSize) {

  m->Draw("AP");
  m->GetXaxis()->SetTitle(x);
  m->GetYaxis()->SetTitle(y);
  m->GetXaxis()->SetTitleOffset(1.4);
  m->GetYaxis()->SetTitleOffset(1.4);
  m->GetXaxis()->SetTitleSize(labelSize);
  m->GetYaxis()->SetTitleSize(labelSize);
  m->GetXaxis()->SetLabelSize(labelSize);
  m->GetYaxis()->SetLabelSize(labelSize);

  TString title = m->GetTitle();

  if (title.Contains("Acceptance") && title.Contains("coupling"))
    m->GetYaxis()->SetRangeUser(1.E-6, 1.E-1);
  else if (title.Contains("Yield per POT") && title.Contains("coupling"))
    m->GetYaxis()->SetRangeUser(1.E-20, 1.E-10);
  else if (title.Contains("Acceptance") && title.Contains("mass"))
    m->GetYaxis()->SetRangeUser(1.E-7, 1.E-4);
  else if (title.Contains("Yield per POT") && title.Contains("mass"))
    m->GetYaxis()->SetRangeUser(1.E-20, 1.E-16);

  if (!title.Contains("momentum"))
    gPad->BuildLegend(0.71, 0.72, 0.98, 0.93);

  gPad->Update();
  gPad->SetLogy();
  gPad->Update();
  c->SaveAs(path + m->GetName() + ".pdf");
  c->SaveAs(path + m->GetName() + ".png");
  
  delete m;
  delete c;

  m = nullptr;
  c = nullptr;
}

TMultiGraph* CreateTMultiGraph(const char* name, const char* title) {

  TMultiGraph *m = new TMultiGraph("m", "");
  m->SetName(name);
  m->SetTitle(title);

  return m;
}

TCanvas* CreateTCanvas() {

  TCanvas *c = new TCanvas();

  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.25);
  c->SetTopMargin(0.15);
  c->SetGrid();
  c->RedrawAxis();
  
  return c;
}

void ParseDir(const char* fName, const char* dirName, TString path, TCanvas* c, TMultiGraph* m, TMultiGraph* m1, Bool_t MCcomp, Bool_t data) {

  TFile *f = TFile::Open(fName);
  TDirectory * dir = (TDirectory*)f->Get(dirName);
  TIter next(dir->GetListOfKeys());
  TKey *key;
  TH1D* hBe   = new TH1D();
  TH1D* hTa   = new TH1D();
  TH1D* hBe1  = new TH1D();
  TH1D* hTa1  = new TH1D();
  TH1D* hTemp = new TH1D();
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;
  
  if (MCcomp) {
    while ((key = (TKey*)next())) {
      if (!strcmp(dirName, "HeavyNeutrinoScan/ToyMC/DS")) {
	TFile *f1 = TFile::Open("/home/li/Desktop/NewHistos/Lorenza_D2Nmu_Npimu_comparison.root");
	TIter next1(f1->GetListOfKeys());
	TKey *key1;
	while ((key1 = (TKey*)next1())) {
	  TString title = key1->GetTitle();
	  if ((title.Contains(" N p prodinacc") && !strcmp(key1->GetName(), "hN_p_prodinacc") && !strcmp(key->GetName(), "pNS1")) || (title.Contains(" N pt prodinacc") && !strcmp(key1->GetName(), "hN_p_prodinacc") && !strcmp(key->GetName(), "ptNS1")) || (!strcmp(key1->GetName(), "hProd_p_prodinacc") && !strcmp(key->GetName(), "pmuS1")) || (!strcmp(key1->GetName(), "hProd_pt_prodinacc") && !strcmp(key->GetName(), "ptmuS1"))) {
	    TH1D* hGaia = (TH1D*)(key1->ReadObj());
	    TH1D* hMio = (TH1D*)(key->ReadObj());
	    TH1Cosmetics(hMio, labelSize, titleSize);
	    hMio->Draw();
	    TString mytitle = key->GetName();
	    if (mytitle.Contains("t"))
	      hMio->GetXaxis()->SetTitle("P_{t} [GeV/c]");
	    Float_t rightmax = hMio->GetMaximum();
	    Float_t scale = hMio->GetMaximum()/hGaia->GetMaximum();
	    //gROOT->ForceStyle();
	    hGaia->Sumw2();
	    hGaia->Scale(scale);
	    hGaia->Draw("same");
	    auto legend = new TLegend(0.71, 0.72, 0.98, 0.93);
	    legend->AddEntry(hGaia, "Toy MC (Gaia)");
	    legend->AddEntry(hMio, "Full MC (Lorenza)");
	    legend->Draw();
	    c->Update();
	    c->SaveAs(path + key->GetName() + "_comp.pdf");
	    c->SaveAs(path + key->GetName() + "_comp.png");
	  }
	}
      }
      else if (!strcmp(dirName, "HeavyNeutrinoScan/ToyMC/D0")) {
	TFile *f1 = TFile::Open("/home/li/Desktop/NewHistos/Lorenza_D2NKmu_Npimu_toys.root");
	TIter next1(f1->GetListOfKeys());
	TKey *key1;
	while ((key1 = (TKey*)next1())) {
	  TString title	= key1->GetTitle();
	  if ((title.Contains(" N p prodinacc") && !strcmp(key1->GetName(), "hN_p_prodinacc") && !strcmp(key->GetName(), "pN01")) || (title.Contains(" N pt prodinacc") && !strcmp(key1->GetName(), "hN_p_prodinacc") && !strcmp(key->GetName(), "ptN01")) || (!strcmp(key1->GetName(), "hProd_p_prodinacc") && !strcmp(key->GetName(), "pmu01")) || (!strcmp(key1->GetName(), "hProd_pt_prodinacc") && !strcmp(key->GetName(), "ptmu01"))) {
	    TH1D* hGaia = (TH1D*)(key1->ReadObj());
	    TH1D* hMio =	(TH1D*)(key->ReadObj());
	    TH1Cosmetics(hMio, labelSize, titleSize);
	    hMio->Draw();
	    TString mytitle = key->GetName();
            if (mytitle.Contains("t"))
              hMio->GetXaxis()->SetTitle("P_{t} [GeV/c]");
	    Float_t rightmax = hMio->GetMaximum();
	    Float_t scale = hMio->GetMaximum()/hGaia->GetMaximum();
	    //gROOT->ForceStyle();
	    hGaia->Sumw2();
	    hGaia->Scale(scale);
	    hGaia->Draw("same");
	    auto legend = new TLegend(0.71, 0.72, 0.98, 0.93);
	    legend->AddEntry(hGaia, "Toy MC (Gaia)");
	    legend->AddEntry(hMio, "Full MC (Lorenza)");
	    legend->Draw();
	    c->Update();
	    c->SaveAs(path + key->GetName() + "_comp.pdf");
	    c->SaveAs(path + key->GetName() + "_comp.png");
	  }
	}
      }
    }
  }
  else {
    while ((key = (TKey*)next())) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      TString Name = key->GetName();
      TString Name1 = Name.Remove(0,Name.First('/')+1);

      if (cl->InheritsFrom("TGraphAsymmErrors")) {
	TGraphAsymmErrors *g = (TGraphAsymmErrors*)(key->ReadObj());

	TGraphCosmetics(g, labelSize, titleSize);

	if (Name1.Contains("Exclusion")) {
	  g->Draw("ALF");
	}

	if (!Name.Contains("Mom")) {                                  
	  if (!Name.Contains("Yield") && !Name.Contains("Acc")) {
	    g->Draw("AL");                                         
	    c->SaveAs(path + Name1 + ".pdf");
	    c->SaveAs(path + Name1 + ".png");
	  }                                                           
	  else {                                                 
	    if (Name.Contains("Yield")) {
	      g->Draw("AP");                                
	      m->Add(g);                                             
	    }                                                     
	    else {                                            
	      g->Draw("AP");                                         
	      m1->Add(g);                                          
	    }                                                              
	  }
	}
	else {
	  TGraphCosmetics(g, labelSize, titleSize);
	  g->Draw("AP");
	  c->SaveAs(path + Name1 + ".pdf");
	  c->SaveAs(path + Name1 + ".png");
	}
      }
      else if (cl->InheritsFrom("TGraph")) {
	TGraph *g = (TGraph*)(key->ReadObj());

	TGraphCosmetics(g, labelSize, titleSize);

	if (Name1.Contains("T10") && !Name1.Contains("POT")) {
	  gStyle->SetOptFit(0);
	}

	if (Name1.Contains("T10") && !Name1.Contains("POT")) 
	  g->Draw("AP*");
	else if (Name1.Contains("POT1") || Name1.Contains("POT2") || Name1.Contains("NK"))
	  g->Draw("AP*");
	else
	  g->Draw("AC");

	c->SaveAs(path + Name1 + ".pdf");
	c->SaveAs(path + Name1 + ".png");
      }
      else if (cl->InheritsFrom("TH2")) {
	TH2 *h2 = (TH2*)key->ReadObj();

	if ((Name.Contains("Coupling") || Name.Contains("Mass")) && !Name.Contains("Energy"))
	  TH2Cosmetics(h2, true, labelSize, titleSize, data);
	else 
	  TH2Cosmetics(h2, false, labelSize, titleSize, data);

	if (data && (TString)h2->GetName() == "SignalRegionTar_Fin")
	  h2->Draw("text");
	else
	  h2->Draw("colz");
      
	c->SaveAs(path + key->GetName() + ".pdf");
	c->SaveAs(path + key->GetName() + ".png");
      }
      else if (!cl->InheritsFrom("TH2") && cl->InheritsFrom("TH1")) {
	TH1 *h1 = (TH1*)key->ReadObj();
	TH1Cosmetics(h1, labelSize, titleSize);

	if (Name1.Contains("POT") || Name1.Contains("T10"))
	  gStyle->SetOptStat(1111111111);

	if (!strcmp(key->GetName(), "InvMassReco")) {
	  h1->Fit("gaus");
	  gStyle->SetStatW(0.2);
	  gStyle->SetStatH(0.2);
	  gStyle->SetOptFit(111111);
	}
	if (!strcmp(key->GetName(), "PhysicsEventsVsCuts")) {
	  h1->Draw("text");
	  h1->SetOption("text");
	}
	if (!strcmp(key->GetName(), "ZDProd")) {
	  hBe = (TH1D*)h1->Clone("hBe");
	  hBe->SetName("ZDProdTarget");
	  hBe->SetTitle("Z of D meson production point in the target");
	  hBe->GetXaxis()->SetRangeUser(-0.25, 0.25);
	  hBe->Draw();
	  c->SaveAs(path + hBe->GetName() + ".pdf");
	  c->SaveAs(path + hBe->GetName() + ".png");
	  hTa = (TH1D*)h1->Clone("hTa");
	  hTa->SetName("ZDProdTAX");
	  hTa->SetTitle("Z of D meson production point in the TAXes");
	  hTa->GetXaxis()->SetRangeUser(23., 25.);
	  hTa->Draw();
	  c->SaveAs(path + hTa->GetName() + ".pdf");
	  c->SaveAs(path + hTa->GetName() + ".png");
	}
	if (!strcmp(key->GetName(), "ZDDecay")) {
	  hBe1 = (TH1D*)h1->Clone("hBe1");
	  hBe1->SetName("ZDDecayTarget");
	  hBe1->SetTitle("Z of D meson decay point in the target");
	  hBe1->GetXaxis()->SetRangeUser(-0.3, 0.3);
	  hBe1->Draw();
	  c->SaveAs(path + hBe1->GetName() + ".pdf");
	  c->SaveAs(path + hBe1->GetName() + ".png");
	  hTa1 = (TH1D*)h1->Clone("hTa1");
	  hTa1->SetName("ZDDecayTAX");
	  hTa1->SetTitle("Z of D meson decay point in the TAXes");
	  hTa1->GetXaxis()->SetRangeUser(23., 25.);
	  hTa1->Draw();
	  c->SaveAs(path + hTa1->GetName() + ".pdf");
	  c->SaveAs(path + hTa1->GetName() + ".png");
	}
	else
	  h1->Draw();
      
	c->SaveAs(path + key->GetName() + ".pdf");
	c->SaveAs(path + key->GetName() + ".png");
      }
    }
  }
}

void Plots(TString dir, TString histo1, Bool_t MCcomp, Bool_t data) {

  //TGaxis::SetMaxDigits(2);

  TCanvas *c = CreateTCanvas();
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;
  
  // FIRST STEP HISTOS
    
  // One value plots for selection analyzer

  TString path = "";
  
  if (dir != "")
    path = dir;
  else
    path = "/home/li/Dropbox/PhD/Talks and papers/Notes/MCnote/images/Plots/";
  
  if (histo1.Contains("1"))
    path += "1/";
  else if (histo1.Contains("2"))
    path += "2/";
  else if (histo1.Contains("3"))
    path += "3/";
  else if (histo1.Contains("Data"))
    path += "Data/";
  else if (histo1.Contains("POT"))
    path += "POT/";
  else
    path += "2/";
  
  TFile *f = TFile::Open(histo1);
  
  if (MCcomp == false) {

    if ((TDirectory*)f->Get("HeavyNeutrino") != nullptr) {
      ParseDir(histo1, "HeavyNeutrino", path, c, nullptr, nullptr, MCcomp, data);
    }

    if ((TDirectory*)f->Get("HeavyNeutrinoScan") != nullptr && data == false) {

      // One value plots for weight quantities

      //ParseDir(histo1, "HeavyNeutrinoScan/SingleValue", path, c, nullptr, nullptr, MCcomp, data);
      
      // Coupling plots
      /*
      TMultiGraph *m  = CreateTMultiGraph("YieldCoupling", "Yield per POT vs coupling");
      TMultiGraph *m1 = CreateTMultiGraph("AccCoupling",   "Acceptance vs coupling");
      
      ParseDir(histo1, "HeavyNeutrinoScan/CouplingScan", path, c, m, m1, MCcomp, data);
      
      TMultiGraphCosmetics(m,  "Log(U^{2})", "Yield per POT", c, path, labelSize, titleSize);
      c = CreateTCanvas();
      TMultiGraphCosmetics(m1, "Log(U^{2})", "Acceptance",    c, path, labelSize, titleSize);
      c = CreateTCanvas();

      // Mass plots

      m  = CreateTMultiGraph("YieldMass", "Yield per POT vs N mass");
      m1 = CreateTMultiGraph("AccMass",   "Acceptance vs N mass");
      
      ParseDir(histo1, "HeavyNeutrinoScan/MassScan", path, c, m, m1, MCcomp, data);
      
      TMultiGraphCosmetics(m,  "N mass [GeV]", "Yield per POT", c, path, labelSize, titleSize);
      c = CreateTCanvas();
      TMultiGraphCosmetics(m1, "N mass [GeV]", "Acceptance",    c, path, labelSize, titleSize);
      c = CreateTCanvas();

      // Total scan plots
      
      ParseDir(histo1, "HeavyNeutrinoScan/TotalScan", path, c, nullptr, nullptr, MCcomp, data);
      */
    }
  }
  else {

    // Toy-MC comparison plots
    
    ParseDir(histo1, "HeavyNeutrinoScan/ToyMC/DS", path, c, nullptr, nullptr, MCcomp, data);
    ParseDir(histo1, "HeavyNeutrinoScan/ToyMC/D0", path, c, nullptr, nullptr, MCcomp, data);
  }
}
