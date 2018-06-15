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
  gPad->Update();
}

void TH2Cosmetics(TH2* h2, Bool_t logScale, Double_t labelSize, Double_t titleSize) {

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

  if ((TString)h2->GetName() == "XYSpec0Mu" || (TString)h2->GetName() == "XYSpec0Pi") {
    h2->GetXaxis()->SetRangeUser(-1., 1.);
    h2->GetYaxis()->SetRangeUser(-1., 1.);
  }
}

void TGraphCosmetics(TGraph* g, Double_t labelSize, Double_t titleSize) {

  TString title = g->GetTitle();

  if (title.Contains("width"))
    g->GetYaxis()->SetTitle("Decay width [MeV]");
  else if (title.Contains("lifetime"))
    g->GetYaxis()->SetTitle("Lifetime [ns]");
  
  if (title.Contains("Sensitivity")) {
    g->GetXaxis()->SetLimits(0.1, 100.);
    g->GetHistogram()->SetMinimum(1.E-12);
    g->GetHistogram()->SetMaximum(100.);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->Update();
  }
  
  gPad->Update();
  g->GetXaxis()->SetTitleOffset(1.4);
  g->GetYaxis()->SetTitleOffset(1.4);
  g->GetXaxis()->SetTitleSize(labelSize);
  g->GetYaxis()->SetTitleSize(labelSize);
  g->GetXaxis()->SetLabelSize(labelSize);
  g->GetYaxis()->SetLabelSize(labelSize);
  gStyle->SetOptStat(0);
  gPad->Update();
}

void TGraphCosmetics(TGraphAsymmErrors* g, Double_t labelSize, Double_t titleSize) {

  TString title = g->GetTitle();

  if (title.Contains("Yield"))
    g->GetYaxis()->SetTitle("Yield per POT");
  else if (title.Contains("Acceptance"))
    g->GetYaxis()->SetTitle("Acceptance");
  
  if (!title.Contains("Sensitivity")) {
    if (title.Contains("coupling"))
      g->GetXaxis()->SetTitle("Log(U^{2})");
    else if (title.Contains("mass"))
      g->GetXaxis()->SetTitle("N mass [GeV/c^{2}]");
    else if (title.Contains("momentum"))
      g->GetXaxis()->SetTitle("N momentum [GeV/c^{2}]");
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
    m->GetYaxis()->SetRangeUser(1.E-12, 1.);
  else if (title.Contains("Yield per POT") && title.Contains("coupling"))
    m->GetYaxis()->SetRangeUser(1.E-28, 1.E-10);
  else if (title.Contains("Acceptance") && title.Contains("mass"))
    m->GetYaxis()->SetRangeUser(1.E-6, 1.E-4);
  else if (title.Contains("Yield per POT") && title.Contains("mass"))
    m->GetYaxis()->SetRangeUser(1.E-19, 1.E-16);

  if (!title.Contains("momentum"))
    gPad->BuildLegend(0.71, 0.72, 0.98, 0.93);

  gPad->Update();
  gPad->SetLogy();
  gPad->Update();
  c->SaveAs(path + m->GetName() + ".pdf");

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

void ParseDir(const char* fName, const char* dirName, TString path, TCanvas* c, TMultiGraph* m, TMultiGraph* m1) {

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

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    TString Name = key->GetName();
    TString Name1 = Name.Remove(0,Name.First('/')+1);

    if (cl->InheritsFrom("TGraphAsymmErrors")) {
      TGraphAsymmErrors *g = (TGraphAsymmErrors*)(key->ReadObj());

      TGraphCosmetics(g, labelSize, titleSize);

      if (!Name.Contains("Mom")) {                                  
        if (!Name.Contains("Yield") && !Name.Contains("Acc")) {
          g->Draw("AL");                                         
          c->SaveAs(path + key->GetName() + ".pdf");                 
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
	g->Draw("AL");
	c->SaveAs(path + "Error" + key->GetName() + ".pdf");
      }
    }
    else if (cl->InheritsFrom("TGraph")) {
       TGraph *g = (TGraph*)(key->ReadObj());

       TGraphCosmetics(g, labelSize, titleSize);

       if (Name1.Contains("Exclusion"))
	 g->Draw("AL");
       else
	 g->Draw("AC");

       c->SaveAs(path + Name1 + ".pdf");
    }
    else if (cl->InheritsFrom("TH2")) {
      TH2 *h2 = (TH2*)key->ReadObj();

      if (Name.Contains("Coupling") || Name.Contains("Mass"))
	TH2Cosmetics(h2, true, labelSize, titleSize);
      else
	TH2Cosmetics(h2, false, labelSize, titleSize);
      
      h2->Draw("colz");
      c->SaveAs(path + key->GetName() + ".pdf");
    }
    else if (!cl->InheritsFrom("TH2") && cl->InheritsFrom("TH1")) {
      TH1 *h1 = (TH1*)key->ReadObj();
      TH1Cosmetics(h1, labelSize, titleSize);

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
        hBe->GetXaxis()->SetRangeUser(-250., 250.);
        hBe->Draw();
	c->SaveAs(path + hBe->GetName() + ".pdf");
        hTa = (TH1D*)h1->Clone("hTa");
        hTa->SetName("ZDProdTAX");
        hTa->SetTitle("Z of D meson production point in the TAXs");
        hTa->GetXaxis()->SetRangeUser(23000., 25000.);
        hTa->Draw();
        c->SaveAs(path + hTa->GetName() + ".pdf");
      }
      if (!strcmp(key->GetName(), "ZDDecay")) {
	hBe1 = (TH1D*)h1->Clone("hBe1");
        hBe1->SetName("ZDDecayTarget");
        hBe1->SetTitle("Z of D meson decay point in the target");
        hBe1->GetXaxis()->SetRangeUser(-250., 300.);
        hBe1->Draw();
	c->SaveAs(path + hBe1->GetName() + ".pdf");
        hTa1 = (TH1D*)h1->Clone("hTa1");
	hTa1->SetName("ZDDecayTAX");
        hTa1->SetTitle("Z of D meson decay point in the TAXs");
	hTa1->GetXaxis()->SetRangeUser(23000., 25000.);
        hTa1->Draw();
        c->SaveAs(path + hTa1->GetName() + ".pdf");
      }
      else
        h1->Draw();
      
      c->SaveAs(path + key->GetName() + ".pdf");
    }
  }
}

void Plots(TString dir, TString histo1, Bool_t MCcomp) {

  TGaxis::SetMaxDigits(2);
  
  TCanvas *c = CreateTCanvas();
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;
  
  // FIRST STEP HISTOS
    
  // One value plots for selection analyzer

  TString path = "";
  
  if (dir != "")
    path = dir;
  else
    path = "/home/li/Dropbox/PhD/Talks, papers, posters, notes/Notes/MCnote/images/Plots/";
  
  if (histo1.Contains("1"))
    path += "1/";
  else if (histo1.Contains("2"))
    path += "2/";
  else if (histo1.Contains("3"))
    path += "3/";
  else {
    path += "2/";
  }

  if (MCcomp == false) {
    ParseDir(histo1, "HeavyNeutrino", path, c, nullptr, nullptr);
    
    // One value plots for weight quantities
    
    ParseDir(histo1, "HeavyNeutrinoScan/SingleValue", path, c, nullptr, nullptr);
    
    // Coupling plots
    
    TMultiGraph *m  = CreateTMultiGraph("YieldCoupling", "Yield per POT vs coupling");
    TMultiGraph *m1 = CreateTMultiGraph("AccCoupling",   "Acceptance vs coupling");
    
    ParseDir(histo1, "HeavyNeutrinoScan/CouplingScan", path, c, m, m1);
    
    TMultiGraphCosmetics(m,  "Log(U^{2})", "Yield per POT", c, path, labelSize, titleSize);
    c = CreateTCanvas();
    TMultiGraphCosmetics(m1, "Log(U^{2})", "Acceptance",    c, path, labelSize, titleSize);
    c = CreateTCanvas();
    
    // Mass plots
    
    m  = CreateTMultiGraph("YieldMass", "Yield per POT vs N mass");
    m1 = CreateTMultiGraph("AccMass",   "Acceptance vs N mass");
    
    ParseDir(histo1, "HeavyNeutrinoScan/MassScan", path, c, m, m1);
    
    TMultiGraphCosmetics(m,  "N mass [GeV]", "Yield per POT", c, path, labelSize, titleSize);
    c = CreateTCanvas();
    TMultiGraphCosmetics(m1, "N mass [GeV]", "Acceptance",    c, path, labelSize, titleSize);
    c = CreateTCanvas();
    
    // Total scan plots
    
    ParseDir(histo1, "HeavyNeutrinoScan/TotalScan", path, c, nullptr, nullptr);
  }
  else {
    // Toy-MC comparison plots
    
    ParseDir(histo1, "HeavyNeutrinoScan/ToyMC/DS", path, c, nullptr, nullptr);
    ParseDir(histo1, "HeavyNeutrinoScan/ToyMC/D0", path, c, nullptr, nullptr);
  }
}
