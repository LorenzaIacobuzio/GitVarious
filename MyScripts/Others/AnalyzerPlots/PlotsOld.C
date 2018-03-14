void Plots() {

  // First step histos
  
  TCanvas *c1 = new TCanvas();
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;

  // One value plots for selection analyzer
  
  TFile *f1 = TFile::Open("/home/li/Desktop/bla.root");
  TDirectory * dir1 = (TDirectory*)f1->Get("HeavyNeutrino");
  TIter next(dir1->GetListOfKeys());
  TKey *key;
  TH1D* hBe  = new TH1D();
  TH1D* hTa  = new TH1D();
  TH1D* hBe1 = new TH1D();
  TH1D* hTa1 = new TH1D();
  TString path = "/home/li/Desktop/HeavyNeutrino/OneValuePlots/";

  TGaxis::SetMaxDigits(2);

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TH2")) {
      TH2 *h2 = (TH2*)key->ReadObj();
      h2->SetTitleSize(titleSize, "t");
      h2->GetXaxis()->SetTitleSize(labelSize);
      h2->GetYaxis()->SetTitleSize(labelSize);
      h2->GetXaxis()->SetLabelSize(labelSize);
      h2->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);
      gPad->Update();
      h2->Draw("colz");
    }
    else if (!cl->InheritsFrom("TH2")) {
      TH1 *h1 = (TH1*)key->ReadObj();
      h1->SetFillColor(38);
      h1->SetTitleSize(titleSize, "t");
      h1->GetXaxis()->SetTitleSize(labelSize);
      h1->GetXaxis()->SetLabelSize(labelSize);
      h1->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetStatW(0.3);
      gStyle->SetStatH(0.3);
      gStyle->SetOptStat(10);
      gPad->Update();
      if (!strcmp(key->GetName(), "InvMassReco")) {
	h1->Fit("gaus");
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);
	gStyle->SetOptFit(1);
      }
      if (!strcmp(key->GetName(), "PhysicsEventsVsCuts"))
	h1->Draw("text");
      else
	h1->Draw();
    }
    
    c1->SaveAs(path + key->GetName() + ".png");
  }

  // One value plots for weight quantities

  dir1 = (TDirectory*)f1->Get("HeavyNeutrinoScan/SingleValue");
  TIter next1(dir1->GetListOfKeys());
  
  while ((key = (TKey*)next1())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TGraph")) {
      TGraph *g = (TGraph*)key->ReadObj();
      g->GetXaxis()->SetTitleSize(labelSize);
      g->GetYaxis()->SetTitleSize(labelSize);
      g->GetXaxis()->SetLabelSize(labelSize);
      g->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);
      c1->SaveAs(path + key->GetName() + ".png");
    }
    if (!cl->InheritsFrom("TH2")) {
      TH1 *h1 = (TH1*)key->ReadObj();
      h1->SetFillColor(38);
      h1->SetTitleSize(titleSize, "t");
      h1->GetXaxis()->SetTitleSize(labelSize);
      h1->GetXaxis()->SetLabelSize(labelSize);
      h1->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetStatW(0.3);
      gStyle->SetStatH(0.3);
      gStyle->SetOptStat(10);
      gPad->Update();
      if (!strcmp(key->GetName(), "ZDProd")) {
	hBe = (TH1D*)h1->Clone("hBe");
	hBe->SetName("ZDProdTarget");
	hBe->SetTitle("Z of D meson production point in the target");
	hBe->GetXaxis()->SetRangeUser(-250., 250.);
	hBe->Draw();
	c1->SaveAs(path + hBe->GetName() + ".png");
	hTa = (TH1D*)h1->Clone("hTa");
	hTa->SetName("ZDProdTAX");
	hTa->SetTitle("Z of D meson production point in the TAXs");
	hTa->GetXaxis()->SetRangeUser(24500., 27000.);
	hTa->Draw();
	c1->SaveAs(path + hTa->GetName() + ".png");
      }
      if (!strcmp(key->GetName(), "ZDDecay")) {
	hBe1 = (TH1D*)h1->Clone("hBe1");
	hBe1->SetName("ZDDecayTarget");
	hBe1->SetTitle("Z of D meson decay point in the target");
	hBe1->GetXaxis()->SetRangeUser(-250., 300.);
	hBe1->Draw();
	c1->SaveAs(path + hBe1->GetName() + ".png");
	hTa1 = (TH1D*)h1->Clone("hTa1");
	hTa1->SetName("ZDDecayTAX");
	hTa1->SetTitle("Z of D meson decay point in the TAXs");
	hTa1->GetXaxis()->SetRangeUser(24500., 27000.);
	hTa1->Draw();
	c1->SaveAs(path + hTa1->GetName() + ".png");
      }
      c1->SaveAs(path + key->GetName() + ".png");
    }
  }

  // Coupling plots

  TMultiGraph *m = new TMultiGraph("m", "");
  m->SetName("YieldCoupling");
  m->SetTitle("Yield per POT vs coupling");
  TMultiGraph *m1 = new TMultiGraph("m1", "");
  m1->SetName("AccCoupling");
  m1->SetTitle("Acceptance vs coupling");

  dir1 = (TDirectory*)f1->Get("HeavyNeutrinoScan/CouplingScan");
  path = "~/Desktop/HeavyNeutrino/ScanPlots/Coupling/";
  TIter next2(dir1->GetListOfKeys());

  while ((key = (TKey*)next2())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TGraph")) {
      TGraph *g = (TGraph*)key->ReadObj();
      g->GetXaxis()->SetTitleSize(labelSize);
      g->GetYaxis()->SetTitleSize(labelSize);
      g->GetXaxis()->SetLabelSize(labelSize);
      g->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);

      TString Name = key->GetName();

      if (!Name.Contains("Yield") && !Name.Contains("Acc")) {
	g->Draw("AC");
	c1->SaveAs(path + key->GetName() + ".png");
      }
      else {
	if (Name.Contains("Yield")) {
	  g->GetYaxis()->SetRangeUser(1.E-30, 1.E-10);
	  g->Draw("AC");
	  m->Add(g);
	}
	else {
	  g->GetYaxis()->SetRangeUser(1.E-10, 1.E-1);
	  g->Draw("AC");
	  m1->Add(g);
	}
      }
    }
    else if (cl->InheritsFrom("TH2")) {
      TH2 *h2 = (TH2*)key->ReadObj();
      h2->SetTitleSize(titleSize, "t");
      h2->GetXaxis()->SetTitleSize(labelSize);
      h2->GetYaxis()->SetTitleSize(labelSize);
      h2->GetXaxis()->SetLabelSize(labelSize);
      h2->GetYaxis()->SetLabelSize(labelSize);
      h2->GetZaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);
      gPad->SetLogy();
      gPad->SetLogz();
      gPad->SetGridx();
      gPad->SetGridy();
      c1->SetRightMargin(0.2);
      h2->Draw("colz");
      c1->SaveAs(path + key->GetName() + ".png");
    }
  }
  m->Draw("AC");
  m->GetXaxis()->SetTitle("Log of coupling");
  m->GetYaxis()->SetTitle("Yield");
  gPad->BuildLegend(0.818, 0.223, 0.984, 0.881);
  gPad->Update();
  gPad->SetLogy();
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->Update();
  c1->SaveAs(path + m->GetName() + ".png");
  m1->Draw("AC");
  m1->GetXaxis()->SetTitle("Log of coupling");
  m1->GetYaxis()->SetTitle("Acceptance");
  gPad->BuildLegend(0.818, 0.223, 0.984, 0.881);
  gPad->Update();
  c1->SaveAs(path + m1->GetName() + ".png");

  // Mass plots

  TMultiGraph *m2 = new TMultiGraph("m2", "");
  m2->SetName("YieldMass");
  m2->SetTitle("Yield per POT vs N mass");
  TMultiGraph *m3 = new TMultiGraph("m3", "");
  m3->SetName("AccMass");
  m3->SetTitle("Acceptance vs N mass");

  dir1 = (TDirectory*)f1->Get("HeavyNeutrinoScan/MassScan");
  path = "~/Desktop/HeavyNeutrino/ScanPlots/Mass/";
  TIter next3(dir1->GetListOfKeys());

  while ((key = (TKey*)next3())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TGraph")) {
      TGraph *g = (TGraph*)key->ReadObj();
      g->GetXaxis()->SetTitleSize(labelSize);
      g->GetYaxis()->SetTitleSize(labelSize);
      g->GetXaxis()->SetLabelSize(labelSize);
      g->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);

      TString Name = key->GetName();

      if (!Name.Contains("Yield") && !Name.Contains("Acc")) {
        g->Draw("AC");
        c1->SaveAs(path + key->GetName() + ".png");
      }
      else {
        if (Name.Contains("Yield")) {
          g->GetYaxis()->SetRangeUser(1.E-30, 1.E-10);
          g->Draw("AC");
          m2->Add(g);
        }
        else {
          g->GetYaxis()->SetRangeUser(1.E-10, 1.E-1);
          g->Draw("AC");
          m3->Add(g);
        }
      }
    }
    else if (cl->InheritsFrom("TH2")) {
      TH2 *h2 = (TH2*)key->ReadObj();
      h2->SetTitleSize(titleSize, "t");
      h2->GetXaxis()->SetTitleSize(labelSize);
      h2->GetYaxis()->SetTitleSize(labelSize);
      h2->GetXaxis()->SetLabelSize(labelSize);
      h2->GetYaxis()->SetLabelSize(labelSize);
      h2->GetZaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);
      gPad->SetLogy();
      gPad->SetLogz();
      gPad->SetGridx();
      gPad->SetGridy();
      c1->SetRightMargin(0.2);
      h2->Draw("colz");
      c1->SaveAs(path + key->GetName() + ".png");
    }
  }
  m2->Draw("AC");
  m2->GetXaxis()->SetTitle("N mass [GeV]");
  m2->GetYaxis()->SetTitle("Yield");
  gPad->BuildLegend(0.818, 0.223, 0.984, 0.881);
  gPad->Update();
  gPad->SetLogy();
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->Update();
  c1->SaveAs(path + m1->GetName() + ".png");
  m3->Draw("AC");
  m3->GetXaxis()->SetTitle("N mass [GeV]");
  m3->GetYaxis()->SetTitle("Acceptance");
  gPad->BuildLegend(0.818, 0.223, 0.984, 0.881);
  gPad->Update();
  c1->SaveAs(path + m3->GetName() + ".png");

  // Total scan plots

  dir1 = (TDirectory*)f1->Get("HeavyNeutrinoScan/TotalScan");
  path = "~/Desktop/HeavyNeutrino/ScanPlots/Total/";
  TIter next4(dir1->GetListOfKeys());

  while ((key = (TKey*)next4())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TGraph")) {
      TGraph *g = (TGraph*)key->ReadObj();
      g->GetXaxis()->SetTitleSize(labelSize);
      g->GetYaxis()->SetTitleSize(labelSize);
      g->GetXaxis()->SetLabelSize(labelSize);
      g->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);
      g->Draw("AC");
      
      c1->SaveAs(path + key->GetName() + ".png");
    }
  }

  // Second step histos (--histomode)

  TFile *f2 = TFile::Open("/home/li/Desktop/minc.root");
  TDirectory * dir2 = (TDirectory*)f2->Get("HeavyNeutrinoScan/SingleValue");
  TIter next5(dir2->GetListOfKeys());

  path = "/home/li/Desktop/HeavyNeutrino/OneValuePlots/";

  while ((key = (TKey*)next5())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TGraph")) {
      TGraph *g = (TGraph*)key->ReadObj();
      g->GetXaxis()->SetTitleSize(labelSize);
      g->GetYaxis()->SetTitleSize(labelSize);
      g->GetXaxis()->SetLabelSize(labelSize);
      g->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);
      c1->SaveAs(path + key->GetName() + ".png");
    }
  }

  delete m;
  delete m1;
  delete m2;
  delete m3;

  m  = nullptr;
  m1 = nullptr;
  m2 = nullptr;
  m3 = nullptr;
  
  dir2 = (TDirectory*)f2->Get("HeavyNeutrinoScan/CouplingScan");
  path = "~/Desktop/HeavyNeutrino/ScanPlots/Coupling/";
  TIter next6(dir2->GetListOfKeys());

  while ((key = (TKey*)next6())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TGraph")) {
      TGraph *g = (TGraph*)key->ReadObj();
      g->GetXaxis()->SetTitleSize(labelSize);
      g->GetYaxis()->SetTitleSize(labelSize);
      g->GetXaxis()->SetLabelSize(labelSize);
      g->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);

      TString Name = key->GetName();
      
      if (Name.Contains("Yield")) {
	g->GetYaxis()->SetRangeUser(1.E-30, 1.E-10);
	g->Draw("AC");
	m->Add(g);
      }
      else {
	g->GetYaxis()->SetRangeUser(1.E-10, 1.E-1);
	g->Draw("AC");
	m1->Add(g);
      }
    }
  }

  m->Draw("AC");
  m->GetXaxis()->SetTitle("Log of coupling");
  m->GetYaxis()->SetTitle("Yield");
  gPad->BuildLegend(0.818, 0.223, 0.984, 0.881);
  gPad->Update();
  gPad->SetLogy();
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->Update();
  c1->SaveAs(path + m->GetName() + ".png");
  m1->Draw("AC");
  m1->GetXaxis()->SetTitle("Log of coupling");
  m1->GetYaxis()->SetTitle("Acceptance");
  gPad->BuildLegend(0.818, 0.223, 0.984, 0.881);
  gPad->Update();
  c1->SaveAs(path + m1->GetName() + ".png");
