void TH1Cosmetics(TH1* h1, Double_t labelSize, Double_t titleSize) {

  TString name = h1->GetName();
  TString title = h1->GetTitle();
  
  if (title.Contains("momentum") && (title.Contains("acceptance") || title.Contains("Yield"))) {
    h1->GetXaxis()->SetTitle("N momentum [GeV/c]");
    //if (title.Contains("Acceptance") && title.Contains("momentum"))
    gPad->SetLogy();
    h1->GetXaxis()->SetTitleOffset(1.4);
    h1->GetYaxis()->SetTitleOffset(1.4);
    gPad->Update();
    h1->GetXaxis()->SetTitleSize(labelSize);
    h1->GetYaxis()->SetTitleSize(labelSize);
    h1->GetXaxis()->SetLabelSize(labelSize);
    h1->GetYaxis()->SetLabelSize(labelSize);
    h1->SetMarkerSize(0.8);
    gStyle->SetOptStat(0);
    gPad->Update();
  }
  else {
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
    
    if (name.Contains("EoP"))
      gPad->SetLogy(1);
    else
      gPad->SetLogy(0);
    
    gPad->Update();
  }
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

  TString name = h2->GetName();
  
  if (name.Contains("EoP") || name.Contains("CDAvsZCDA") || name.Contains("SR"))
    gPad->SetLogz();
  else
    gPad->SetLogz(0);
  
  if (logScale == true && !name.Contains("Prob") && !name.Contains("WeightMass") && !name.Contains("WeightCoupling")) {
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
  
  if (name.Contains("Res2D"))
    gPad->SetLogy(0);

  if (name.Contains("Prob")) {
    gPad->SetLogz();
    gPad->SetLogy(0);
  }
}

void TGraphCosmetics(TGraph* g, Double_t labelSize, Double_t titleSize) {

  TString title = g->GetTitle();
  TString name = g->GetName();

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

  if (title.Contains("Contour"))
    gPad->SetLogy(0);
  else if (name.Contains("CL") || name.Contains("MeanMass") || name.Contains("MeanCoupling"))
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

  if (!title.Contains("Sensitivity") && !title.Contains("momentum")) {
    if (!title.Contains("Target") && !title.Contains("TAX"))
      g->SetMarkerStyle(24);
  }
      
  g->GetXaxis()->SetTitleOffset(1.4);
  g->GetYaxis()->SetTitleOffset(1.4);
  gPad->Update();
  g->GetXaxis()->SetTitleSize(labelSize);
  g->GetYaxis()->SetTitleSize(labelSize);
  g->GetXaxis()->SetLabelSize(labelSize);
  g->GetYaxis()->SetLabelSize(labelSize);
  g->SetMarkerSize(0.8);
  gStyle->SetOptStat(0);
  gPad->Update();
}

void TMultiGraphCosmetics(TMultiGraph *m, const char* x, const char* y, TCanvas* c, TString path, Double_t labelSize, Double_t titleSize) {

  TString name = m->GetName();

  if (!name.Contains("El") && !name.Contains("Tau"))
    m->Draw("AP");
  else 
    m->Draw("AL");

  if (name.Contains("Coupling") && name.Contains("Sel"))
    m->GetYaxis()->SetRangeUser(1.E-2, 1.);

  if (name.Contains("Coupling") && (name.Contains("Yield") || name.Contains("FV"))) {
    m->GetYaxis()->SetRangeUser(1.E-34, 1.E-6);
  }

  if (name.Contains("Yield") && name.Contains("Mass")) {
    m->GetYaxis()->SetRangeUser(1.E-20, 1.E-15);
  }

  if (name.Contains("Mass") && name.Contains("FV")) {
    m->GetYaxis()->SetRangeUser(1.E-16, 1.E-11);
  }

  if (name.Contains("Mass") && name.Contains("Reg")) {
    m->GetYaxis()->SetRangeUser(1.E-4, 1.E-3);
  }

  gPad->Update();
  m->GetXaxis()->SetTitle(x);
  m->GetYaxis()->SetTitle(y);
  m->GetXaxis()->SetTitleOffset(1.4);
  m->GetYaxis()->SetTitleOffset(1.4);
  m->GetXaxis()->SetTitleSize(labelSize);
  m->GetYaxis()->SetTitleSize(labelSize);
  m->GetXaxis()->SetLabelSize(labelSize);
  m->GetYaxis()->SetLabelSize(labelSize);

  if (!name.Contains("MergedContours")) {
    gPad->BuildLegend(0.77, 0.77, 0.98, 0.93);
    gPad->Update();
    if (!name.Contains("El") && !name.Contains("Tau")) {
      gPad->SetLogy();
      gPad->Update();
    }
  }

  gPad->Update();

  if (!name.Contains("El") && !name.Contains("Tau")) {
    c->SaveAs(path + m->GetName() + ".pdf");
    c->SaveAs(path + m->GetName() + ".png");
  }
  else {
    c->SaveAs("/home/li/Desktop/Histos/" + name + ".pdf");
    c->SaveAs("/home/li/Desktop/Histos/" + name + ".png");
  }

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

TGraph* ParseDir(const char* fName, const char* dirName, TString path, TCanvas* c, TMultiGraph* m, TMultiGraph* m1, TMultiGraph* m2, TMultiGraph* m3) {

  TFile *f = TFile::Open(fName);
  TString histo1 = (TString)(fName);
  TString histo1Name = (TString)(fName);
  TDirectory * dir = (TDirectory*)f->Get(dirName);
  TIter next(dir->GetListOfKeys());
  TKey *key;
  TH1D* hBe   = new TH1D();
  TH1D* hTa   = new TH1D();
  TH1D* hBe1  = new TH1D();
  TH1D* hTa1  = new TH1D();
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;
  TGraph *retG = new TGraph();

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    TString Name = key->GetName();
    TString Name1 = Name.Remove(0, Name.First('/') + 1);

    if (cl->InheritsFrom("TGraphAsymmErrors")) {
      TGraphAsymmErrors *g = (TGraphAsymmErrors*)(key->ReadObj());
      
      TGraphCosmetics(g, labelSize, titleSize);

      if (Name1.Contains("Exclusion")) {
	g->Draw("AP*");
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
	    m->Draw("AP");
	  }                                                     
	  else if (Name.Contains("Acc")) {
	    if (Name.Contains("Sel")) {
	      g->Draw("AP");                                         
	      m1->Add(g); 
	      m1->Draw("AP");            
	    }
	    else if (Name.Contains("Reg")) {                             
	      g->Draw("AP");
	      m2->Add(g);
	      m2->Draw("AP");
	    }
	    else if (Name.Contains("FV")) {      
	      g->Draw("AP");
	      m3->Add(g);
	      m3->Draw("AP");
	    }
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
      /*
      if (Name1.Contains("T10") && !Name1.Contains("POT")) {
	gStyle->SetOptFit(0);
      }
      
      if (Name1.Contains("T10") && !Name1.Contains("POT")) 
	g->Draw("AP*");
      else if (Name1.Contains("POT1") || Name1.Contains("POT2") || Name1.Contains("NK"))
	g->Draw("AP*");
      else if (Name1.Contains("CL")) {
	g->Draw("AC");
	m->Add(g);
      }
      */
      if (Name1.Contains("Contours")) {
	g->Draw("AL");
	retG = g;
	Int_t N = 0;
	if (histo1.Contains("Users"))
	  N = 35;
	else if (histo1.Contains("li"))
	  N = 29;
	TString modName = histo1Name.Remove(0, histo1Name.First('/') + N);
	modName.Remove(modName.Last('.'), modName.Last('.') + 4);
	histo1Name = histo1;
	c->SaveAs("/home/li/Desktop/Histos/" + modName + ".png");
	g->SetTitle(modName);
      }
      else
	g->Draw("AC");

      c->SaveAs(path + Name1 + ".pdf");
      c->SaveAs(path + Name1 + ".png");
    }
    else if (cl->InheritsFrom("TH2")) {
      TH2 *h2 = (TH2*)key->ReadObj();

      if (Name1.Contains("DThetaMom"))
	h2->RebinX(4);

      if (Name.Contains("Exclusion"))
	h2->Draw("cont z");
      else
	h2->Draw("colz");
      
      if ((Name.Contains("Coupling") || Name.Contains("Mass")) && !Name.Contains("Energy"))
	TH2Cosmetics(h2, true, labelSize, titleSize);
      else 
	TH2Cosmetics(h2, false, labelSize, titleSize);

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
      else {
	h1->Draw();
      }
      
      c->SaveAs(path + key->GetName() + ".pdf");
      c->SaveAs(path + key->GetName() + ".png");
    }
  }

  return retG;
}

TGraph* SingleHisto(TString dir, TString histo1, TString mode, Int_t val, Int_t suf, TCanvas *c) {

  // dir = output dir, histo1 = histo to do cosmetics on, mode: all = all dirs, hn = HeavyNeutrino dir, hnss = HeavyNeutrinoScan/SingleValue, hnsc = HeavyNeutrinoScan/Coupling, hnsm = HeavyNeutrinoScan/Mass, hnst = HeavyNeutrinoScan/Total

  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;
  TString path = "";
  TGraph* retG = new TGraph();

  if (dir != "")
    path = dir;
  else
    path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/";
  
  if (histo1.Contains("52"))
    path += "ShapoOld/1/";
  else if (histo1.Contains("16"))
    path += "ShapoOld/2/";
  else if (histo1.Contains("0.061"))
    path += "ShapoOld/3/";
  else if (histo1.Contains("1-1-1"))
    path += "ShapoNew/All/";
  else if (histo1.Contains("0-1-0") && val == 0)
    path += "ShapoNew/Muon/";
  else if (histo1.Contains("-1-0") && val != 0) {
    TString prefix = Form("%d", val);
    path += "ShapoNew/El/" + prefix + "/";
  }
  else if (histo1.Contains("0-1-") && suf != 0) {
    TString prefix = Form("%d", suf);
    path += "ShapoNew/Tau/" + prefix + "/";
  }
  else {
    cout<<"Don't know which directory to put " << histo1 << " into!"<<endl;
    _exit(1);
  }
  
  TFile *f = TFile::Open(histo1);
  TString hnDir = "HeavyNeutrino";
  TString hnsDir = "HeavyNeutrinoScan";
  
  if (path.Contains("El")) {
    hnDir += "NewEl";
    hnsDir += "NewEl";
  }
  else if (path.Contains("Tau")) {
    hnDir += "NewTau";
    hnsDir += "NewTau";
  }
  
  if ((TDirectory*)f->Get(hnDir) != nullptr && (mode == "all" || mode == "hn")) {
    
    ParseDir(histo1, hnDir, path + hnDir + "/", c, nullptr, nullptr, nullptr, nullptr);
    //ParseDir(histo1, hnDir + "/" + hnDir, path + hnDir + "/", c, nullptr, nullptr, nullptr, nullptr); //old format where hn dir is in hn dir (to be dismissed with next grid production
  }

  if (f->Get(hnsDir) != nullptr) {

    TMultiGraph *m  = CreateTMultiGraph("YieldCoupling", "Yield per POT vs coupling");
    TMultiGraph *m1 = CreateTMultiGraph("AccSelCoupling", "Selection acceptance vs coupling");
    TMultiGraph *m2 = CreateTMultiGraph("AccRegCoupling", "Regeneration acceptance vs coupling");
    TMultiGraph *m3 = CreateTMultiGraph("AccFVCoupling", "FV acceptance vs coupling");
    
    if (mode == "all" || mode == "hnss") {
      
      // One value plots for weight quantities
      
      ParseDir(histo1, hnsDir + "/SingleValue", path + hnsDir + "/SingleValue/", c, nullptr, nullptr, nullptr, nullptr);
    }

    if (mode == "all" || mode == "hnsc") {
      
      // Coupling plots
      
      ParseDir(histo1, hnsDir + "/CouplingScan", path + hnsDir + "/CouplingScan/", c, m, m1, m2, m3);

      c = CreateTCanvas();      
      TMultiGraphCosmetics(m, "Log(U^{2})", "Yield per POT", c, path + hnsDir + "/CouplingScan/", labelSize, titleSize);
      c = CreateTCanvas();
      TMultiGraphCosmetics(m1, "Log(U^{2})", "Selection acceptance", c, path + hnsDir + "/CouplingScan/", labelSize, titleSize);
      c = CreateTCanvas();
      TMultiGraphCosmetics(m2, "Log(U^{2})", "Regeneration acceptance", c, path + hnsDir + "/CouplingScan/", labelSize, titleSize);
      c = CreateTCanvas();
      TMultiGraphCosmetics(m3, "Log(U^{2})", "FV acceptance", c, path + hnsDir + "/CouplingScan/", labelSize, titleSize);
      c = CreateTCanvas();
    }

    if (mode == "all" || mode == "hnsm") {
      
      // Mass plots

      c = CreateTCanvas();      
      m = CreateTMultiGraph("YieldMass", "Yield per POT vs N mass");
      m1 = CreateTMultiGraph("AccSelMass", "Selection acceptance vs N mass");
      m2 = CreateTMultiGraph("AccRegMass", "Regeneration acceptance vs N mass");
      m3 = CreateTMultiGraph("AccFVMass", "FV acceptance vs N mass");
      
      ParseDir(histo1, hnsDir + "/MassScan", path + hnsDir + "/MassScan/", c, m, m1, m2, m3);
      
      TMultiGraphCosmetics(m, "N mass [GeV/c^{2}]", "Yield per POT", c, path + hnsDir + "/MassScan/", labelSize, titleSize);
      c = CreateTCanvas();      
      TMultiGraphCosmetics(m1, "N mass [GeV/c^{2}]", "Selection acceptance", c, path + hnsDir + "/MassScan/", labelSize, titleSize);
      c = CreateTCanvas();      
      TMultiGraphCosmetics(m2, "N mass [GeV/c^{2}]", "Regeneration acceptance", c, path + hnsDir + "/MassScan/", labelSize, titleSize);
      c = CreateTCanvas();      
      TMultiGraphCosmetics(m3, "N mass [GeV/c^{2}]", "FV acceptance", c, path + hnsDir + "/MassScan/", labelSize, titleSize);
      c = CreateTCanvas();
    }

    if (mode == "all" || mode == "hnst") {
      
      // Total scan plots

      retG = nullptr;
      retG = ParseDir(histo1, hnsDir + "/TotalScan", path + hnsDir + "/TotalScan/", c, nullptr,nullptr, nullptr, nullptr);
    }
  }

  return retG;
}

void MultiHisto(TString direc, TString mode) {

  TMultiGraph *mEl = CreateTMultiGraph("El", "Electron model");
  TMultiGraph *mTau = CreateTMultiGraph("Tau", "Tau model");
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;
  TString path;
  Bool_t multi = false;
  Int_t colors[14] = {602, 434, 887, 861, 623, 632, 797, 803, 402, 419, 416, 1, 922, 886}; //blue+2, cyan+2, violet+7, azure+1, red-9, red, orange-3, orange+3, yellow+2, green+3, green, black, grey+2, violet+6
  Int_t colC = 0;

  if (direc != "")
    path = direc;
  else
    path = "/home/li/Desktop/Histos/";

  TSystemDirectory dir("/home/li/Desktop/FinalHistos/", "/home/li/Desktop/FinalHistos/");
  TList *files = dir.GetListOfFiles();
  TCanvas *c = CreateTCanvas();

  if (files) {
    TSystemFile *file;
    TString fname;
    TString prefix;
    TString suffix;
    Int_t val;
    Int_t suf;
    TGraph *g = new TGraph();
    TIter next(files);

    while ((file = (TSystemFile*)next())) {
      fname = file->GetName();
      TString name = fname;
      if (!file->IsDirectory() && fname.EndsWith(".root")) {
        prefix = name.Remove(name.First('-'), name.Last('.') + 4);
        val = prefix.Atoi();
	name = fname;
	suffix = name.Remove(name.Last('.'), name.Last('.') + 4);
	suffix.Remove(0, suffix.Last('-') + 1);
	suf = suffix.Atoi();
	name = fname;
	cout << "****************** I am processing file " << fname << endl;
	if (mode.Contains("all") || mode.Contains("hnst")) {
	  if (!fname.Contains("52-1-1") && !fname.Contains("0.061-1-4.3") && !fname.Contains("1-16-3.8") && !fname.Contains("1-1-1")) {
	    if (fname.Contains("-1-0") && val != 0) {
	      multi = true;
	      g = SingleHisto("", "/home/li/Desktop/FinalHistos/" + fname, mode, val, suf, c);
	      g->SetLineColor(colors[colC]);
	      g->Draw("AL");
	      mEl->Add(g);
	      mEl->Draw("AL");
	      colC++;
	    }
	    else if (fname.Contains("0-1-") && suf != 0) {
	      multi = true;
	      g = SingleHisto("", "/home/li/Desktop/FinalHistos/" + fname, mode, val, suf, c);
	      g->SetLineColor(colors[colC]);
	      g->Draw("AL");
	      mTau->Add(g);
	      mTau->Draw("AL");
	      colC++;
	    }
	    else {
	      SingleHisto("", "/home/li/Desktop/FinalHistos/" + fname, mode, val, suf, c);
	    }
	  }
	  else {
	    SingleHisto("", "/home/li/Desktop/FinalHistos/" + fname, mode, val, suf, c);
	  }
	}
	else {
	  SingleHisto("", "/home/li/Desktop/FinalHistos/" + fname, mode, val, suf, c);
	}
      }
    }
  }
  
  if (multi) {
    c = CreateTCanvas();
    TMultiGraphCosmetics(mEl, "N mass [GeV/c^{2}]", "Log(U^{2})", c, path, labelSize, titleSize);
    c = CreateTCanvas();
    TMultiGraphCosmetics(mTau, "N mass [GeV/c^{2}]", "Log(U^{2})", c, path, labelSize, titleSize);
    c = CreateTCanvas();
  }
}

void Plots(TString mode) {

  MultiHisto("", mode);
  exit(0);
}
