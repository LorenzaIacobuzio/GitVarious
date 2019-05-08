void TMultiGraphCosmetics(TMultiGraph *m, const char* x, const char* y, TCanvas* c, TString path, Double_t labelSize, Double_t titleSize) {

  m->Draw("AC");
  m->GetXaxis()->SetTitle(x);
  m->GetYaxis()->SetTitle(y);
  m->GetXaxis()->SetTitleOffset(1.4);
  m->GetYaxis()->SetTitleOffset(1.4);
  m->GetXaxis()->SetTitleSize(labelSize);
  m->GetYaxis()->SetTitleSize(labelSize);
  m->GetXaxis()->SetLabelSize(labelSize);
  m->GetYaxis()->SetLabelSize(labelSize);
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

void LaunchAllModels(TString direc, Int_t in, Int_t step, Int_t fin) {
  
  // dir = output dir, in = first value of the increasing ratio in the new Shapo models (set in = 999 if it is not a new Shapo model), step = step of the increasing ratio, fin = last value of the increasing ratio (ex. if plots are made for 1:1:0, 10:1:0 and 20:1:0, in = 1, step = 10, fin = 20)

  TCanvas *c = CreateTCanvas();
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;
  TString path = "";
  TMultiGraph *mEl = CreateTMultiGraph("Sensitivity", "Sensitivity");
  TMultiGraph *mTau = CreateTMultiGraph("Sensitivity", "Sensitivity");
  TGraph *g = new TGraph();

  if (direc != "")
    path = direc;
  else
    path = "/Users/lorenza/Desktop/Histos/";

  TSystemDirectory dir("/Users/lorenza/Desktop/FinalHistos/", "/Users/lorenza/Desktop/FinalHistos/");
  TList *files = dir.GetListOfFiles(); 

  if (files) { 
    TSystemFile *file; 
    TString fname; 
    TString prefix;
    Int_t val;
    const char * mode = "hnst";
    const char * initDir = "";

    TIter next(files); 
    while ((file = (TSystemFile*)next())) { 
      fname = file->GetName(); 
      TString name = fname;
      if (!file->IsDirectory() && fname.EndsWith(".root")) { 
	prefix = name.Remove(name.First('-'), name.Last('.') + 4);
	val = prefix.Atoi();
	if (fname.Contains(prefix + "-1-0") && !fname.Contains("52-1-1") && !fname.Contains("0.061-1-4.3") && !fname.Contains("1-16-4.3") && !fname.Contains("1-1-1") && !fname.Contains("0-1-0")) {
	  TString source = "Users/lorenza/Desktop/FinalHistos/" + fname;
	  TString a = Form(".x Plots(\"%s\",\"%s\",\"%s\",%d,%d,%d,%d)", initDir, source.Data(), mode, in, step, fin, val);
	  g = gROOT->ProcessLine(a);
	  mEl->Add(g);
	}
	else if (fname.Contains("0-1-" + prefix) && !fname.Contains("52-1-1") && !fname.Contains("0.061-1-4.3") && !fname.Contains("1-16-4.3") && !fname.Contains("1-1-1") && !fname.Contains("0-1-0")) {
	  g = gROOT->ProcessLine(a);
	  mTau->Add(g);
	}
	else {
	  TMultiGraph *m = CreateTMultiGraph("Sensitivity", "Sensitivity");
	  g = gROOT->ProcessLine(a);
	  m->Add(g);
	  m->SetTitle(g->GetTitle);
	  TMultiGraphCosmetics(m, "N mass [GeV/c^{2}]", "Log(U^{2})", c, path, labelSize, titleSize);
	}
      }
    }
  }
      
  TMultiGraphCosmetics(mEl, "N mass [GeV/c^{2}]", "Log(U^{2})", c, path, labelSize, titleSize);
  TMultiGraphCosmetics(mTau, "N mass [GeV/c^{2}]", "Log(U^{2})", c, path, labelSize, titleSize);
}
