void ScanPlots() {

  TFile *f1 = TFile::Open("/home/li/Desktop/.root");
  TDirectory * dir1 = (TDirectory*)f1->Get("HeavyNeutrinoScan/SingleValue");  
  TIter next(dir1->GetListOfKeys());
  TKey *key;
  TCanvas *c1 = new TCanvas();
  TMultiGraph *m = new TMultiGraph("m", "");
  m->SetName("m");
  m->SetTitle("Yield per POT vs coupling");
  TMultiGraph *m1 = new TMultiGraph("m1", "");
  m1->SetName("m1");
  m1->SetTitle("Yield per POT vs N mass");
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;

  TGaxis::SetMaxDigits(2);

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH2") && cl->InheritsFrom("TH1")) {
      TH1 *h1 = (TH1*)key->ReadObj();
      h1->Draw();
      h1->SetFillColor(38);
      h1->SetTitleSize(titleSize, "t");
      h1->GetXaxis()->SetTitleSize(labelSize);
      h1->GetXaxis()->SetLabelSize(labelSize);
      h1->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetStatW(0.3);
      gStyle->SetStatH(0.3);
      gStyle->SetOptStat(10);
      gPad->Update();

      TString path = "~/Desktop/HeavyNeutrino/OneValuePlots/";
      c1->SaveAs(path + key->GetName() + ".png");
    }
  }
  
  dir1 = (TDirectory*)f1->Get("HeavyNeutrinoScan/CouplingScan");  

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());    
    if (cl->InheritsFrom("TGraph")) {
      TGraph *g = (TGraph*)key->ReadObj();
      g->GetXaxis()->SetTitleSize(labelSize);
      g->GetYaxis()->SetTitleSize(labelSize);
      g->GetXaxis()->SetLabelSize(labelSize);
      g->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);
      g->Draw("AC");
      if (key->GetName().Contains("Yield")) {
	g->GetYaxis()->SetRangeUser(1.E-30, 1.E-10);
	g->Draw("AC");
	m->Add(gr);
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
    }
   
    TString path = "~/Desktop/HeavyNeutrino/ScanPlots/Coupling/";
    c1->SaveAs(path + key->GetName() + ".png");
  }

  dir1 = (TDirectory*)f1->Get("HeavyNeutrinoScan/MassScan");

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TGraph")) {
      TGraph *g = (TGraph*)key->ReadObj();
      g->GetXaxis()->SetTitleSize(labelSize);
      g->GetYaxis()->SetTitleSize(labelSize);
      g->GetXaxis()->SetLabelSize(labelSize);
      g->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);
      g->Draw("AC");
      if (key->GetName().Contains("Yield")) {
        g->GetYaxis()->SetRangeUser(1.E-30, 1.E-10);
        g->Draw("AC");
        m->Add(gr);
      }
      m->Draw("AC");
      m->GetXaxis()->SetTitle("N mass [GeV]");
      m->GetYaxis()->SetTitle("Yield");
      gPad->BuildLegend(0.818, 0.223, 0.984, 0.881);
      gPad->Update();
      gPad->SetLogy();
      gPad->SetLogz();
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->Update();
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
    }

    TString path = "~/Desktop/HeavyNeutrino/ScanPlots/Mass/";
    c1->SaveAs(path + key->GetName() + ".png");
  }

  dir1 = (TDirectory*)f1->Get("HeavyNeutrinoScan/TotalScan");

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TGraph")) {
      TGraph *g = (TGraph*)key->ReadObj();
      g->GetXaxis()->SetTitleSize(labelSize);
      g->GetYaxis()->SetTitleSize(labelSize);
      g->GetXaxis()->SetLabelSize(labelSize);
      g->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);
      g->Draw("AC");

      TString path = "~/Desktop/HeavyNeutrino/ScanPlots/Total/";
      c1->SaveAs(path + key->GetName() + ".png");
    }
  }
}
