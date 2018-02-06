void CouplingPlots() {

  TFile *f1 = TFile::Open("/home/li/Desktop/HNL1GeVCoupling.root");
  TDirectory * dir1 = (TDirectory*)f1->Get("HeavyNeutrinoCouplingScan");  
  TIter next(dir1->GetListOfKeys());
  TKey *key;
  TCanvas *c1 = new TCanvas();

  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;

  TGaxis::SetMaxDigits(2);

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
      if (!strcmp(key->GetName(), "Yield"))
	g->GetYaxis()->SetRangeUser(1.E-30, 1.E-10);
      g->Draw("AC");
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
}
