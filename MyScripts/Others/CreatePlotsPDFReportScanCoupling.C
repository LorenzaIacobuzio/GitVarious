void CreatePlotsPDFReportScanCoupling() {

  TFile *f1 = TFile::Open("~/NA62AnalysisTool/bla.root");
  TDirectory * dir1 = (TDirectory*)f1->Get("HeavyNeutrinoPiMuSelectionScanCoupling");  
  TIter next(dir1->GetListOfKeys());
  TKey *key;
  TCanvas *c1 = new TCanvas();

  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TGraph")) {
      TGraph *g = (TGraph*)key->ReadObj();
      //g->SetTitleSize(titleSize, "t");
      g->GetXaxis()->SetTitleSize(labelSize);
      g->GetYaxis()->SetTitleSize(labelSize);
      g->GetXaxis()->SetLabelSize(labelSize);
      g->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);
      g->Draw();
    }
    else if (cl->InheritsFrom("TH2")) {
      TH2 *h2 = (TH2*)key->ReadObj();
      h2->SetTitleSize(titleSize, "t");
      h2->GetXaxis()->SetTitleSize(labelSize);
      h2->GetYaxis()->SetTitleSize(labelSize);
      h2->GetXaxis()->SetLabelSize(labelSize);
      h2->GetYaxis()->SetLabelSize(labelSize);
      gStyle->SetOptStat(0);
      h2->Draw("colz");    
      gPad->SetLogy();
      gPad->Update();
    }
   
    TString path = "~/NA62AnalysisTool/MyScripts/Others/";
    c1->SaveAs(path + key->GetName() + ".png");
  }
}
