void OneValuePlots() {
  TFile *f1 = TFile::Open("/home/li/Desktop/.root");
  TDirectory * dir1 = (TDirectory*)f1->Get("HeavyNeutrino");  
  TIter next(dir1->GetListOfKeys());
  TKey *key;
  TCanvas *c1 = new TCanvas();
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;
  TH1D* hh1 = new TH1D();
  TH1D* hhh1 = new TH1D();      
  TString path = "/home/li/Desktop/HeavyNeutrino/OneValuePlots/";
  
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1") && !cl->InheritsFrom("TH2")) continue;
    else {
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
	if (!strcmp(key->GetName(), "ZDProd")) {
	  hh1 = (TH1D*)h1->Clone("hh1");
	  hh1->SetName("ZDProdTarget");
	  hh1->SetTitle("Z of D meson production point in the target");
	  hh1->GetXaxis()->SetRangeUser(-250., 250.);
	  hh1->Draw();
	  c1->SaveAs(path + hh1->GetName() + ".png");
	  hhh1 = (TH1D*)h1->Clone("hhh1");
	  hhh1->SetName("ZDProdTAX");
	  hhh1->SetTitle("Z of D meson production point in the TAXs");
	  hhh1->GetXaxis()->SetRangeUser(24500., 27000.);
	  hhh1->Draw();
	  c1->SaveAs(path + hhh1->GetName() + ".png");
	}
	if (!strcmp(key->GetName(), "ZDDecay")) {
          hhhh1 = (TH1D*)h1->Clone("hhhh1");
	  hhhh1->SetName("ZDDecayTarget");
	  hhhh1->SetTitle("Z of D meson decay point in the target");
          hhhh1->GetXaxis()->SetRangeUser(-250., 300.);
          hhhh1->Draw();
	  c1->SaveAs(path + hhhh1->GetName() + ".png");
	  hhhhh1 = (TH1D*)h1->Clone("hhhhh1");
	  hhhhh1->SetName("ZDDecayTAX");
	  hhhhh1->SetTitle("Z of D meson decay point in the TAXs");
          hhhhh1->GetXaxis()->SetRangeUser(24500., 27000.);
	  hhhhh1->Draw();
	  c1->SaveAs(path + hhhhh1->GetName() + ".png");
        }
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
  }
}
