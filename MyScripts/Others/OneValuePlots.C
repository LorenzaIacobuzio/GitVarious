void OneValuePlots() {
  TFile *f1 = TFile::Open("/home/li/Desktop/HNL1GeVAnalysis.root");
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
	if (!strcmp(key->GetName(),   "CDAvsZVertex_TotMomToBeamlineAfterGeomCuts") || 
	    !strcmp(key->GetName(),           "ZvertexvsBeamlineDistAfterGeomCuts") || 
	    !strcmp(key->GetName(), "BeamlineDistvsTargetDist_TotMomAfterGeomCuts") ||
	    !strcmp(key->GetName(),   "CDAvsZVertex_TotMomToBeamlineAfterVetoes") || 
	    !strcmp(key->GetName(),           "ZvertexvsBeamlineDistAfterVetoes") || 
	    !strcmp(key->GetName(), "BeamlineDistvsTargetDist_TotMomAfterVetoes") ||
	    !strcmp(key->GetName(),   "CDAvsZVertex_TotMomToBeamlineFinal") || 
	    !strcmp(key->GetName(),           "ZvertexvsBeamlineDistFinal") || 
	    !strcmp(key->GetName(), "BeamlineDistvsTargetDist_TotMomFinal") || 
	    !strcmp(key->GetName(),    "CDAvsZVertex_TrackToBeamlineAfterCut") || 
	    !strcmp(key->GetName(),       "CDAvsZVertex_TrackToTrackAfterCut") || 
	    !strcmp(key->GetName(),    "CDAvsZVertex_TrackToBeamlineFinal") || 
	    !strcmp(key->GetName(),       "CDAvsZVertex_TrackToTrackFinal") || 
	    !strcmp(key->GetName(),                    "SignalRegion"))
	  h2->Draw("box");
	else
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
          hh1 = (TH1D*)h1->Clone("hh1");
	  hh1->SetName("ZDDecayTarget");
	  hh1->SetTitle("Z of D meson decay point in the target");
          hh1->GetXaxis()->SetRangeUser(-250., 300.);
          hh1->Draw();
	  c1->SaveAs(path + hh1->GetName() + ".png");
	  hhh1 = (TH1D*)h1->Clone("hhh1");
	  hhh1->SetName("ZDDecayTAX");
	  hhh1->SetTitle("Z of D meson decay point in the TAXs");
          hhh1->GetXaxis()->SetRangeUser(24500., 27000.);
	  hhh1->Draw();
	  c1->SaveAs(path + hhh1->GetName() + ".png");
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
  /*  
  TFile *f2 = TFile::Open("/home/li/Desktop/HNLMCAnalysisAccYield.root");
  TDirectory * dir2 = (TDirectory*)f2->Get("HeavyNeutrino");
  TIter next2(dir2->GetListOfKeys());
  TKey *key2;
  TCanvas *c2 = new TCanvas();

  while ((key2 = (TKey*)next2())) {
    TClass *cl2 = gROOT->GetClass(key2->GetClassName());
    if (!cl2->InheritsFrom("TH1") && !cl2->InheritsFrom("TH2")) continue;
    else {
      if (cl2->InheritsFrom("TH2")) continue;
      else if (!cl2->InheritsFrom("TH2")) {
	TH1 *h12 = (TH1*)key2->ReadObj();
        if (!strcmp(key2->GetName(), "Acc") || !strcmp(key2->GetName(), "Yield")) {
          //h12->Fit("gaus");
          gStyle->SetStatW(0.2);
          gStyle->SetStatH(0.2);
          //gStyle->SetOptFit(1);
	  gStyle->SetOptStat(1110);
	  h12->SetFillColor(38);
	  h12->SetTitleSize(titleSize, "t");
	  h12->GetXaxis()->SetTitleSize(labelSize);
	  h12->GetXaxis()->SetLabelSize(labelSize);
	  h12->GetYaxis()->SetLabelSize(labelSize);
	  h12->GetXaxis()->SetRangeUser(0., 30.E-6);
	  gPad->Update();
	  h12->Draw();
	}

	TString path = "/home/li/Desktop/HeavyNeutrino/OneValuePlots/";
	c2->SaveAs(path + key2->GetName() + ".png");
      }
    }
    }*/
}
