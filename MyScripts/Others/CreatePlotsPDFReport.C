void CreatePDF() {
  TFile *f1 = TFile::Open("~/NA62AnalysisTool/HNLAnalysis.root");
  TDirectory * dir1 = (TDirectory*)f1->Get("HeavyNeutrinoPiMuSelection");  
  TIter next(dir1->GetListOfKeys());
  TKey *key;
  TCanvas *c1 = new TCanvas();
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;
  Double_t x = 0.07;
  Double_t y = 0.065;

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
	gPad->Update();
	if (!strcmp(key->GetName(), "InvMassReco")) {
	  h1->Fit("gaus");
	  gStyle->SetOptFit(1);
	}
        if (!strcmp(key->GetName(), "PhysicsEventsVsCuts"))
          h1->Draw("text");
        else
          h1->Draw();
      }
      
      TString path = "/home/lorenza/GitVarious/PhD/HeavyNeutrino/Plots/";
      c1->SaveAs(path + key->GetName() + ".png");
    }
  }
}
