void CreateGraphsPDFReport() {
  TFile *f1 = TFile::Open("/home/li/Desktop/Plots.root");
  TIter next(f1->GetListOfKeys());
  TKey *key;
  TCanvas *c1 = new TCanvas();

  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TCanvas")) {
      TCanvas *c = (TCanvas*)key->ReadObj();
      c->Draw();    
    }
    TString path = "/home/li/GitVarious/HeavyNeutrino/Graphs/";
    c1->SaveAs(path + key->GetName() + ".png");
  }
}
