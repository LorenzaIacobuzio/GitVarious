void CreateGraphsPDFReport() {
  TFile *f1 = TFile::Open("/home/li/Desktop/Plots.root");
  TString path = "/home/li/GitVarious/HeavyNeutrino/Graphs/";
  TGaxis::SetMaxDigits(2);
  TCanvas *c = (TCanvas *)f1->Get("c");
  c->Draw();
  c->SaveAs(path + "BRF.png");
  TCanvas *c1 = (TCanvas *)f1->Get("c1");
  c1->Draw();
  c1->SaveAs(path + "Fe.png");
  TCanvas *c1b = (TCanvas *)f1->Get("c1b");
  c1b->Draw();
  c1b->SaveAs(path + "Fmu.png");
  TCanvas *c2 = (TCanvas *)f1->Get("c2");
  c2->Draw();
  c2->SaveAs(path + "BR.png");
  TCanvas *c3 = (TCanvas *)f1->Get("c3");
  c3->Draw();
  c3->SaveAs(path + "Gamma.png");
}
