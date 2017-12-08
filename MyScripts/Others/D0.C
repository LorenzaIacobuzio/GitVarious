void D0() {
  TFile *f1 = TFile::Open("/home/li/Desktop/charmed_10MEvts_tot_8.2_ok.root");
  f1->cd();
  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  Double_t Pl, Pt, val, sum, N;
  ofstream outfile("PtVsPlD0.txt");
  
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!strcmp(key->GetName(), "D0PlPt")) {
      TH2 *h2 = (TH2*)key->ReadObj();
      N = h2->GetEntries();
      for (Int_t i = 1; i < h2->GetNbinsX(); i++) {
	for (Int_t j = 1; j < h2->GetNbinsY(); j++) {
	  Pl = h2->GetXaxis()->GetBinCenter(i);
	  Pt = h2->GetYaxis()->GetBinCenter(j);
	  val = h2->GetBinContent(i, j)/N;
	  outfile << Pl << " " << Pt << " " << val << "\n";
	}
      }
    }
  }
}
