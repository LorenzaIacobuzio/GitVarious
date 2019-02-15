void TreePlots(TString dir, TString histo1) {

  // dir = output dir, histo1 = histo to do cosmetics on, data = data or mc, mode: all = all dirs, hn = HeavyNeutrino dir, hnss = HeavyNeutrinoScan/SingleValue, hnsc = HeavyNeutrinoScan/Coupling, hnsm = HeavyNeutrinoScan/Mass, hnst = HeavyNeutrinoScan/Total

  TCanvas *c = new TCanvas();  
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;  
  TString path = "";

  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.25);
  c->SetTopMargin(0.15);
  c->SetGrid();
  c->RedrawAxis();

  if (dir != "")
    path = dir;
  else
    path = "/home/li/cernbox/PhD/Talks and papers/Notes/MCnote/images/Plots/"; ///Users/lorenza/cernbox/PhD/Talks and papers/Notes/MCnote/images/Plots/
  
  if (histo1.Contains("1"))
    path += "1/";
  else if (histo1.Contains("2"))
    path += "2/";
  else if (histo1.Contains("3"))
    path += "3/";
  else if (histo1.Contains("Data"))
    path += "Data/";
  else if (histo1.Contains("POT"))
    path += "POT/";
  else
    path += "2/";

  TFile *f = TFile::Open(histo1);

  if (f == 0) {
    cout << "Error: cannot open " << histo1 << endl;
    return;
  }

  TTree* tree = (TTree*)f->Get("HeavyNeutrino/Passed");
  Double_t Weight;
  Double_t CHODTime1;
  Double_t CHODTime2;
  Double_t CDA;
  Double_t Zvertex;
  Double_t CDALine;
  Double_t ZCDALine;
  Double_t BeamlineDist;
  Double_t xSR;
  Double_t ySR;
  Double_t MuEoP;
  Double_t PiEoP;
  Double_t R;
  Double_t energyPi;
  Double_t energyMu;
  Double_t invMass;
  TVector3 *Mom1 = new TVector3();
  TVector3 *Mom2 = new TVector3();
  TVector3 *TotMom = new TVector3();
  TVector3 *Vertex = new TVector3();
  TVector3 *threeMomPi = new TVector3();
  TVector3 *threeMomMu = new TVector3();
  Bool_t Target;
  Int_t Assoc;

  tree->SetBranchAddress("Weight", &Weight);
  tree->SetBranchAddress("CHODTime1", &CHODTime1);
  tree->SetBranchAddress("CHODTime2", &CHODTime2);
  tree->SetBranchAddress("CDA", &CDA);
  tree->SetBranchAddress("Zvertex", &Zvertex);
  tree->SetBranchAddress("CDALine", &CDALine);
  tree->SetBranchAddress("ZCDALine", &ZCDALine);
  tree->SetBranchAddress("BeamlineDist", &BeamlineDist);
  tree->SetBranchAddress("xSR", &xSR);
  tree->SetBranchAddress("ySR", &ySR);
  tree->SetBranchAddress("MuEoP", &MuEoP);
  tree->SetBranchAddress("PiEoP", &PiEoP);
  tree->SetBranchAddress("R", &R);
  tree->SetBranchAddress("energyPi", &energyPi);
  tree->SetBranchAddress("energyMu", &energyMu);
  tree->SetBranchAddress("invMass", &invMass);
  tree->SetBranchAddress("Mom1", &Mom1);
  tree->SetBranchAddress("Mom2", &Mom2);
  tree->SetBranchAddress("TotMom", &TotMom);
  tree->SetBranchAddress("Vertex", &Vertex);
  tree->SetBranchAddress("threeMomPi", &threeMomPi);
  tree->SetBranchAddress("threeMomMu", &threeMomMu);
  tree->SetBranchAddress("Target", &Target);
  tree->SetBranchAddress("Assoc", &Assoc);
    
  TH2D *h = new TH2D("h", "h", 100, 1.E-14, 1.E-12, 100, 0., 100.);
  
  for(Int_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    h->Fill(Weight, Mom1->X()/1000.);
  }

  h->Draw("text");
  gPad->SetLogx();
  c->SaveAs(path + "HeavyNeutrino/" + h->GetName() + ".pdf");
  c->SaveAs(path + "HeavyNeutrino/" + h->GetName() + ".png");
  tree->ResetBranchAddresses();
}
