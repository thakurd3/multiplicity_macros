
void SetStyle(Bool_t graypalette=kFALSE);
TH1D* ScaleXaxisHisto(TH1D* histogram, Double_t scaleFactor);
void Ntrk_Dist()
{
  SetStyle();
  
  TCanvas *c1 = new TCanvas("c1", "c1",207,89,625,491);
  gStyle->SetOptStat(0);
  c1->SetLogy();
  c1->Draw();

  TH2F *r1 = new TH2F("r1","",14,0.0,14.0,100,0.1,100000000.0);
  r1->GetYaxis()->SetTitle("Events");
  r1->GetXaxis()->SetTitle("N_{tracks}/<N_{tracks}>");
  r1->Draw();

  TFile *file = new TFile("../AnalysisResults.root");

  TList *l1 = (TList*) file->Get("table-maker/output");
  TList *l_MEPM1 = (TList*)l1->FindObject("Event_AfterCuts");
  TH1D *Mult_ITS = (TH1D*)l_MEPM1->FindObject("MultNTracksHasITS");
  TH1D *Mult_TPC = (TH1D*)l_MEPM1->FindObject("MultNTracksHasTPC");
  TH1D *Mult_ITS_TPC = (TH1D*)l_MEPM1->FindObject("MultNTracksITSTPC");
  //MultNTracksHasITS

  double mean_its = Mult_ITS->GetMean();
  double mean_tpc = Mult_TPC->GetMean();
  double mean_its_tpc = Mult_ITS_TPC->GetMean();
  cout<<" mean_its "<<mean_its<< " mean_tpc "<<mean_tpc <<" mean_its_tpc "<<mean_its_tpc<<endl;
  TH1D* scaleHist_its = ScaleXaxisHisto(Mult_ITS, mean_its);
  TH1D* scaleHist_tpc = ScaleXaxisHisto(Mult_TPC, mean_tpc);
  TH1D* scaleHist_its_tpc = ScaleXaxisHisto(Mult_ITS_TPC, mean_its_tpc);

  scaleHist_its->Scale(1.0/Mult_ITS->GetEntries());
  scaleHist_tpc->Scale(1.0/Mult_TPC->GetEntries());
  scaleHist_its_tpc->Scale(1.0/Mult_ITS_TPC->GetEntries());

  scaleHist_its->GetYaxis()->SetTitle("Events");
  scaleHist_its->GetXaxis()->SetTitle("N_{tracks}/<N_{tracks}>");
  scaleHist_its->SetMarkerStyle(71);
  scaleHist_its->SetLineColor(kRed);
  scaleHist_its->SetMarkerColor(kRed);

  scaleHist_tpc->SetMarkerStyle(33);
  scaleHist_tpc->SetLineColor(kBlack);
  scaleHist_tpc->SetMarkerColor(kBlack);
  
  scaleHist_its_tpc->SetMarkerStyle(24);
  scaleHist_its_tpc->SetLineColor(kViolet);
  scaleHist_its_tpc->SetMarkerColor(kViolet);

  scaleHist_its->Draw("p");
  scaleHist_tpc->Draw("psame");
  scaleHist_its_tpc->Draw("psame");

  TLegend *leg1 = new TLegend(0.53,0.62,0.88,0.91);
  leg1->SetTextSize(0.04);
  leg1->SetFillStyle(0);
  leg1->SetHeader("pp, 13.6 TeV, |#eta| < 0.8");
  leg1->AddEntry(scaleHist_its, "MultNTracksHasITS","p");
  leg1->AddEntry(scaleHist_tpc, "MultNTracksHasTPC","p");
  leg1->AddEntry(scaleHist_its_tpc, "MultNTracksITSTPC","p");
  leg1->Draw();
  

}


TH1D* ScaleXaxisHisto(TH1D* histogram, Double_t scaleFactor) {
TH1D *newhist = new TH1D("newhist","",1500,0.0,15);
  for (int i = 1; i <= histogram->GetNbinsX(); i++) {
        Double_t binContent = histogram->GetBinContent(i);
        Double_t binCenter = histogram->GetBinCenter(i);
	Double_t new_bincCenter = binCenter/scaleFactor;
	Int_t new_bin = newhist->GetXaxis()->FindBin(new_bincCenter);
        newhist->SetBinContent(new_bin, binContent);
    }
  return newhist;
}

void SetStyle(Bool_t graypalette) {
  //cout << "Setting style!" << endl;
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kBlack);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.009,"y");
  gStyle->SetLabelOffset(0.007,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(1.08,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetTickLength(0.03,"X");
  gStyle->SetTickLength(0.03,"Y");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
}



 
# multiplicity_macros
