#include "TProfile.h"
#include <fstream>
#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TGraphErrors.h"

void LoadStyle();
void SetLegend(TLegend *);

TH1D* projectHistogram(TH2F *, double , double);
TH1D* projectHistogramX(TH2F *, double , double);
TProfile* projectProfile(TH2F *, double , double);

const int nMultBins = 6;
double minMultBins[] = {0,10,20,30,40,50};
double maxMultBins[] = {10,20,30,40,50,120};
double nMultBin[] = {0,10,20,30,40,50,120};

void Mult_estimation_run3_pp()
{

  TFile *fInData  = new TFile("AnalysisResults_Data.root", "READ");
    // Same event pairing
    TList *listData = (TList*) fInData -> Get("table-maker/output");
    TList *listMultData = (TList*) listData -> FindObject("Event_AfterCuts");
    listMultData->ls();
   if (!listMultData) {
        std::cout << "Failed to retrieve TList: Event_AfterCuts" << std::endl;
        return;
    }

   TH2F* hist_VtxZ_NcontribReal = (TH2F*) listMultData->FindObject("VtxZ_VtxNContribReal");
   //TCanvas *cData = new TCanvas("cData", "Response Matrix", 800, 600);
   //TH1F *histMultData = (TH1F*) projectHistogram(hist_VtxZ_NcontribReal,0,100);
   //histMultData->Draw();

   
    TFile *fIn  = new TFile("AnalysisResults_Mult.root", "READ");
    // Same event pairing
    TList *list1 = (TList*) fIn -> Get("table-maker-m-c/output");
    TList *listMult = (TList*) list1 -> FindObject("Event_AfterCuts");

   if (!listMult) {
        std::cout << "Failed to retrieve TList: Event_AfterCuts" << std::endl;
        return;
    }

   // Create a vector to store slope parameters
    TVectorD slopes(nMultBins);
    TVectorD slopes_err(nMultBins);

    // Stores the mean of each projection
    TVectorD meanNch(nMultBins); 
    TVectorD meanNch_err(nMultBins);

    TVectorD meanNTrk(nMultBins); // Stores the mean of each projection
    TVectorD meanNTrk_err(nMultBins);


    // Define colors for the fits (one for each histogram)
    int colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan,kPink,kYellow,kOrange,kViolet,kRed+2,kBlue+2};

   TH1F* hist_recNtr = (TH1F*) listMult->FindObject("VtxNContribReal");
   TH1F* hist_Nch = (TH1F*) listMult->FindObject("MultMCNParticlesEta08");
   TH2F* hist_recNtr_Nch = (TH2F*) listMult->FindObject("VtxNContribReal_MultMCNParticlesEta08");
   TH2F *hist_recNtr_Nch2 = (TH2F*)hist_recNtr_Nch->Clone();
   hist_recNtr_Nch ->SetName("recNtr_Nch");
   hist_recNtr_Nch -> SetTitle("recNtr_Nch");
   hist_recNtr_Nch->GetXaxis()->SetTitle(" recN_{Trk} ");
   hist_recNtr_Nch->GetYaxis()->SetTitle(" genN_{ch} ");

  
   TH1D *histMult[nMultBins];
   TH1D *histMultData[nMultBins];
   TProfile *profMult[nMultBins];
   TF1 *fitFuncMult[nMultBins];
    for (int imult = 0;imult < nMultBins;imult++) {
            histMultData[imult] = (TH1D*) projectHistogramX(hist_VtxZ_NcontribReal, minMultBins[imult], maxMultBins[imult]);
	    histMultData[imult] -> SetName(Form("histMultData_%0.1f_%0.1f", minMultBins[imult], maxMultBins[imult]));
            histMultData[imult] ->SetMarkerStyle(20);
	    histMultData[imult] ->SetMarkerSize(0.5);
	    histMultData[imult] ->SetMarkerColor(colors[imult]);

	    double mNTrk = histMultData[imult]->GetMean();
	    double mNTrk_er = histMultData[imult]->GetMeanError();

	    // Store the results in vectors
	    meanNTrk[imult] = mNTrk;
	    meanNTrk_err[imult] = mNTrk_er;
	    


	    histMult[imult] = (TH1D*) projectHistogram(hist_recNtr_Nch2, minMultBins[imult], maxMultBins[imult]);
	    histMult[imult] -> SetName(Form("histMult_%0.1f_%0.1f", minMultBins[imult], maxMultBins[imult]));
            histMult[imult] ->SetMarkerStyle(20);
	    histMult[imult] ->SetMarkerSize(0.5);
	    histMult[imult] ->SetMarkerColor(colors[imult]);

	    double mNch = histMult[imult]->GetMean();
	    double mNch_er = histMult[imult]->GetMeanError();

	    // Store the results in vectors
	    meanNch[imult] = mNch;
	    meanNch_err[imult] = mNch_er;


	    
	    profMult[imult] = (TProfile*) projectProfile(hist_recNtr_Nch2, minMultBins[imult], maxMultBins[imult]);
            profMult[imult] -> SetName(Form("profMult_%0.1f_%0.1f", minMultBins[imult], maxMultBins[imult]));
	    
	    // Define the linear function for fitting y = mx
	    fitFuncMult[imult] = new TF1(Form("fitMult_%0.1f_%0.1f", minMultBins[imult], maxMultBins[imult]), "[0]*x", profMult[imult]->GetXaxis()->GetXmin(), profMult[imult]->GetXaxis()->GetXmax());
            fitFuncMult[imult]->SetLineColor(colors[imult]);
	    // Perform the fit
	    profMult[imult]->Fit(fitFuncMult[imult], "REMQ");  // "Q" stands for quiet mode, suppress fit output
	    // Retrieve the slope parameter (m)
	    double m = fitFuncMult[imult]->GetParameter(0);
	    double m_err = fitFuncMult[imult]->GetParError(0);
	    // Store the slope parameter in the vector
	    slopes[imult] = m;
	    slopes_err[imult] = m_err;
    }

    
    // Optional: Draw the TProfile
    TCanvas *c1 = new TCanvas("c1", "Response Matrix", 800, 600);
    TProfile *profileNch = hist_recNtr_Nch->ProfileX("profileNch");
    hist_recNtr_Nch->Draw("colz");
    profileNch->Draw("psame");
    for (int imult = 0;imult < nMultBins;imult++) {
      fitFuncMult[imult]->Draw("same");
	}

    // Create a canvas to draw all the TH1F projections on the same canvas
    TCanvas *c2 = new TCanvas("c2", "Fits of TH1F from TH2F", 800, 600);
    c2->SetLogy();
    hist_recNtr->Draw("p");
     // Loop over the Ntrk intervals and create projections for each
    for (int imult = 0;imult < nMultBins;imult++) {
        // Draw the projection on the canvas
          histMult[imult]->Draw("psame");
    }

    
    // Output the slope parameters stored in the vector
    int nBins = sizeof(nMultBin) / sizeof(nMultBin[0]) - 1;  // Number of bins
    TH1F *hist_alpha = new TH1F("hist_alpha", "Ntrkk vs alpha", nBins, nMultBin);
    hist_alpha->GetXaxis()->SetTitle(" N_{Trk} ");
    hist_alpha->GetYaxis()->SetTitle(" #alpha ");
    hist_alpha ->SetMarkerStyle(20);
    hist_alpha ->SetMarkerSize(0.5);
    hist_alpha ->SetMarkerColor(kRed);
    std::cout << "Slope parameters stored in the vector:" << std::endl;
    for (int i = 0; i < nMultBins; ++i) {
      hist_alpha->SetBinContent(i+1,slopes[i]);
      hist_alpha->SetBinError(i+1,slopes_err[i]);
      std::cout << "Slope for multiplicity " << minMultBins[i]<<"-"<< maxMultBins[i] << ": " << slopes[i] <<"+/-"<<slopes_err[i]<<std::endl;
    }


    // Output the mean and error values
    TH1F *hist_mNch = new TH1F("hist_mNch", "Ntrk vs <Nch>", nBins, nMultBin);
    hist_mNch->GetXaxis()->SetTitle(" N_{Trk} ");
    hist_mNch->GetYaxis()->SetTitle(" <N_{ch}>");
    hist_mNch ->SetMarkerStyle(20);
    hist_mNch ->SetMarkerSize(0.5);
    hist_mNch ->SetMarkerColor(kRed);
    std::cout << "Mean and Error for each Ntrk interval:" << std::endl;
    for (int i = 0; i < nMultBins; ++i) {
	hist_mNch->SetBinContent(i+1,meanNch[i]);
        hist_mNch->SetBinError(i+1,meanNch_err[i]);
        std::cout << " Nch Interval " << minMultBins[i]<<"-"<< maxMultBins[i] << ": <Nch> = " << meanNch[i] << ", <Nch> Error = " << meanNch_err[i] << std::endl;
    }

    // Output the mean Ntrk from data and error values
    TH1F *hist_mTrk = new TH1F("hist_mTrk", "Ntrk vs alpha*<NTrk>", nBins, nMultBin);
    hist_mTrk->GetXaxis()->SetTitle(" N_{Trk} ");
    hist_mTrk->GetYaxis()->SetTitle(" #alpha #time <N_{Trk}> / <N_{ch}>");
    hist_mTrk ->SetMarkerStyle(20);
    hist_mTrk ->SetMarkerSize(0.5);
    hist_mTrk ->SetMarkerColor(kBlue);
    std::cout << "Mean and Error for each Ntrk interval:" << std::endl;
    for (int i = 0; i < nMultBins; ++i) {
	hist_mTrk->SetBinContent(i+1,slopes[i]*meanNTrk[i]);
        hist_mTrk->SetBinError(i+1,slopes[i]*meanNTrk_err[i]);
        std::cout << " NTrk * alpha : " << minMultBins[i]<<"-"<< maxMultBins[i] << ": <Nch> = " << slopes[i]*meanNTrk[i] << ", <NTrk> Error = " << slopes[i]*meanNTrk_err[i] << std::endl;
    }


    
    TCanvas *cNch = new TCanvas("cNch", "alpha and <Nch>",220,130,1188,398);
    cNch->Divide(3,1);
    cNch->cd(1);
    hist_alpha->Draw("p");
    cNch->cd(2);
    hist_mNch->Draw("p");
    cNch->cd(3);
    hist_mTrk->Draw("p");
    hist_mNch->Draw("psame");
    
     
}

TH1D* projectHistogramX(TH2F *hist2D, double minMultRange, double maxMultRange) {
    double minMultBin = hist2D -> GetXaxis() -> FindBin(minMultRange);
    double maxMultBin = hist2D -> GetXaxis() -> FindBin(maxMultRange - 0.00001);

    Printf("minMultBin = %0.1f, maxMultBin = %0.1f", minMultBin, maxMultBin);
    hist2D -> GetXaxis() -> SetRange(minMultBin, maxMultBin);
   
    TH1D *histProj = (TH1D*) hist2D -> ProjectionX(Form("histProj_%0.1f_%0.1f", minMultRange, maxMultRange),minMultBin, maxMultBin);
    return histProj;
}

TH1D* projectHistogram(TH2F *hist2D, double minMultRange, double maxMultRange) {
    double minMultBin = hist2D -> GetXaxis() -> FindBin(minMultRange);
    double maxMultBin = hist2D -> GetXaxis() -> FindBin(maxMultRange - 0.00001);

    Printf("minMultBin = %0.1f, maxMultBin = %0.1f", minMultBin, maxMultBin);
    hist2D -> GetXaxis() -> SetRange(minMultBin, maxMultBin);
   
    TH1D *histProj = (TH1D*) hist2D -> ProjectionY(Form("histProj_%0.1f_%0.1f", minMultRange, maxMultRange),minMultBin, maxMultBin);
    return histProj;
}


TProfile* projectProfile(TH2F *hist2D, double minMultRange, double maxMultRange) {
  double minMultBin = hist2D -> GetXaxis() -> FindBin(minMultRange);
    double maxMultBin = hist2D -> GetXaxis() -> FindBin(maxMultRange - 0.00001);
    
    Printf("minMultBin = %0.1f, maxMultBin = %0.1f", minMultBin, maxMultBin);
    hist2D -> GetXaxis() -> SetRange(minMultBin, maxMultBin);

    TProfile *histProj = (TProfile*) hist2D -> ProfileX(Form("histProj__%0.1f_%0.1f", minMultRange, maxMultRange),minMultBin, maxMultBin);
    return histProj;
}
