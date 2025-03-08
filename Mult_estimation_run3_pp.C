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
TProfile* projectProfileY(TH2F *, double , double);
double EvalError(TF1* , double );

const int nMultBins = 10;
double minMultBins[] = {1,9,15,21,26,34,42,51,61,81};
double maxMultBins[] = {8,14,20,25,33,41,50,60,80,115};
double nMultBin[] = {1,9,15,21,26,34,42,51,61,81,115};

bool isLinear = true;
bool isPol3 = false;

//for plotting
std::string bins[] = {"1-8", "9-14", "15-20", "21-25", "26-33", 
		      "34-41", "42-50", "51-60", "61-80", "81-115"};

struct StoreMeanNch {
  double mean;
  double error;
  std::string bin;
};

void Mult_estimation_run3_pp()
{

  std::vector<StoreMeanNch> dataNch;
  std::vector<StoreMeanNch> MCNch;
  std::vector<StoreMeanNch> alphaNch;
  std::vector<StoreMeanNch> pol3Nch;
  std::vector<StoreMeanNch> pol1Nch;

    
  TFile *fInData  = new TFile("AnalysisResults_Data.root", "READ");
  // Same event pairing
  TList *listData = (TList*) fInData -> Get("table-maker/output");
  TList *listMultData = (TList*) listData -> FindObject("Event_AfterCuts");
  listMultData->ls();
  if (!listMultData) {
    std::cout << "Failed to retrieve TList: Event_AfterCuts" << std::endl;
    return;
  }

  TH2F* hist_VtxZ_NcontribReal = (TH2F*) listMultData->FindObject("VtxZ_VtxNcontribReal");
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

  //Get the linear TProfile and fit polynomials
  TCanvas *c0 = new TCanvas("c0", "Integrated response matrix Fit", 800, 600);
  TProfile *profMultData_integrated = (TProfile*) projectProfile(hist_recNtr_Nch2,1,120);
  TF1 *fitMult_pol1 = new TF1("fitMult_pol1","pol1",1,80);
  fitMult_pol1->SetLineColor(kBlack);
  TF1 *fitMult_pol3 = new TF1("fitMult_pol3","pol3",1,80);
  fitMult_pol3->SetLineColor(kBlack);
  fitMult_pol3->SetLineStyle(6);
  profMultData_integrated->Fit(fitMult_pol1, "REM"," ",1,80);
  profMultData_integrated->Fit(fitMult_pol3, "REM"," ",1,80);

  
  TH1D *histMult[nMultBins];
  TH1D *histMultData[nMultBins];
  TProfile *profMult[nMultBins];
  TF1 *fitFuncMult[nMultBins];
  for (int imult = 0;imult < nMultBins;imult++) {
    histMultData[imult] = (TH1D*) projectProfileY(hist_VtxZ_NcontribReal, minMultBins[imult], maxMultBins[imult]);
    histMultData[imult] -> SetName(Form("histMultData_%0.1f_%0.1f", minMultBins[imult], maxMultBins[imult]));
    histMultData[imult] ->SetMarkerStyle(20);
    histMultData[imult] ->SetMarkerSize(0.5);
    histMultData[imult] ->SetMarkerColor(colors[imult]);
    dataNch.push_back({histMultData[imult]->GetMean(2), histMultData[imult]->GetMeanError(2), bins[imult]}); //<Ntrk>, <Ntrk> error


	    
    histMult[imult] = (TH1D*) projectHistogram(hist_recNtr_Nch2, minMultBins[imult], maxMultBins[imult]);
    histMult[imult] -> SetName(Form("histMult_%0.1f_%0.1f", minMultBins[imult], maxMultBins[imult]));
    histMult[imult] ->SetMarkerStyle(20);
    histMult[imult] ->SetMarkerSize(0.5);
    histMult[imult] ->SetMarkerColor(colors[imult]);
    MCNch.push_back({histMult[imult]->GetMean(), histMult[imult]->GetMeanError(), bins[imult]}); // MC <Nch> , error
	    


	    
    profMult[imult] = (TProfile*) projectProfile(hist_recNtr_Nch2, minMultBins[imult], maxMultBins[imult]);
    profMult[imult] -> SetName(Form("profMult_%0.1f_%0.1f", minMultBins[imult], maxMultBins[imult]));	    
    // Define the linear function for fitting y = mx
    if(isLinear) fitFuncMult[imult] = new TF1(Form("fitMult_%0.1f_%0.1f", minMultBins[imult], maxMultBins[imult]), "[0]*x", profMult[imult]->GetXaxis()->GetXmin(), profMult[imult]->GetXaxis()->GetXmax());	    
    fitFuncMult[imult]->SetLineColor(colors[imult]);
    profMult[imult]->Fit(fitFuncMult[imult], "REM"," ",minMultBins[imult], maxMultBins[imult]);  // "Q" stands for quiet mode, suppress fit output
    // Retrieve the slope parameter (m)
    double m = fitFuncMult[imult]->GetParameter(0);
    double m_err = fitFuncMult[imult]->GetParError(0);
    slopes[imult] = m;
    slopes_err[imult] = m_err;

    //alpha X <Ntrk>
    alphaNch.push_back({histMultData[imult]->GetMean(2),histMultData[imult]->GetMeanError(2), bins[imult]}); //alpha * <NTrk>
    pol1Nch.push_back({fitMult_pol1->Eval(histMultData[imult]->GetMean(2)), EvalError(fitMult_pol1,histMultData[imult]->GetMean(2)), bins[imult]}); //Pol1.Evaluation, error
    pol3Nch.push_back({fitMult_pol3->Eval(histMultData[imult]->GetMean(2)), EvalError(fitMult_pol3,histMultData[imult]->GetMean(2)), bins[imult]}); //Pol1.Evaluation, error

  }

    
  // Optional: Draw the TProfile
  TCanvas *c1 = new TCanvas("c1", "Response Matrix", 800, 600);
  TProfile *profileNch = hist_recNtr_Nch->ProfileX("profileNch");
  hist_recNtr_Nch->Draw("colz");
  profileNch->Draw("psame");
  fitMult_pol1->Draw("lsame");
  fitMult_pol3->Draw("lsame");
   
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
  cout<<endl;
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
  for (int i = 0; i < nMultBins; ++i) {
    hist_mNch->SetBinContent(i+1,MCNch[i].mean);
    hist_mNch->SetBinError(i+1,MCNch[i].error);
  }

    
  // Output the mean Ntrk from data and error values
  TH1F *hist_mTrk = new TH1F("hist_mTrk", "<Nch> vs alpha*<NTrk>", nBins, nMultBin);
  hist_mTrk->GetXaxis()->SetTitle(" N_{Trk} ");
  hist_mTrk->GetYaxis()->SetTitle(" #alpha #time <N_{Trk}> / <N_{ch}>");
  hist_mTrk ->SetMarkerStyle(20);
  hist_mTrk ->SetMarkerSize(0.5);
  hist_mTrk ->SetMarkerColor(kBlue);
  cout<<endl;
  for (int i = 0; i < nMultBins; ++i) {
    hist_mTrk->SetBinContent(i+1,alphaNch[i].mean);
    hist_mTrk->SetBinError(i+1,alphaNch[i].error);
  }


  cout<<" MC Truth <Nch> "<<endl;
  for (const auto& point : MCNch) {
    std::cout << " Bin: " << point.bin << "  Mean: " << point.mean << ", Error: " << point.error << std::endl;
  }

  cout<<" alpha times <Ntrk> "<<endl;
  for (const auto& point : alphaNch) {
    std::cout << " Bin: " << point.bin << "  Mean: " << point.mean << ", Error: " << point.error << std::endl;
  }


  TCanvas *c3 = new TCanvas("c3", "<Ntrk> profile in mult bins", 800, 600);
  c3->Divide(nMultBins,1);
  for (int i = 0; i < nMultBins; ++i) {
    histMultData[i]->SetTitle(Form("Mult_%0.1f_%0.1f", minMultBins[i], maxMultBins[i]));
    c3->cd(i+1);
    histMultData[i]->Draw("p");
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


 // Arrays to store means and variances
  std::vector<double> means(nMultBins, 0.0);
  std::vector<double> systematics(nMultBins, 0.0);
  int nChTrials = 3;
  // Calculate mean for each bin
  cout<<endl;
  for (int i = 0; i < nMultBins; ++i) {
    double sum = alphaNch[i].mean + pol3Nch[i].mean + pol1Nch[i].mean;
    means[i] = sum / nChTrials;  // Mean of the four methods
    // Calculation of variance
    double sum_squared_diff = 
      pow(alphaNch[i].mean - means[i], 2) + 
      pow(pol3Nch[i].mean - means[i], 2) + 
      pow(pol1Nch[i].mean - means[i], 2);
        
     systematics[i] = TMath::Sqrt(sum_squared_diff / nChTrials);  // Variance
  }

  // Output results
  std::cout << "Bin\tMean\tSystematics\n";
  for (int i = 0; i < nMultBins; ++i) {
    std::cout << "Bin " << i << "\t" << means[i] << "\t" << systematics[i] << "\n";
  }


   std::cout << "\n Ratio MC<Nch>/<Eval<Nch> \n";
  for (int i = 0; i < nMultBins; ++i) {
    std::cout << " Bin " << i <<"\t"<< MCNch[i].mean/means[i] << "\n";
  }


    
    
     
}

TH1D* projectHistogramX(TH2F *hist2D, double minMultRange, double maxMultRange) {
  double minMultBin = hist2D -> GetXaxis() -> FindBin(minMultRange);
  double maxMultBin = hist2D -> GetXaxis() -> FindBin(maxMultRange);

  Printf("minMultBin = %0.1f, maxMultBin = %0.1f", minMultBin, maxMultBin);
  hist2D -> GetXaxis() -> SetRange(minMultBin, maxMultBin);
   
  TH1D *histProj = (TH1D*) hist2D -> ProjectionX(Form("histProj_%0.1f_%0.1f", minMultRange, maxMultRange),minMultBin, maxMultBin);
  return histProj;
}

TH1D* projectHistogram(TH2F *hist2D, double minMultRange, double maxMultRange) {
  double minMultBin = hist2D -> GetXaxis() -> FindBin(minMultRange);
  double maxMultBin = hist2D -> GetXaxis() -> FindBin(maxMultRange);

  Printf("minMultBin = %0.1f, maxMultBin = %0.1f", minMultBin, maxMultBin);
  hist2D -> GetXaxis() -> SetRange(minMultBin, maxMultBin);
   
  TH1D *histProj = (TH1D*) hist2D -> ProjectionY(Form("histProj_%0.1f_%0.1f", minMultRange, maxMultRange),minMultBin, maxMultBin);
  return histProj;
}


TProfile* projectProfile(TH2F *hist2D, double minMultRange, double maxMultRange) {
  double minMultBin = hist2D -> GetXaxis() -> FindBin(minMultRange);
  double maxMultBin = hist2D -> GetXaxis() -> FindBin(maxMultRange);
    
  Printf("minMultBin = %0.1f, maxMultBin = %0.1f", minMultBin, maxMultBin);
  hist2D -> GetXaxis() -> SetRange(minMultBin, maxMultBin);

  TProfile *histProj = (TProfile*) hist2D -> ProfileX(Form("histProj__%0.1f_%0.1f", minMultRange, maxMultRange),minMultBin, maxMultBin);
  return histProj;
}


TProfile* projectProfileY(TH2F *hist2D, double minMultRange, double maxMultRange) {
  double minMultBin = hist2D -> GetYaxis() -> FindBin(minMultRange);
  double maxMultBin = hist2D -> GetYaxis() -> FindBin(maxMultRange);
    
  Printf("minMultBin = %0.1f, maxMultBin = %0.1f", minMultBin, maxMultBin);
  hist2D -> GetYaxis() -> SetRange(minMultBin, maxMultBin);

  TProfile *histProj = (TProfile*) hist2D -> ProfileX(Form("histProj__%0.1f_%0.1f", minMultRange, maxMultRange),minMultBin, maxMultBin);
  return histProj;
}


// Function to calculate the error in the evaluated function value at x_value
double EvalError(TF1* f, double x_value) {
  // Get the number of parameters of the TF1 function
  int n_params = f->GetNpar();
    
  // Initialize the error to 0
  double error = 0.0;
    
  // Loop over all the parameters to calculate the error
  for (int i = 0; i < n_params; ++i) {
    // Get the current parameter value and its error
    double value = f->GetParameter(i);
    double par_error = f->GetParError(i);
        
    // Temporarily modify the parameter to calculate the derivative
    double original_value = value; // Save the original parameter value
        
    // Compute the function's value with a small perturbation
    double perturb = 1e-6; // Small perturbation
    f->SetParameter(i, original_value + perturb);
    double f_plus = f->Eval(x_value);
        
    f->SetParameter(i, original_value - perturb);
    double f_minus = f->Eval(x_value);
        
    f->SetParameter(i, original_value); // Restore the original value
        
    // Calculate the numerical derivative
    double partial_derivative = (f_plus - f_minus) / (2 * perturb);
        
    // Add the squared contribution to the total error
    error += TMath::Power(partial_derivative * par_error, 2);
  }
    
  // Return the square root of the accumulated error
  return TMath::Sqrt(error);
}

