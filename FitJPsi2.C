#include "TProfile.h"
#include <fstream>
#include <iostream>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TH2D.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include <TLatex.h>
#include <TStyle.h>
#include"TMatrixDSym.h"
#include"TFitResult.h"
Bool_t reject;
void SetStyle(Bool_t graypalette=kTRUE);

Double_t BackgroundVWG(Double_t *x, Double_t *par);
Double_t CrystalBallExtended(Double_t *x,Double_t *par);
Double_t CrystalBallExtended2(Double_t *x,Double_t *par);
Double_t fitFunctionVWG(Double_t *x, Double_t *par);
Double_t fitFunctionCB2VWG(Double_t *x, Double_t *par);
Double_t FuncBck1(Double_t *x, Double_t *par);
Double_t FuncBck2(Double_t *x, Double_t *par);
Double_t BackgroundBCK(Double_t *x, Double_t *par);
Double_t fitFunctionBCK(Double_t *x, Double_t *par);
Double_t fitFunctionPol3(Double_t *x, Double_t *par);
Double_t fitFunction2CB2VWG(Double_t *x, Double_t *par);
Double_t BackgroundPol3(Double_t *x, Double_t *par);
Double_t fitFunction2CB2Pol3(Double_t *x, Double_t *par);

void FitJPsi2()
{
  Double_t FitLow = 2.2;
  Double_t FitHigh = 4.7;
  TCanvas *c = new TCanvas("c","c",20,20,600,600);
  //c->SetLogy();
  SetStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);


  TPaveText *display1 = new TPaveText(4.2867031,3915.316,5.0608323,6109.354,"BLARC");
  display1->SetTextFont(62);
  display1->SetTextSize(0.03);
  display1->SetTextColor(kBlack);
  display1->SetBorderSize(0);
  display1->SetFillColor(0);

  TPaveText *display2 = new TPaveText(4.3318969,963.6719,4.8850661,1293.534,"BLARC");
  display2->SetTextFont(62);
  display2->SetTextSize(0.03);
  display2->SetTextColor(kBlack);
  display2->SetBorderSize(0);
  display2->SetFillColor(0);


  TPaveText *display3 = new TPaveText(4.03318969,1300.6719,4.850661,1893.534,"BLARC");
  display3->SetTextFont(62);
  display3->SetTextSize(0.03);
  display3->SetTextColor(kBlack);
  display3->SetBorderSize(0);
  display3->SetFillColor(0);

  
  TH2F *f0=new TH2F("f","M_{#mu#mu}(GeV/c^{2})",120,2.0,5.0,1.e+07,1.6,1.e+07);
  f0->GetXaxis()->SetTitle("M_{#mu#mu}(GeV/c^{2})");
  f0->GetYaxis()->SetTitle("dN/dM_{#mu#mu} (c^{2}/GeV)");
  f0->GetXaxis()->CenterTitle();
  f0->GetYaxis()->CenterTitle();


  double minCentrBins[] = {10};
  double maxCentrBins[] = {50};
  
  double minPtBins[] = {0};
  double maxPtBins[] = {5};

  double mixDep[] = {6};
  double F_fact[5];
  double Norm_fact[5];
  
  TCanvas* c0 = new TCanvas("c0", "c0", 500, 500);
  c0->Divide(2,1,0);
  TFile *fOut = new TFile("Fit_Function_Cent10-50_Pt0_5Bin_new.root","RECREATE");
  TFile *f = new TFile(Form("Histograms_MixingDep_%1.0fmuonLowPt210SigmaPDCA_centr_10_50.root", mixDep[0]));
  TFile *f2 = new TFile("Histograms_matchedMchMid.root","READ");
  if (f->IsOpen()) cout << "File opened successfully" << endl;
  if (f2->IsOpen()) cout << "File opened successfully" << endl;

  TH1D *h_Mass_SEPM = (TH1D*) f2->Get("histMassSEPM_0_100__0_90");
  TH1D *h_Mass_MEPM = (TH1D*) f2->Get("histMassMEPM_0_100__0_90");
  TH1D *h_Mass = (TH1D*) f2->Get("histMassBkgSubtrSEPM_0_100__0_90");
  h_Mass_SEPM->GetXaxis()->SetRangeUser(2,5);
  //h_Mass_MEPM->GetXaxis()->SetRangeUser(2,5);  
  //h_Mass_SEPM->Write("Mass_SEPM",TObject::kWriteDelete);
 
  
  
  //h_Mass_SEPM->Draw("p");
  //h_Mass_MEPM->Draw("psame");

  
  double SE =h_Mass_SEPM->Integral(h_Mass_SEPM->GetXaxis()->FindBin(1.0),h_Mass_SEPM->GetXaxis()->FindBin(2.5)) + h_Mass_SEPM->Integral(h_Mass_SEPM->GetXaxis()->FindBin(3.72),h_Mass_SEPM->GetXaxis()->FindBin(5.0));

  double ME = h_Mass_MEPM->Integral(h_Mass_MEPM->GetXaxis()->FindBin(1.0),h_Mass_MEPM->GetXaxis()->FindBin(2.5)) + h_Mass_MEPM->Integral(h_Mass_MEPM->GetXaxis()->FindBin(3.72),h_Mass_MEPM->GetXaxis()->FindBin(5.0));

  cout<<" Normalization factor "<<ME/SE <<endl;

  h_Mass_MEPM->Scale(SE/ME);

  //h_Mass_SEPM->Draw("p");
  //h_Mass_MEPM->Draw("psame");

  
  h_Mass_SEPM->GetXaxis()->SetRangeUser(2,5);
  h_Mass_MEPM->GetXaxis()->SetRangeUser(2,5);
  h_Mass_SEPM->Add(h_Mass_MEPM,-1);

  //TH1D *h_Mass_MEPM = (TH1D*) f->Get(Form("histMassMEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[0], maxPtBins[0], minCentrBins[0], maxCentrBins[0]));
  //TH1D *h_Mass_MEPP = (TH1D*) f->Get(Form("histMassMEPP_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[0], maxPtBins[0], minCentrBins[0], maxCentrBins[0]));
  //TH1D *h_Mass_MEMM = (TH1D*) f->Get(Form("histMassMEMM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[0], maxPtBins[0], minCentrBins[0], maxCentrBins[0]));

  //TH1D *h_Mass_SEPM = (TH1D*) f->Get(Form("histMassSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[0], maxPtBins[0], minCentrBins[0], maxCentrBins[0]));
  //TH1D *h_Mass_SEPP = (TH1D*) f->Get(Form("histMassSEPP_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[0], maxPtBins[0], minCentrBins[0], maxCentrBins[0]));
  //TH1D *h_Mass_SEMM = (TH1D*) f->Get(Form("histMassSEMM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[0], maxPtBins[0], minCentrBins[0], maxCentrBins[0]));

   
  TH1F *hT = (TH1F*)h_Mass_SEPM->Clone();
 
  cout<<" Maximum  "<<hT->GetBinContent(hT->FindBin(3.096))<<endl;
  Double_t PeakPos = hT->GetBinContent(hT->FindBin(3.096));
  		       
  
  ///////////////Fitting //////////////////////////////////
  
  //TF1 *fitFctCB2VWG;
  TF1 *fitFctCB2VWG = new TF1("fitFctCB2VWG", fitFunction2CB2Pol3,FitLow,FitHigh, 12);
  fitFctCB2VWG->SetLineColor(kBlue);
  fitFctCB2VWG->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");
  fitFctCB2VWG->SetParName(11, "kPsi'");


  fitFctCB2VWG->SetParameter(0,200000); // background parameters
  fitFctCB2VWG->SetParameter(1,-0.0150326); // background parameters
  fitFctCB2VWG->SetParameter(2,0.00757168); // background parameters
  fitFctCB2VWG->SetParameter(3,-0.00112113); // background parameters
  

  
  //fitFctCB2VWG->SetParameter(0, 5000000.);
  //fitFctCB2VWG->SetParameter(1, 1.9);
  //fitFctCB2VWG->SetParameter(2, 0.5);
  //fitFctCB2VWG->SetParLimits(2, 0., 100.);
  //fitFctCB2VWG->SetParameter(3, 0.3);
  //fitFctCB2VWG->SetParLimits(3, 0., 100.);
  //fitFctCB2VWG->SetParameter(4, 100.);
  
  
  fitFctCB2VWG->SetParameter(4,240000);
  fitFctCB2VWG->SetParLimits(4,236000, 245001);
  //fitFctCB2VWG->SetParameter(4,3.2);
  fitFctCB2VWG->SetParameter(5, 3.09);
  fitFctCB2VWG->SetParLimits(5, 3.082, 3.1);
  fitFctCB2VWG->SetParameter(6, 0.072);
  fitFctCB2VWG->SetParLimits(6, 0.06, 0.15);

  fitFctCB2VWG->FixParameter(7,1.0169);
  //fitTotal->SetParLimits(7,0.1,10.0);
 
  fitFctCB2VWG->FixParameter(8,3.7064);
  //fitTotal->SetParLimits(8,0.0,10.0);
  
  fitFctCB2VWG->FixParameter(9, 2.4183);
  //fitTotal->SetParLimits(9,0.1,10.0);
  
  fitFctCB2VWG->FixParameter(10,7.9655);
  //fitTotal->SetParLimits(10,0.0,10.0);
  fitFctCB2VWG->SetParameter(11, 4000.0);
  //fitFctCB2VWG->SetParLimits(11,250,950);
  
  TFitResultPtr r = hT->Fit(fitFctCB2VWG,"RESLI");
  
  ///////////////////////////////////////////////////////////////
  

  TMatrixDSym cov = r->GetCovarianceMatrix();

  hT->Sumw2();
  //hT->GetXaxis()->SetRangeUser(2.0,4.2);
  //hT->GetYaxis()->SetRangeUser(1.7,2.4);
  //hT->SetMinimum(1.);
  //gPad->SetLogy(1);
  //hT->Draw("E1");
  //hT->SetMarkerStyle(20);
  //hT->SetMarkerColor(kBlack);
  //hT->SetMarkerSize(0.5);
  

 
  TF1 *fitFctCB2 = new TF1("fitFctCB2", CrystalBallExtended,FitLow,FitHigh, 7);
  fitFctCB2->SetLineColor(kRed);
  fitFctCB2->SetLineStyle(2);
  fitFctCB2->SetParNames("kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");
  fitFctCB2->FixParameter(0, fitFctCB2VWG->GetParameter(4)); 
  fitFctCB2->FixParameter(1, fitFctCB2VWG->GetParameter(5)); 
  fitFctCB2->FixParameter(2, fitFctCB2VWG->GetParameter(6));
  fitFctCB2->FixParameter(3, fitFctCB2VWG->GetParameter(7));   
  fitFctCB2->FixParameter(4, fitFctCB2VWG->GetParameter(8));   
  fitFctCB2->FixParameter(5, fitFctCB2VWG->GetParameter(9));   
  fitFctCB2->FixParameter(6, fitFctCB2VWG->GetParameter(10));
  fitFctCB2->Draw("same");
  //printf("\n --> nJPsi = %f\n\n", fitFctCB2->Integral(0., 100.)/hT->GetBinWidth(1));


  TF1 *fitFctCB21 = new TF1("fitFctCB21", CrystalBallExtended,FitLow,FitHigh, 7);
  fitFctCB21->SetLineColor(kRed);
  fitFctCB21->SetLineStyle(2);
  fitFctCB21->SetParNames("kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");
  fitFctCB21->SetParameter(0, fitFctCB2VWG->GetParameter(4)); 
  fitFctCB21->SetParameter(1, fitFctCB2VWG->GetParameter(5)); 
  fitFctCB21->SetParameter(2, fitFctCB2VWG->GetParameter(6));
  fitFctCB21->FixParameter(3, fitFctCB2VWG->GetParameter(7));   
  fitFctCB21->FixParameter(4, fitFctCB2VWG->GetParameter(8));   
  fitFctCB21->FixParameter(5, fitFctCB2VWG->GetParameter(9));   
  fitFctCB21->FixParameter(6, fitFctCB2VWG->GetParameter(10));
  TMatrixDSym covPeak = cov.GetSub(4,10,4,10);
  
  
  TF1 *fitFctCB22 = new TF1("fitFctCB22", CrystalBallExtended,FitLow,FitHigh, 7);
  fitFctCB22->SetLineColor(kRed);
  fitFctCB22->SetLineStyle(2);
  fitFctCB22->SetParNames("kPsi'","mPsi'","sPsi'","alPsi'","nlPsi'","auPsi'","nuPsi'");
  fitFctCB22->SetParameter(0, fitFctCB2VWG->GetParameter(11)); 

  //fitFctCB22->FixParameter(1, 3.68609+(fitFctCB2VWG->GetParameter(5)-3.096916)/3.096916*3.68609); 
  //fitFctCB22->FixParameter(2, fitFctCB2VWG->GetParameter(6)/3.096916*3.68609);

  fitFctCB22->SetParameter(1, (fitFctCB2VWG->GetParameter(5)+(3.68609-3.096916))); 
  fitFctCB22->SetParameter(2, fitFctCB2VWG->GetParameter(6)*(3.68609/3.096916));
  
  fitFctCB22->FixParameter(3, fitFctCB2VWG->GetParameter(7));   
  fitFctCB22->FixParameter(4, fitFctCB2VWG->GetParameter(8));   
  fitFctCB22->FixParameter(5, fitFctCB2VWG->GetParameter(9));   
  fitFctCB22->FixParameter(6, fitFctCB2VWG->GetParameter(10));
  TMatrixDSym covPeak2 = cov.GetSub(4,10,4,10);
  fitFctCB22->Draw("same");
  //printf(" --> nJPsi' = %f\n\n", fitFctCB22->Integral(0., 100.)/hT->GetBinWidth(1));
  

  TF1 *exp2 = new TF1("exp2",BackgroundPol3,2.,5.0,4);
  for(int i=0;i<4;i++) exp2->SetParameter(i,fitFctCB2VWG->GetParameter(i));
  exp2->SetLineColor(kBlack);
  exp2->SetLineStyle(2);
  TMatrixDSym covBck = cov.GetSub(0,3,0,3);
  exp2->Draw("same");

  cout<<"mass Jpsi" << fitFctCB2VWG->GetParameter(5) << " +/-" << fitFctCB2VWG->GetParError(5) << endl;
  cout<<"sigma Jpsi" << fitFctCB2VWG->GetParameter(6) << " +/-" << fitFctCB2VWG->GetParError(6) << endl;;
  cout << " mass of psi' "<< (3.68609+(fitFctCB2VWG->GetParameter(5)-3.096916)/3.096916*3.68609) << endl;
  cout << "sigma psi' " << (fitFctCB2VWG->GetParameter(6)/3.096916*3.68609) << endl;


  //====================J/Psi Signal Estimation==========================================================
 
  // Set("NofJPsi",njpsi,nerr);
  
  double m =fitFctCB2VWG->GetParameter(5);
  double m_err =fitFctCB2VWG->GetParError(5);
  double s =fitFctCB2VWG->GetParameter(6);
  double s_err =fitFctCB2VWG->GetParError(6);
  double njpsi3s = fitFctCB2->Integral(m-3*s,m+3*s)/hT->GetBinWidth(1);
  double nerr3s = fitFctCB21->IntegralError(m-3*s,m+3*s,fitFctCB21->GetParameters(),covPeak.GetMatrixArray())/hT->GetBinWidth(1);
  cout<< "No of J/Psi == " <<  njpsi3s <<"+/-" << nerr3s<<endl;
  Double_t chi=fitFctCB2VWG->GetChisquare();
  Double_t ndf=fitFctCB2VWG->GetNDF();
  cout<<" Chi2 "<<chi<<" NDF "<<ndf<<endl;
  Double_t ch_ndf=chi/ndf;

  //================= psi2S signal estimation =================
  double m_psi2 = fitFctCB22->GetParameter(5);
  double m_psi2_er = fitFctCB22->GetParError(5);
  double s_psi2 = fitFctCB22->GetParameter(6);
  double s_psi2_er = fitFctCB22->GetParError(6);
  double npsi2s3s = fitFctCB22->Integral(m_psi2-3*s_psi2,m_psi2+3*s_psi2)/hT->GetBinWidth(1);
  double nerrpsi2s3s = fitFctCB22->IntegralError(m_psi2-3*s_psi2,m_psi2+3*s_psi2,fitFctCB22->GetParameters(),covPeak2.GetMatrixArray())/hT->GetBinWidth(1);
  //double nerrpsi2s3s = fitFctCB22->Integral(3.55,3.8);
  cout<<" nerrpsi2s3s "<<nerrpsi2s3s<<endl;
  cout<< "No of Psi(2S) == " <<  npsi2s3s <<"+/-" << nerrpsi2s3s<<endl;
  //===============================================================================
  //Computation of bin significance and signal over background
  
  double nbck3s = exp2->Integral(m-3.*s,m+3.*s)/hT->GetBinWidth(1);
  //double nbck3sErr = exp2->IntegralError(m-3*s,m+3*s,r2->GetParams(),r2->GetCovarianceMatrix().GetMatrixArray() )/hT->GetBinWidth(1);
  double nbck3sErr = exp2->IntegralError(m-3*s,m+3*s,exp2->GetParameters(),covBck.GetMatrixArray())/hT->GetBinWidth(1);
  double sOverB3s = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));

  cout<<" S/B ==="<< sOverB3s << "+/-" << sOverB3sErr <<endl;
  //Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);
  
  double sig = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                               TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );


  cout<<" s/sqrt(s+b) " <<sig << " +/- "<<sigErr<<endl;  

  //==========================================================================================================

  
  TText *text0 = display1->AddText(Form("#chi^{2}/ndf =%f",ch_ndf));
  TText *text1 = display1->AddText(Form("M_{J/#psi}(3#sigma) =%f +/- %f GeV",m,m_err));
  TText *text2 = display1->AddText(Form("#sigma_{J/#psi}(3#sigma) = %f +/- %f GeV",s,s_err));
  TText *text3 = display2->AddText(Form("N_{J/#psi}(3#sigma) = %f +/- %f", njpsi3s, nerr3s));
  TText *text4 = display2->AddText(Form("N_{#psi2S}(3#sigma) = %f +/- %f", npsi2s3s, nerrpsi2s3s));
  TText *text5 = display2->AddText(Form("#frac{S}{B}_{J/#psi}(3#sigma) = %f +/- %f",sOverB3s,sOverB3sErr));
  //TText *text5 = display2->AddText(Form("#frac{S}{#sqrt{S+B}}_{J/#psi}(3#sigma) = %f +/- %f",sig,sigErr));



  display1->Draw("same");
  display2->Draw("same");
  
  //TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
  //leg->AddEntry("fitFctCB2","J/#psi Signal","l");
  //leg->AddEntry("fitFctCB22","#psi Signal","l");
  //leg->AddEntry("exp2","Background","l");
  //leg->Draw();

  

  c->Update();
  c->SaveAs("Spectra2.root");  
  
 
}

//------------------------------------------------------------------------------
Double_t BackgroundPol3(Double_t *x, Double_t *par)
{
  // pol2 3 params
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0];
}


//------------------------------------------------------------------------------
Double_t BackgroundVWG(Double_t *x, Double_t *par)
{
  // gaussian variable width
  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
  return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
  
}

//------------------------------------------------------------------------------
Double_t CrystalBallExtended(Double_t *x,Double_t *par)
{
  //par[0] = Normalization
  //par[1] = mean
  //par[2] = sigma
  //par[3] = alpha
  //par[4] = n
  //par[5] = alpha'
  //par[6] = n'  
  
  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  
  if (t >= -absAlpha && t < absAlpha2) // gaussian core
    {
      return par[0]*(exp(-0.5*t*t));
    }
  
  if (t < -absAlpha) //left tail
    {
      Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
      Double_t b = par[4]/absAlpha - absAlpha;
      return par[0]*(a/TMath::Power(b - t, par[4]));
    }
  
  if (t >= absAlpha2) //right tail
    {
    
      Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
      Double_t d = par[6]/absAlpha2 - absAlpha2;
      return par[0]*(c/TMath::Power(d + t, par[6]));
    }
  
  return 0. ; 
} 

//------------------------------------------------------------------------------

Double_t CrystalBallExtended2(Double_t *x,Double_t *par)
{
  //par[0] = Normalization
  //par[1] = mean
  //par[2] = sigma
  //par[3] = alpha
  //par[4] = n
  //par[5] = alpha'
  //par[6] = n'  
  
  Double_t t = (x[0]-par[8])/par[9];
  if (par[10] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[10]);
  Double_t absAlpha2 = fabs((Double_t)par[12]);
  
  if (t >= -absAlpha && t < absAlpha2) // gaussian core
    {
      return par[7]*(exp(-0.5*t*t));
    }
  
  if (t < -absAlpha) //left tail
    {
      Double_t a =  TMath::Power(par[11]/absAlpha,par[11])*exp(-0.5*absAlpha*absAlpha);
      Double_t b = par[11]/absAlpha - absAlpha;
      return par[7]*(a/TMath::Power(b - t, par[11]));
    }
  
  if (t >= absAlpha2) //right tail
    {
    
      Double_t c =  TMath::Power(par[13]/absAlpha2,par[13])*exp(-0.5*absAlpha2*absAlpha2);
      Double_t d = par[13]/absAlpha2 - absAlpha2;
      return par[7]*(c/TMath::Power(d + t, par[13]));
    }
  
  return 0. ; 
} 

//------------------------------------------------------------------------------



Double_t fitFunctionVWG(Double_t *x, Double_t *par)
{
  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  return BackgroundVWG(x, par);
}

Double_t fitFunctionPol3(Double_t *x, Double_t *par)
{
  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  return BackgroundPol3(x, par);
}

//------------------------------------------------------------------------------
Double_t fitFunctionCB2VWG(Double_t *x, Double_t *par)
{
  return BackgroundVWG(x, par) + CrystalBallExtended(x, &par[4]);
}

//--------------------------------------------------------------------------------
Double_t FuncBck1(Double_t *x, Double_t *par)
{ //exponential
  Double_t FitBck1 = exp(par[0]+par[1]*x[0]);
  return FitBck1;
}


Double_t FuncBck2(Double_t *x, Double_t *par)
{ //exponential
  Double_t FitBck2 = exp(par[2]+par[3]*x[0]);
  return FitBck2;
}

Double_t BackgroundBCK(Double_t *x, Double_t *par)
{ //exponential
  return FuncBck1(x,par)+FuncBck2(x,par);
}
//.......................................................................................
Double_t fitFunctionBCK(Double_t *x, Double_t *par)
{
  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();                      
  if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();                       
  return BackgroundBCK(x, par);                                           
}                                                                         
//................................................................        
                                                              

Double_t fitFunction2CB2VWG(Double_t *x, Double_t *par)
{
  Double_t par2[7] = { par[11],
		       //3.68609+(par[5]-3.096916)/3.096916*3.68609,
		       //par[6]/3.096916*3.68609,
		       par[5]+(3.68609-3.096916),
		       par[6]*(3.68609/3.096916),
		       par[7],
		       par[8],
		       par[9],
		       par[10] };
  return fitFunctionVWG(x, par) + CrystalBallExtended(x, &par[4]) + CrystalBallExtended(x, par2);
}

//------------------------------------------------------------------------------
Double_t fitFunction2CB2Pol3(Double_t *x, Double_t *par)
{
  Double_t par2[7] = { par[11],
		       //3.68609+(par[5]-3.096916)/3.096916*3.68609,
		       //par[6]/3.096916*3.68609,
		       par[5]+(3.68609-3.096916),
		       par[6]*(3.68609/3.096916),
		       par[7],
		       par[8],
		       par[9],
		       par[10] };
  return fitFunctionPol3(x, par) + CrystalBallExtended(x, &par[4]) + CrystalBallExtended(x, par2);
}



void SetStyle(Bool_t graypalette) {
  //cout << "Setting style!" << endl;
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(1);
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
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kBlack);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.019,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(1.31,"y");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(28);
  gStyle->SetTextFont(42);
  gStyle->SetTickLength(0.03,"X");
  gStyle->SetTickLength(0.03,"Y"); 
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
}
