#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TString.h"
#include <TStopwatch.h>
#include <TComplex.h>
#include <vector>
#include "../include/toyflowinputs.h"
#include "../include/rootcommon.h"

using namespace std;

Int_t gMarkers[]= {20,24,21,25,22,26,23,27,32,28};
Int_t gColors[]={kRed+1, kOrange+2, kCyan+2, kSpring-6, kRed-7, kOrange+1,kCyan-6,kGreen+7,kRed-9,kOrange-9,kAzure+6,kGreen-9};
Int_t gStyles[]={1,2,3,4,5,6,7,8,9,10};

TH1D *hDeltaPhiSum[NC];
Double_t vn[NH][NC]={-999};
Double_t vnError[NH][NC]={-999};
void LoadData();
void FitData();
void DrawTwoPArt();

void FitTwoPart(TString infile="output.root")
{
	LoadData();
	//Functions below should be in a for loop over NC
	FitDrawTwo();
}
//-------Member functions------------
void LoadData(TString inputname)
{
	TFile *fIn = TFile::Open(inputname,"read");
	for(int ic=0; ic<NC; ic++){
		hDeltaPhiSum[ic]=(TH1D*)fIn->Get(Form("hDeltaPhiSum_C%02d",ic));

	}

}

void FitDrawTwo(int ic=0)
{
	TF1 *fFit = new TF1("fFit", "[0]*(1+2*TMath::Power([1],2)*TMath::Cos(1*x) + 2*TMath::Power([2],2)*TMath::Cos(2*x) + 2*TMath::Power([3],2)*TMath::Cos(3*x) + 2*TMath::Power([4],2)*TMath::Cos(4*x) + 2*TMath::Power([5],2)*TMath::Cos(5*x))", 0, 2.0*TMath::Pi());
	fFit->SetParameter(0,1E4);
	fFit->SetParameter(1,0);
	fFit->SetParameter(2,0.10);
	fFit->SetParameter(3,0.06);
	fFit->SetParameter(4,0.06);
	fFit->SetParameter(5,0.06);
	fFit->SetParNames("const","v1","v_2","v_3","v_4","v_5");
	TF1 *fFitvn[NH];

	gStyle->SetOptStat(0);

	TCanvas *can = new TCanvas("C","canvas",1024,740);
	can->SetFillStyle(4000);
	//fit
	hDeltaPhiSum[ic]->Fit("fFit"); 

	double lowx = 0.,highx=2*TMath::Pi();
  	double ly=hDeltaPhiSum[ic]->GetMinimum()*0.99,hy=hDeltaPhiSum[ic]->GetMaximum()*1.01;
  	TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
  	hset( *hfr, "#Delta#phi=#phi_{1}-#phi_{2}", "dN/d#Delta#phi",0.9,0.9, 0.05,0.05, 0.01,0.01, 0.03,0.03, 510,505);//settings of the upper pad: x-axis, y-axis
  	hfr->Draw();
  	
	hDeltaPhiSum[ic]->SetMarkerStyle(20);
	hDeltaPhiSum[ic]->Draw("psame");
	fFit->SetLineColor(1);
	fFit->Draw("same");

	/*fourier->SetLineColor(kRed);
	fourier->SetLineStyle(10);
	fourier->Draw("same");*/
    
    for (Int_t n=0; n<=(NH-1); n++){
    	TString formula = Form("[0]*(1 + 2*TMath::Power([1],2)*TMath::Cos(%d*x))",n+1);
		fFitvn[n]= new TF1(Form("fFitvn%02d",n+1),formula, 0, 2.0*TMath::Pi());//Jus for drawing
	}

	//get vn's from fit

    TLegend *legendPhi = new TLegend(0.55,0.7,0.75,0.9,"","brNDC");
    legendPhi->SetTextSize(0.03);legendPhi->SetBorderSize(0);legendPhi->SetFillStyle(0);//legend settings;
	//for loop 
	for (Int_t n=1; n<=(NH-1); n++) {
		vn[n][ic]=fFit->GetParameter(n+1);
		vnError[n][ic]=fFit->GetParError(n+1);
		fFitvn[n]->SetParameter(1,vn[n][ic]);// Setting individual component vn
		fFitvn[n]->SetParameter(0, fFit->GetParameter(0));//Normalization
		fFitvn[n]->SetLineColor(gColors[n+1]);
		fFitvn[n]->Draw("same");
		legendPhi->AddEntry(fFitvn[n],Form("n = %d, v_n = %.3f #pm %.4f ",n+1, vn[n][ic], vnError[n][ic]),"l");
	}	
	legendPhi->Draw();
	gPad->GetCanvas()->SaveAs(Form("figs/TwoPartDecomposeC%02d.pdf",ic));

}

void SaveVns()
{
	double Cent[NC] = {0.,1.,2.};
	double eCent[NC] = {0.,0.,0.};
	TGraphErrors *gr_fittwo[NH];
	for(int ih=0; ih<NH; ih++)gr_fittwo[ih]= new TGraphErrors(NC,Cent,vn[ih],eCent,vnError[ih]);
	TFile *output = new TFile("out_VnFitTwo.root");
	for(int ih=0; ih<NH; ih++)gr_fittwo[ih]->Write(Form("gr_fittwo_H%02d_cent",ih));
	output->Write();
	output->Close();
}

