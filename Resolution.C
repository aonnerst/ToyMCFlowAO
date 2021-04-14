#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TFile.h"
#include "TRatioPlot.h"
#include "TString.h"
#include "include/rootcommon.h"

Int_t gMarkers[]= {20,24,21,25,22,26,23,27,32,28};
Int_t gColors[]={kRed+1, kOrange+2, kCyan+2, kSpring-6, kRed-7, kOrange+1,kCyan-6,kGreen+7,kRed-9,kOrange-9,kAzure+6,kGreen-9};
Int_t gStyles[]={1,2,3,4,5,6,7,8,9,10};

void Resolution()
{
	TFile *fIn = TFile::Open("outputDJ.root","read");
    const Int_t NH = 5;

    TH1D *hResolutionDist[NH];

    gStyle->SetOptStat(0);
    TCanvas *cRes = new TCanvas("cRes","cRes");

    for (Int_t i=0; i<NH; i++){
    	hResolutionDist[i]=(TH1D*)fIn->Get(Form("hResolutionDist%02d",i+1));
    	//hResolutionDist[i]->Sumw2();
	   	//hResolutionDist[i] -> SetLineColor(i+1);
    }

    double lowx = -10,highx=10;
  	double ly=hResolutionDist[0]->GetMinimum()*0.1,hy=hResolutionDist[0]->GetMaximum()*1.4;
  	TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
  	hset( *hfr, "#psi_{n}-#psi_{EP}", "dN/d#phi",0.7,0.7, 0.07,0.07, 0.01,0.01, 0.03,0.03, 510,505);//settings of the upper pad: x-axis, y-axis
  	hfr->Draw();

  	/*fourier->SetLineColor(kBlack);
   	fourier->SetLineStyle(10);
   	fourier->Draw("same");*/


  	TLegend *legendRes = new TLegend(0.5,0.6,0.8,0.85,"","brNDC");
    legendRes->SetTextSize(0.04);legendRes->SetBorderSize(0);legendRes->SetFillStyle(0);//legend settings;

    for(Int_t i=0; i<NH; i++){
   		hResolutionDist[i]->Draw("lpsame");
   		hResolutionDist[i]->SetLineColor(gColors[i]);
   		hResolutionDist[i]->SetLineWidth(5);
   		//hResolutionDist[i]->SetLineStyle(gStyles[i]);
   		//hResolutionDist[i]->SetMarkerStyle(gMarkers[i]);
   		//hResolutionDist[i]->SetMarkerColor(gColors[i]);
   		legendRes->AddEntry(hResolutionDist[i],Form("n=%d",i+1));
   	}
   	legendRes -> Draw("same");
   	gPad->GetCanvas()->SaveAs("figs/Resolution.pdf");

}