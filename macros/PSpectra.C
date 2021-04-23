#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TGraphPainter.h"
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
const int NMethod=3;
TGraphErrors *gr_pvn[NC][NMethod];
TString gr_Names[NMethod]={"SP","TP","EP"};

void LoadData(TString); //Loading TGraphs
void DrawPSpectra(int);

//---Main Function------
void PSpectra(TString infile="../output/toymcflowao_1d0eb42_10k.root")
{
	LoadData(infile);
	for(int ic=0; ic<NC; ic++) DrawPSpectra(ic);
}

//------Member Functions-------
void LoadData(TString inputname)
{
	TFile *fIn = TFile::Open(inputname,"read");

	for(int i=0; i<NMethod; i++){
		for (int ic=0; ic<NC; ic++){
			gr_pvn[ic][i]=(TGraphErrors*)fIn->Get(Form("gr_pv%02d_%s",ic+1, gr_Names[i].Data()));
		}
	}
}

void DrawPSpectra(int ic=0)
{
	gStyle->SetOptStat(0);
	TCanvas *can = new TCanvas("C","canvas",1024,740);
	can->SetFillStyle(4000);
	can->SetLeftMargin(0.15);
   	can->SetBottomMargin(0.15);
   	can->SetLogy(1);
	TLegend *legend = new TLegend(0.5,0.6,0.8,0.85,"","brNDC");
    legend->SetTextSize(0.04);legend->SetBorderSize(0);legend->SetFillStyle(0);//legend settings;
	double lowx = -0.5,highx=7.5;
  	double ly=1.2e-5,hy=4e-1;
  	TH2F *hfr = new TH2F("hfr",Form("Cent%02d",ic), 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
  	hset( *hfr, "n", "v_{n}",0.7,0.7, 0.07,0.07, 0.01,0.01, 0.03,0.03, 510,505);//settings of the upper pad: x-axis, y-axis
  	hfr->Draw();

  	for(int i=0; i<NMethod; i++){
  	
		gr_pvn[ic][i]->SetLineColor(gColors[i]);
		gr_pvn[ic][i]->SetMarkerStyle(gMarkers[i]);
		gr_pvn[ic][i]->SetMarkerColor(gColors[i]);
		gr_pvn[ic][i]->Draw("psame");
		legend->AddEntry(gr_pvn[ic][i],Form("%s", gr_Names[i].Data()));
  		
  	}
  	legend -> Draw("same");
   	gPad->GetCanvas()->SaveAs(Form("../figs/PSpectraC%02d.pdf",ic));
}
