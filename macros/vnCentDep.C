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

Int_t gMarkers[]= {20,25,22,26,23,27,32,28};
Int_t gColors[]={kRed+1, kCyan+2, kSpring-6, kRed-7, kOrange+1,kCyan-6,kGreen+7,kRed-9,kOrange-9,kAzure+6,kGreen-9};
Int_t gStyles[]={1,2,3,4,5,6,7,8,9,10};
const int NMethod=3;
TGraphErrors *gr_vn_cent[NMethod][NH]; // Will need to add to array after reading in from all inputfiles??
TString gr_Names[NMethod]={"SP","TP","EP"};
TGraphErrors *gr_vnFit_cent[2][NH]; // Will need to add to array after reading in from all inputfiles??
TString gr_FitNames[2]={"single","two"};

void LoadData(TString/*,TString,TString*/); //Loading TGraphs
void DrawVnCent(int);

//---Main Function------
void vnCentDep(TString infile="../output/toymcflowao_1d0eb42_10k.root")
{
	LoadData(infile);
	for (int ih=1; ih<NH; ih++) DrawVnCent(ih);
}

//-----Member Functions-------
void LoadData(TString inputname)
{
	TFile *fIn = TFile::Open(inputname,"read");
	for(int i=0; i<NMethod; i++){
		for (int ih=0; ih<NH; ih++){
			gr_vn_cent[i][ih]=(TGraphErrors*)fIn->Get(Form("gr_v%02d_%s_cent",ih+1, gr_Names[i].Data()));
		}
	}
	// Reading in extra files from fit
	TString extrafiles[2] = {"out_VnFitSingle.root","out_VnFitTwo.root"};
	TFile *fext[2];
	for(int i=0; i<2; i++){
		fext[i]=TFile::Open(extrafiles[i]);
	}
	for(int i=0; i<2; i++){
		for (int ih=0; ih<NH; ih++){
			cout<< Form("gr_fit%s_H%02d_cent", gr_FitNames[i].Data(),ih+1)<<endl;
			gr_vnFit_cent[i][ih]=(TGraphErrors*)fext[i]->Get(Form("gr_fit%s_H%02d_cent", gr_FitNames[i].Data(),ih+1));
			gr_vnFit_cent[i][ih]->Print();
		}
	}
}

void DrawVnCent(int ih=1)
{	
	gStyle->SetOptStat(0);
	TCanvas *can = new TCanvas("C","canvas",1024,740);
	can->SetLeftMargin(0.15);
	can->SetBottomMargin(0.15);
	can->SetFillStyle(4000);
	TLegend *legend = new TLegend(0.45,0.3,0.65,0.5,"","brNDC");
	legend->SetTextSize(0.04);legend->SetBorderSize(0);legend->SetFillStyle(0);//legend settings;
	double lowx = -0.5,highx=2.5;
	double ly=-0.01,hy=inputVn[ih][NC-2]*1.8;
	TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
	hset( *hfr, "Centrality [%]", "v_{n}",0.7,0.7, 0.07,0.07, 0.01,0.01, 0.03,0.03, 510,505);//settings of the upper pad: x-axis, y-axis
	hfr->Draw();

	//legend->AddEntry((TObjArray*)NULL,Form("Centrality %s",strCentrality[ic].Data()));

	for(int i=0; i<NMethod-1; i++){
			gr_vn_cent[i][ih]->SetLineColor(gColors[i]);
			gr_vn_cent[i][ih]->SetMarkerStyle(gMarkers[i]);
			gr_vn_cent[i][ih]->SetMarkerColor(gColors[i]);
			gr_vn_cent[i][ih]->Draw("psame");
			legend->AddEntry(gr_vn_cent[i][ih],Form("n=%d %s", ih+1, gr_Names[i].Data()), "p");
	}
	for(int i=0; i<2; i++){
			//cout << Form("n=%02d, FitMethod=%02d", ih, i)<<endl;
			gr_vnFit_cent[i][ih]->SetLineColor(gColors[i+NMethod]);
			gr_vnFit_cent[i][ih]->SetMarkerStyle(gMarkers[i+NMethod]);
			gr_vnFit_cent[i][ih]->SetMarkerColor(gColors[i+NMethod]);
			gr_vnFit_cent[i][ih]->Draw("psame");
			legend->AddEntry(gr_vnFit_cent[i][ih],Form("n=%d %s", ih+1, gr_FitNames[i].Data()),"p");
	}
	legend -> Draw("same");
	gPad->GetCanvas()->SaveAs(Form("../figs/CentDep_H%02d.pdf",ih+1));
}
