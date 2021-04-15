#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TFile.h"
#include "TRatioPlot.h"
#include "TString.h"
#include <TStopwatch.h>
#include <TComplex.h>
#include "toyflowinputs.h"

std::string prd(const double x, const int decDigits, const int width);
std::string center(const string s, const int w);
double DeltaPhi(double phi1, double phi2); // relative angle

void toymc()
{
	// Declare variables
	cout<< strCentrality[0]<<endl;
	Int_t Nevt = 10000;
	Int_t NPhiHist = 10;
	Double_t Psi_n[NH]={0.0};
	Double_t vn_psi[NH][NC]={0.0};
	Double_t vn_obs_EP[NH][NC]={0.0};
	Double_t vn_phi[NH][NC]={0.0};
	Double_t weight = 1.0;
	//Double_t Qn_x[NH];//={0.0};
	//Double_t Qn_y[NH];//={0.0};
	Double_t MeanArrayTwoParticle[NH][NC]={0.0};
	Double_t MeanArrayEventPlane[NH][NC]={0.0};
	Double_t MeanArrayEventPlaneQVec[NH][NC]={0.0};
	Double_t MeanArrayEvtPlError[NH][NC]={0.0};
	Double_t MeanArrayEvtPlErrorQvec[NH][NC]={0.0};
	Double_t MeanArrayTwoPartError[NH][NC]={0.0};
	Double_t MeanArrayResolution[NH][NC]={0.0};
	Double_t MeanArrayResolutionError[NH][NC]={0.0};
	Double_t vn_obs_ERROR[NH][NC]={0.0};
	Double_t vn_TwoPart[NH][NC]={0.0};
	Double_t vn_EvtPl[NH][NC]={0.0};
	Double_t vn_EvtPlQvec[NH][NC]={0.0};



	
	TString outfile = "output.root";
	TFile *output = new TFile(outfile,"recreate");
	output->cd();
	
	//Define uniform function for option B
	TF1 *centSamp = new TF1("centSamp", "[0]",0.0,0.9);
	centSamp->SetParameter(0,1.0);

	//Define histogram for option B
	TH1D *hCentSample = new TH1D("hCentSample","hCentSample",3,0.0,0.9);

	TF1 *uniform[NH];
	//TF1 *fourier[NH]; 

	TH1D *hEventPlane[NH][NC];
	TH1D *hEventPlaneEP[NH][NC];
	TH1D *hTwoParticle[NH][NC];
	TH1D *hPhiPsi[NH][NC];
	TH1D *hPhiPsiQ[NH][NC];
	TH1D *hDeltaPhi[NH][NC];
	TH1D *hPhiEvent[NPhiHist][NC];
	TH1D *hResolution[NH][NC];
	TH1D *hResolutionDist[NH][NC];

	TF1 *fourier = new TF1("Fourier", "[0]*(1+2*[1]*TMath::Cos(1*(x-[6])) + 2*[2]*TMath::Cos(2*(x-[7])) + 2*[3]*TMath::Cos(3*(x-[8])) + 2*[4]*TMath::Cos(4*(x-[9])) + 2*[5]*TMath::Cos(5*(x-[10])))", 0.0, 2.0*TMath::Pi());
	//range 0 to 2*pi
	for (Int_t ih=0; ih<NH; ih++){

		//------Symmetry planes random ------
		uniform[ih]= new TF1(Form("uniform%02d",ih+1),"[0]", -1*TMath::Pi()/(ih+1), 1.0*TMath::Pi()/(ih+1));
		uniform[ih]->SetParameter(0,1.);

		for (Int_t ic=0; ic<NC; ic++){
			//-----Histograms---------
			hEventPlane[ih][ic]=new TH1D(Form("hEventPlane%02d_%s",ih+1,strCentrality[ic].Data()),Form("hEventPlane%02d_%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
			hEventPlaneEP[ih][ic]=new TH1D(Form("hEventPlaneEP%02d_%s",ih+1,strCentrality[ic].Data()),Form("hEventPlaneEP%02d_%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
			hTwoParticle[ih][ic]=new TH1D(Form("hTwoParticle%02d_%s",ih+1,strCentrality[ic].Data()),Form("hTwoParticle%02d_%s",ih+1,strCentrality[ic].Data()),200,-1.0, 1.0);
			hPhiPsi[ih][ic] = new TH1D(Form("hPhiPsi%02d_%s",ih+1,strCentrality[ic].Data()),Form("hPhiPsi%02d_%s",ih+1,strCentrality[ic].Data()),200,0.0, 2.0*TMath::Pi());
			hPhiPsiQ[ih][ic] = new TH1D(Form("hPhiPsiQ%02d_%s",ih+1,strCentrality[ic].Data()),Form("hPhiPsiQ%02d_%s",ih+1,strCentrality[ic].Data()),200,0.0, 2.0*TMath::Pi());
			hDeltaPhi[ih][ic] = new TH1D(Form("hDeltaPhi%02d_%s",ih+1,strCentrality[ic].Data()),Form("hDeltaPhi%02d_%s",ih+1,strCentrality[ic].Data()),200,0.0, 2.0*TMath::Pi());
			hResolution[ih][ic] = new TH1D(Form("hResolution%02d_%s",ih+1,strCentrality[ic].Data()),Form("hResolution%02d_%s",ih+1,strCentrality[ic].Data()),200,-100, 100);
			hResolutionDist[ih][ic] = new TH1D(Form("hResolutionDist%02d_%s",ih+1,strCentrality[ic].Data()),Form("hResolutionDist%02d_%s",ih+1,strCentrality[ic].Data()),200,-10, 10);
		}
		
		
	}

	TH1D *hDeltaPhiSum[NC];
	for (Int_t ic=0; ic<NC; ic++){
		hDeltaPhiSum[ic] = new TH1D(Form("hDeltaPhiSum_%s",strCentrality[ic].Data()),Form("hDeltaPhiSum_%s",strCentrality[ic].Data()),200, 0.0, 2.0*TMath::Pi());
	}
	

	for (Int_t iPhiEvt=0; iPhiEvt<NPhiHist; iPhiEvt++){
		for (Int_t ic=0; ic<NC; ic++){
			hPhiEvent[iPhiEvt][ic] = new TH1D(Form("hPhiEvent%02d_%s",(iPhiEvt+1),strCentrality[ic].Data()),Form("hPhiEvent%02d_%s",(iPhiEvt+1),strCentrality[ic].Data()),100,0.0, 2.0*TMath::Pi());
		}
	}
	int ieout = Nevt/20;
    if (ieout<1) ieout=1;
    TStopwatch timer;
    timer.Start();
	//Event loop
	for (Int_t iEvent=0; iEvent<Nevt; iEvent++)
	{
		if(iEvent % ieout == 0) { cout << iEvent << "\t" << int(float(iEvent)/Nevt*100) << "%" << endl ;}
		//Sample randomly from centSamp
    	Double_t dice = centSamp->GetRandom();
    	Int_t Nch=0;
    	Int_t ic=0;
    	if(dice>= 0.0 && dice<0.3) ic=0;
    	if(dice>=0.3 && dice<0.6) ic=1;
    	if(dice>=0.6 && dice<=0.9) ic=2;
    	Nch=inputNch[ic];
		//Get Psi for different harmonics
		for (Int_t n=0; n<NH; n++) Psi_n[n]=uniform[n]->GetRandom();//harmonic loop
		fourier->SetParameter(0,Nch); 
		for (Int_t i=0; i<NH-1; i++)fourier->SetParameter(i+1,inputVn[NH][ic]); //Setting the vn parameters
		for (Int_t i=NH; i<2*NH; i++)fourier->SetParameter(i+1,Psi_n[i-NH]); //Setting the Psi parameters

		//Initializing 
		Double_t phiarray[Nch]={0.0};
		Double_t Qn_x[NH][NC] = {0.0};//Should be 0 since we sum the qvectors
		Double_t Qn_y[NH][NC] = {0.0};//
		TComplex QvectorsEP[NH][NC];
		Double_t Psi_n_EP[NH]={0.0};
		Double_t Psi_n_EPQ[NH]={0.0};
		Double_t Resolution[NH]={0.0};
    	for(int iH=0;iH<NH;iH++) QvectorsEP[iH][ic] = TComplex(0,0);

    	
	
		for (Int_t t=0; t<Nch; t++)//track loop
		{
			phiarray[t] = fourier->GetRandom();

			if(iEvent<NPhiHist) {
				hPhiEvent[iEvent][ic]->Fill(phiarray[t]);
				fourier->Write(Form("fourier%02d",iEvent));
			}
			
			//Harmonic loop
			for (Int_t n=0; n<NH; n++)
			{
				hPhiPsi[n][ic]->Fill(DeltaPhi(phiarray[t],Psi_n[n]));
				vn_psi[n][ic] = TMath::Cos((n+1)*(DeltaPhi(phiarray[t], Psi_n[n]))); //Change to symmetry plane
				hEventPlane[n][ic]->Fill(vn_psi[n][ic]);
				
				// calculating eventplane with Q-vectors
				Qn_x[n][ic] += weight*TMath::Cos((n+1)*phiarray[t]);
				Qn_y[n][ic] += weight*TMath::Sin((n+1)*phiarray[t]);
				QvectorsEP[n][ic] += TComplex(TMath::Cos((n+1)*phiarray[t]),TMath::Sin((n+1)*phiarray[t]));
			}

		}//End of track loop

		//Only after the track loop, must sum over the tracks first
		for (Int_t n=0; n<NH; n++) Psi_n_EP[n]=(1/double(n+1))*TMath::ATan2(Qn_y[n],Qn_x[n]); 
		for (Int_t n=0; n<NH; n++) Psi_n_EPQ[n] = QvectorsEP[n].Theta()/double(n+1);

		//Two Particle correlation
		for (Int_t i=0; i<Nch; i++){
			//Evenplane method calculated vn
			for (Int_t n=0; n<NH; n++) {
				hEventPlaneEP[n][ic]->Fill(TMath::Cos((n+1)*(DeltaPhi(phiarray[i], Psi_n_EP[n])))); 
				hPhiPsiQ[n][ic]->Fill(DeltaPhi(phiarray[i], Psi_n_EP[n]));
			}
			for (Int_t j=0; j<Nch;j++){
				if(i==j) continue;
				hDeltaPhiSum[ic]->Fill(DeltaPhi(phiarray[i], phiarray[j]));//For fitting
				for (Int_t n=0; n<NH; n++){
					vn_phi[n][ic] = TMath::Cos((n+1)*(DeltaPhi(phiarray[i], phiarray[j])));//Analytic solution
					hTwoParticle[n][ic]->Fill(vn_phi[n][ic]);

				}
			}		
		}
		
		//Resolution for every event
		for (Int_t n=0; n<NH; n++)
		{
			Resolution[n][ic] = TMath::Cos((n+1)*(DeltaPhi(Psi_n[n], Psi_n_EP[n]))); //Analystical resultion
			//if(n==1) cout <<  Form("n=%d,psi=%.3f, %.3f, %.3f",n,Psi_n[n], Psi_n_EP[n], Psi_n_EPQ[n]) << endl;
			hResolution[n][ic]->Fill(Resolution[n][ic]);
			hResolutionDist[n][ic]->Fill(DeltaPhi(Psi_n[n], Psi_n_EP[n]));
		}
		
	}// End of event loop



	// Calculating the avarage over event
	for (Int_t n=0; n<NH; n++)
	{
		for (Int_t ic=0; ic<NC; ic++){
			MeanArrayTwoParticle[n][ic]=hTwoParticle[n][ic]->GetMean();
			MeanArrayTwoPartError[n][ic]=hTwoParticle[n][ic]->GetMeanError()/((n+1)*MeanArrayTwoParticle[n][ic]);
			MeanArrayEventPlane[n][ic]=hEventPlane[n][ic]->GetMean();
			MeanArrayEvtPlError[n][ic]=hEventPlane[n][ic]->GetMeanError();
			MeanArrayEventPlaneQVec[n][ic]=hEventPlaneEP[n][ic]->GetMean();
			MeanArrayEvtPlErrorQvec[n][ic]=hEventPlaneEP[n][ic]->GetMeanError();
			MeanArrayResolution[n][ic]=hResolution[n][ic]->GetMean();
			MeanArrayResolutionError[n][ic]=hResolution[n][ic]->GetMeanError();



			vn_TwoPart[n][ic]=TMath::Sqrt(TMath::Abs(MeanArrayTwoParticle[n][ic]));
			vn_EvtPl[n][ic]=MeanArrayEventPlane[n][ic];
			vn_EvtPlQvec[n][ic]=MeanArrayEventPlaneQVec[n][ic]/MeanArrayResolution[n][ic];
			vn_obs_ERROR[n][ic]=TMath::Power(((1/MeanArrayResolution[n][ic])*MeanArrayEvtPlErrorQvec[n][ic]),2)+TMath::Power(((MeanArrayEventPlaneQVec[n][ic]/TMath::Power(MeanArrayResolution[n][ic],2))*MeanArrayResolutionError[n][ic]),2);	
		}

	}

	
	output->Write();
	output->Close();

	timer.Print();
/*
	std::cout << center("harmonic",10)   << " | "
          << center("input",10)     << " | "
          << center("vn{Psi_n}",20) << " | "
          << center("vn_obs/R",20)    << " | " 
          << center("vn{n}",20)     << " | "
          << center("fit",20)       << " | "
          << center("R",10)         << "\n";

	std::cout << std::string(10*7 + 2*7, '-') << "\n";

	for(Int_t n=1; n<NH; n++) 
		{
    		std::cout << prd(n+1,0,10)             << " | "
              		<< prd(vn[n],2,8)          << " | "
              		<< prd(vn_EvtPl[n],5,9)     << " +-"<<prd(MeanArrayEvtPlError[n],5,9)    << " | "
              		<< prd(vn_EvtPlQvec[n],5,9) << " +-"<<prd(vn_obs_ERROR[n],5,9)           << " | "
              		<< prd(vn_TwoPart[n],5,9)   << " +-"<<prd(MeanArrayTwoPartError[n],5,9)  << " | "
              	
              		<< prd(MeanArrayResolution[n],5,10)   <<  " +-"<<prd(MeanArrayResolutionError[n],5,9)  <<"\n";
		}	

	std::cout << center("harmonic",10)   << " | "
      	<< center("input",10)     << " | "
      	<< center("|input-vn{Psi_n}|-error",20) << " | "
      	<< center("|input-vn_obs/R|-error",20)    << " | " 
      	<< center("|input-vn{n}|-error",20)     << " | "
      	<< center("fit",20)       << " | "
      	<< center("R",10)         << "\n";

	std::cout << std::string(10*7 + 2*7, '-') << "\n";

	for(Int_t n=1; n<NH; n++) 
		{
    		std::cout << prd(n+1,0,10)         << " | "
              		<< prd(vn[n],2,8)          << " | "
              		<< prd((TMath::Abs(vn[n] - vn_EvtPl[n])-MeanArrayEvtPlError[n] ),5,9) << " | "
              		<< prd((TMath::Abs(vn[n] - vn_EvtPlQvec[n])-vn_obs_ERROR[n]),5,9) << " | "
              		<< prd((TMath::Abs(vn[n] - vn_TwoPart[n])-MeanArrayTwoPartError[n]),5,9) << " | "

              		<< prd(MeanArrayResolution[n],5,10)   <<  " +-"<<prd(MeanArrayResolutionError[n],5,9)  <<"\n";
		}
*/
}

///////////////////////////////////////////////////////

	/* Convert double to string with specified number of places after the decimal
   and left padding. */
std::string prd(const double x, const int decDigits, const int width) {
    stringstream ss;
    ss << fixed << right;
    ss.fill(' ');        // fill space around displayed #
    ss.width(width);     // set  width around displayed #
    ss.precision(decDigits); // set # places after decimal
    ss << x;
    return ss.str();
}

/*! Center-aligns string within a field of width w. Pads with blank spaces
    to enforce alignment. */
std::string center(const string s, const int w) {
    stringstream ss, spaces;
    int padding = w - s.size();                 // count excess room to pad
    for(int i=0; i<padding/2; ++i)
        spaces << " ";
    ss << spaces.str() << s << spaces.str();    // format with padding
    if(padding>0 && padding%2!=0)               // if odd #, add 1 space
        ss << " ";
    return ss.str();
}

double DeltaPhi(double phi1, double phi2) {
  // dphi
  double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
  return res>0 ? res : 2.*TMath::Pi()+res ;
}