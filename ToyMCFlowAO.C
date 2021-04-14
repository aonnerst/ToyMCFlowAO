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

std::string prd(const double x, const int decDigits, const int width);
std::string center(const string s, const int w);
double DeltaPhi(double phi1, double phi2); // relative angle

void toymc()
{
	// Declare variables
	 //Number of Centralities
	Int_t Nevt = 10000;
	
	Int_t NPhiHist = 10;
	Double_t phiarray[Nch];
	Double_t Psi_n[NH]={0.0};
	Double_t vn_psi[NC][NH]={0.0};
	Double_t vn_obs_EP[NC][NH]={0.0};
	Double_t vn_phi[NC][NH]={0.0};
	Double_t weight = 1.0;
	//Double_t Qn_x[NH];//={0.0};
	//Double_t Qn_y[NH];//={0.0};
	Double_t MeanArrayTwoParticle[NC][NH]={0.0};
	Double_t MeanArrayEventPlane[NC][NH]={0.0};
	Double_t MeanArrayEventPlaneQVec[NC][NH]={0.0};
	Double_t MeanArrayEvtPlError[NC][NH]={0.0};
	Double_t MeanArrayEvtPlErrorQvec[NC][NH]={0.0};
	Double_t MeanArrayTwoPartError[NC][NH]={0.0};
	Double_t MeanArrayResolution[NC][NH]={0.0};
	Double_t MeanArrayResolutionError[NC][NH]={0.0};
	Double_t vn_obs_ERROR[[NC]NH]={0.0};
	Double_t vn_TwoPart[NC][NH]={0.0};
	Double_t vn_EvtPl[NC][NH]={0.0};
	Double_t vn_EvtPlQvec[NC][NH]={0.0};



	
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

	TH1D *hEventPlane[NH];
	TH1D *hEventPlaneEP[NH];
	TH1D *hTwoParticle[NH];
	TH1D *hPhiPsi[NH];
	TH1D *hPhiPsiQ[NH];
	TH1D *hDeltaPhi[NH];
	TH1D *hPhiEvent[NPhiHist];
	TH1D *hResolution[NH];
	TH1D *hResolutionDist[NH];

	TF1 *fourier = new TF1("Fourier", "[0]*(1+2*[1]*TMath::Cos(1*(x-[6])) + 2*[2]*TMath::Cos(2*(x-[7])) + 2*[3]*TMath::Cos(3*(x-[8])) + 2*[4]*TMath::Cos(4*(x-[9])) + 2*[5]*TMath::Cos(5*(x-[10])))", 0.0, 2.0*TMath::Pi());
	//range 0 to 2*pi
	for (Int_t n=0; n<NH; n++){

		//------Symmetry planes random ------
		uniform[n]= new TF1(Form("uniform%02d",n+1),"[0]", -1*TMath::Pi()/(n+1), 1.0*TMath::Pi()/(n+1));
		uniform[n]->SetParameter(0,1.);

		//-----Histograms---------
		hEventPlane[n]=new TH1D(Form("hEventPlane%02d",n+1),Form("hEventPlane%02d",n+1),200,-1.0, 1.0);
		hEventPlaneEP[n]=new TH1D(Form("hEventPlaneEP%02d",n+1),Form("hEventPlaneEP%02d",n+1),200,-1.0, 1.0);
		hTwoParticle[n]=new TH1D(Form("hTwoParticle%02d",n+1),Form("hTwoParticle%02d",n+1),200,-1.0, 1.0);
		hPhiPsi[n] = new TH1D(Form("hPhiPsi%02d",n+1),Form("hPhiPsi%02d",n+1),200,0.0, 2.0*TMath::Pi());
		hPhiPsiQ[n] = new TH1D(Form("hPhiPsiQ%02d",n+1),Form("hPhiPsiQ%02d",n+1),200,0.0, 2.0*TMath::Pi());
		hDeltaPhi[n] = new TH1D(Form("hDeltaPhi%02d",n+1),Form("hDeltaPhi%02d",n+1),200,0.0, 2.0*TMath::Pi());
		hResolution[n] = new TH1D(Form("hResolution%02d",n+1),Form("hResolution%02d",n+1),200,-100, 100);
		hResolutionDist[n] = new TH1D(Form("hResolutionDist%02d",n+1),Form("hResolutionDist%02d",n+1),200,-10, 10);
	}

	//decided to first define option B, hence need to be finished
	TH1D *hDeltaPhiSum[NC] = new TH1D(Form("hDeltaPhiSumCent%d",),"hDeltaPhiSum",200, 0.0, 2.0*TMath::Pi());

	for (Int_t iPhiEvt=0; iPhiEvt<NPhiHist; iPhiEvt++){
		hPhiEvent[iPhiEvt] = new TH1D(Form("hPhiEvent%02d",(iPhiEvt+1)),Form("hPhiEvent%02d",(iPhiEvt+1)),100,0.0, 2.0*TMath::Pi());
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
    	if(dice>= 0.0 && dice<0.3){ ic=0;}
    	if else(dice>=0.3 && dice<0.6){ic=1;}
    	if else(dice){ic=2;}
    	Nch=inputNch[ic];
		//Get Psi for different harmonics
		for (Int_t n=0; n<NH; n++) Psi_n[n]=uniform[n]->GetRandom();//harmonic loop
		fourier->SetParameter(0,Nch); 
		for (Int_t i=0; i<NH-1; i++)fourier->SetParameter(i+1,inputVn[NH][ic]); //Setting the vn parameters
		for (Int_t i=NH; i<2*NH; i++)fourier->SetParameter(i+1,Psi_n[i-NH]); //Setting the Psi parameters

		//Initializing 
		phiarray[Nch]={-999.0};
		Double_t Qn_x[NH] = {0.0};//Should be 0 since we sum the qvectors
		Double_t Qn_y[NH] = {0.0};//
		TComplex QvectorsEP[NH];
		Double_t Psi_n_EP[NH]={0.0};
		Double_t Psi_n_EPQ[NH]={0.0};
		Double_t Resolution[NH]={0.0};
    	for(int iH=0;iH<NH;iH++) QvectorsEP[iH] = TComplex(0,0);

    	
	
		for (Int_t t=0; t<Nch; t++)//track loop
		{
			phiarray[t] = fourier->GetRandom();

			if(iEvent<NPhiHist) {
				hPhiEvent[iEvent]->Fill(phiarray[t]);
				fourier->Write(Form("fourier%02d",iEvent));
			}
			
			//Harmonic loop
			for (Int_t n=0; n<NH; n++)
			{
				hPhiPsi[n]->Fill(DeltaPhi(phiarray[t],Psi_n[n]));
				vn_psi[n] = TMath::Cos((n+1)*(DeltaPhi(phiarray[t], Psi_n[n]))); //Change to symmetry plane
				hEventPlane[n]->Fill(vn_psi[n]);
				
				// calculating eventplane with Q-vectors
				Qn_x[n] += weight*TMath::Cos((n+1)*phiarray[t]);
				Qn_y[n] += weight*TMath::Sin((n+1)*phiarray[t]);
				QvectorsEP[n] += TComplex(TMath::Cos((n+1)*phiarray[t]),TMath::Sin((n+1)*phiarray[t]));
			}

		}//End of track loop

		//Only after the track loop, must sum over the tracks first
		for (Int_t n=0; n<NH; n++) Psi_n_EP[n]=(1/double(n+1))*TMath::ATan2(Qn_y[n],Qn_x[n]); 
		for (Int_t n=0; n<NH; n++) Psi_n_EPQ[n] = QvectorsEP[n].Theta()/double(n+1);

		//Two Particle correlation
		for (Int_t i=0; i<Nch; i++){
			//Evenplane method calculated vn
			for (Int_t n=0; n<NH; n++) {
				hEventPlaneEP[n]->Fill(TMath::Cos((n+1)*(DeltaPhi(phiarray[i], Psi_n_EP[n])))); 
				hPhiPsiQ[n]->Fill(DeltaPhi(phiarray[i], Psi_n_EP[n]));
			}
			for (Int_t j=0; j<Nch;j++){
				if(i==j) continue;
				hDeltaPhiSum->Fill(DeltaPhi(phiarray[i], phiarray[j]));//For fitting
				for (Int_t n=0; n<NH; n++){
					vn_phi[n] = TMath::Cos((n+1)*(DeltaPhi(phiarray[i], phiarray[j])));//Analytic solution
					hTwoParticle[n]->Fill(vn_phi[n]);

				}
			}		
		}
		
		//Resolution for every event
		for (Int_t n=0; n<NH; n++)
		{
			Resolution[n] = TMath::Cos((n+1)*(DeltaPhi(Psi_n[n], Psi_n_EP[n]))); //Analystical resultion
			//if(n==1) cout <<  Form("n=%d,psi=%.3f, %.3f, %.3f",n,Psi_n[n], Psi_n_EP[n], Psi_n_EPQ[n]) << endl;
			hResolution[n]->Fill(Resolution[n]);
			hResolutionDist[n]->Fill(DeltaPhi(Psi_n[n], Psi_n_EP[n]));
		}
		
	}// End of event loop



	// Calculating the avarage over event
	for (Int_t n=0; n<NH; n++)
	{
		MeanArrayTwoParticle[n]=hTwoParticle[n]->GetMean();
		MeanArrayTwoPartError[n]=hTwoParticle[n]->GetMeanError()/((n+1)*MeanArrayTwoParticle[n]);
		MeanArrayEventPlane[n]=hEventPlane[n]->GetMean();
		MeanArrayEvtPlError[n]=hEventPlane[n]->GetMeanError();
		MeanArrayEventPlaneQVec[n]=hEventPlaneEP[n]->GetMean();
		MeanArrayEvtPlErrorQvec[n]=hEventPlaneEP[n]->GetMeanError();
		MeanArrayResolution[n]=hResolution[n]->GetMean();
		MeanArrayResolutionError[n]=hResolution[n]->GetMeanError();



		vn_TwoPart[n]=TMath::Sqrt(TMath::Abs(MeanArrayTwoParticle[n]));
		vn_EvtPl[n]=MeanArrayEventPlane[n];
		vn_EvtPlQvec[n]=MeanArrayEventPlaneQVec[n]/MeanArrayResolution[n];
		vn_obs_ERROR[n]=TMath::Power(((1/MeanArrayResolution[n])*MeanArrayEvtPlErrorQvec[n]),2)+TMath::Power(((MeanArrayEventPlaneQVec[n]/TMath::Power(MeanArrayResolution[n],2))*MeanArrayResolutionError[n]),2);
	}

	
	output->Write();
	output->Close();

	timer.Print();

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
              		<< prd(/*vn_fit[n]*/0,5,10)       << " | "
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
              		<< prd(/*vn_fit[n]*/0,5,10)       << " | "
              		<< prd(MeanArrayResolution[n],5,10)   <<  " +-"<<prd(MeanArrayResolutionError[n],5,9)  <<"\n";
		}

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