#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TNtuple.h"
#include "Math/Vector4D.h"
#include <math.h>
#include "TRandom.h"
#include <time.h>
#include <cstdio>
#include <fstream>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TH1D.h"
#include "RooGlobalFunc.h"
#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "TStyle.h"


using namespace RooFit;
using namespace ROOT::Math;
using namespace std;

//------------------------------------------------------------Function-to-choose-a-case-----------------------------------------------------------------------

int affiche() {

	int Case;

 	cout<<"Choisir le mode de fonctionnement : "<<endl;
 	cout<<"1. Etude d'une seul simulation à une certaine énergie."<<endl;
 	cout<<"2. Etude de plusieurs simulations à différentes énergies."<<endl;
 	cout<<"3. Etude des effets de la variation d'un paramétre de résolution."<<endl;
	cout<<"4. Détermination de l'incertitude sur la résolution à une certaine énergie"<<endl;
	do {
 	cin>>Case;
	} while (Case<1 || Case>4);
	
	return Case;
}


//-------------------------------------------------------Function to calculate the sum with the right resol---------------------------------------------------


template <class F> double calculSommeEnergie(F func,std::vector<ROOT::Math::XYZTVector> *vec, int m)         
{
	double sum=0;
	if (m==1) {                     //If Resolution on energy
		for (std::vector<ROOT::Math::XYZTVector>::iterator it=vec->begin(); it != vec->end(); ++it){
	    		sum+=func(it->energy());
		}
	}
	if (m==2) {                    //If Resolution on momentum
		for (std::vector<ROOT::Math::XYZTVector>::iterator it=vec->begin(); it != vec->end(); ++it){
			
	    		sum+=sqrt((it->mass2())+pow(it->Pz(),2)+pow(func(it->Pt()),2));
		}
	}
  	return sum;
}


//-------------------------------------------------------------------Resolutions-definitions----------------------------------------------------------------- 


double identity(double E) {                           //If no resolution
	return E;
}

double ResolCPTrack(double Pt) {                            //ResolCP on transverse momentum
	
	return gRandom->Gaus(Pt,5e-5*Pt*Pt);
}

double ResolCPMorg (double E) {                            //ResolCP on energy (Morgunov approximation)

	return gRandom->Gaus(E,1e-4*E*E);
}


class ResolutionStochastiqueEM                        //ResolEM
{
	public:

 	ResolutionStochastiqueEM(double f=1.0) : facteur(f) {}
  	double operator()(double E) {return gRandom->Gaus(E,sqrt((facteur*sqrt(E))*(facteur*sqrt(E))+((0.55/100)*E)*((0.55/100)*E)+(15.5/100)*(15.5/100)));}

	private:

  	double facteur;
};


class ResolutionStochastiqueNH                         //ResolNH
{
	public:

 	ResolutionStochastiqueNH(double f=1.0) : facteur(f) {}
  	double operator()(double E) {return gRandom->Gaus(E,sqrt((facteur*sqrt(E))*(facteur*sqrt(E))+((1.6/100)*E)*((1.6/100)*E)/*+(0.18/100)*(0.18/100)*/));}

	private:

  	double facteur;
};


//---------------------------------------------------------Class-definition-for-energies-histogramms----------------------------------------------------------


struct MesHistogrammes
{
  TH1F EnergyTotaleVraie;
  TH1F EnergyTotaleSmeared;
  TH1F ESmearWithoutNeut;
  TH1F HNEnergySmeasured;
  TH1F CPEnergySmeasured;
  TH1F EMEnergySmeasured;
  MesHistogrammes(float eCDM);
};

MesHistogrammes::MesHistogrammes(float eCDM) : EnergyTotaleVraie("Evraie","Energie Totale Vraie",200,0,2*eCDM),
				     EnergyTotaleSmeared("Esmear","Energie Totale avec resolution",200,0,2*eCDM),
				     ESmearWithoutNeut("EsmearwoutN","Energie totale avec résolution pour les evenements sans neutrinos",200,0,2*eCDM),
				     HNEnergySmeasured("HNmeas","Energy sum distribution of neutral hadrons (measured)",200,0,eCDM),
				     CPEnergySmeasured("CPmeas","Energy sum distribution of charged particles (measured)",200,0,eCDM),
				     EMEnergySmeasured("EMmeas","Energy sum distribution of photons (measured)",200,0,eCDM)		   
{}


//------------------------------------Function-wich-read-the-tree-from-main-then-do-sums-and-fill-histogramms-------------------------------------------------


MesHistogrammes* PlotEnergy(double Ecdm, double ResolVar, int Case, std::ofstream &out)    
{

  //ResolVar is the parameter that will evolve in resolution (NH/EM/CP) formula

  //ROOT file opening
  TFile f1("Simulation.root");
  //get tree
  TTree *t1 = (TTree*)f1.Get("event");

  MesHistogrammes *histos=new MesHistogrammes(Ecdm);   
  
  gRandom->SetSeed(time(0));

  std::vector<int> * pChargedPDGid = 0;
  t1->SetBranchAddress("CP_PDG",&pChargedPDGid);
  std::vector<ROOT::Math::XYZTVector> * pChargedTracks = 0;
  t1->SetBranchAddress("CP_LV",&pChargedTracks);
  std::vector<int> * phnPDGid = 0;
  t1->SetBranchAddress("HN_PDG",&phnPDGid);
  std::vector<ROOT::Math::XYZTVector> * pHNTracks = 0;
  t1->SetBranchAddress("HN_LV",&pHNTracks);
  std::vector<ROOT::Math::XYZTVector> * pGammaTracks = 0;
  t1->SetBranchAddress("GAM_LV",&pGammaTracks);
  std::vector<ROOT::Math::XYZTVector> * pNeutrinoTracks =0;
  t1->SetBranchAddress("NU_LV",&pNeutrinoTracks);

  int n = (int) t1->GetEntries();
  
  for (int i = 0; i < n; ++i)
    {
      t1->GetEntry(i);

      double EtNH=calculSommeEnergie(ResolutionStochastiqueNH(0.847),pHNTracks,1);                  //SUM NH
      if(EtNH>0) {
      	histos->HNEnergySmeasured.Fill(EtNH);
      }

      double EtEM=calculSommeEnergie(ResolutionStochastiqueEM(0.027),pGammaTracks,1);            //SUM EM
      if(EtEM>0) {
      	histos->EMEnergySmeasured.Fill(EtEM);
      }

      double EtCPTrack=calculSommeEnergie(ResolCPTrack,pChargedTracks,2);                                     //SUM CP (Momentum Res)
      if(EtCPTrack>0) {
      	histos->CPEnergySmeasured.Fill(EtCPTrack);
      }
	
      double EtCPMorg=calculSommeEnergie(ResolCPMorg,pChargedTracks,1);                                     //SUM CP (Energy Res)

      double Etot=calculSommeEnergie(identity,pHNTracks,1)+calculSommeEnergie(identity,pGammaTracks,1)+calculSommeEnergie(identity,pChargedTracks,1);

      double Esmear=EtNH+EtEM+EtCPTrack;				//Total Energy Measured (with CP resol on momentum)
	
      histos->EnergyTotaleVraie.Fill(Etot);
      histos->EnergyTotaleSmeared.Fill(Esmear);
	
      if (Etot<Ecdm+1e-3&& Etot>Ecdm-1e-3) {                            //If you want to consider only event with Etot=Ecdm (to have a gaussian)
     	histos-> ESmearWithoutNeut.Fill(Esmear);
	}
    }      

 
     

  if (Case==3) {
  	out<<ResolVar<<" "<<(histos->EnergyTotaleSmeared.GetStdDev())/histos->EnergyTotaleSmeared.GetMean()<<endl;   //Use if you need to draw Restot=f(k);
  }

  return histos;
}


//-----------------------------------------------------------Functions-called-in-main-------------------------------------------------------------------------


void plotEnergy(double Ecdm, int d, int Case, MesHistogrammes *R, float &RESO, float &RESO2) {                        
 
  if (Case==1 || Case==2 || Case==4) {

	std::ofstream f;
	std::ofstream u;
	
  	R=PlotEnergy(Ecdm,0, Case,f);


	//Finding crystall ball parameters to fit to the measured energy histogramm (with RooFit)
	TH1F *H = &R->EnergyTotaleSmeared;
	RooRealVar x("x","x",H->GetBinCenter(H->FindFirstBinAbove(0,1)),H->GetBinCenter(H->FindLastBinAbove(0,1)));
	RooDataHist data("data","data set of x1", x, H);
	RooRealVar cbMean("cbMean", "cbMean", H->GetMean(), H->GetMean()-10, H->GetMean()+10);
	RooRealVar cbSigma("cbSigma", "cbSigma", H->GetStdDev(),0,3*H->GetStdDev());
	RooRealVar alpha("alpha","",1,-5,5);
	RooRealVar n("n","",3,0,5);
	RooCBShape CBall("CBall", "Crystal Ball fit", x, cbMean, cbSigma, alpha, n);
	CBall.fitTo(data);
	//Fitting the histogramm with a crystal ball (with Root)
	TF1 *Fiti=new TF1("Fiti","crystalball");
	Fiti->SetParameters(0,cbMean.getValV(),cbSigma.getValV(),alpha.getValV(),n.getValV());
	H->Fit("Fiti");
	gStyle->SetOptFit();
	
	
  	if (d==0) {                                                //New file Results.root if (Case=1 or 4) or if(Case=2 and d=0)
		
                TFile f1("Results.root","RECREATE");              //Sending new histogramms to a new file Results.root
  		R->EnergyTotaleVraie.Write();
  		H->Write();
		R->ESmearWithoutNeut.Write();
  		R->HNEnergySmeasured.Write();
  		R->EMEnergySmeasured.Write();
  		R->CPEnergySmeasured.Write();
  		f1.Close();
		f.open("output.txt");                       //Usefull to clear the file before writing in
 		f.close();
		u.open("output4.txt");                       //Usefull to clear the file before writing in
 		u.close();
         }

         else {	
	
		TFile f1("Results.root","UPDATE");            //Sending new histogramms by updating the existing file Results.root
  		R->EnergyTotaleVraie.Write();
  		H->Write();
		R->ESmearWithoutNeut.Write();
  		R->HNEnergySmeasured.Write();
  		R->EMEnergySmeasured.Write();
  		R->CPEnergySmeasured.Write();
	 }

   	if (Case==2) {                                     
    		f.open ("output.txt", std::ofstream::out | std::ofstream::app);    //To draw f(EtotMeasured)=sigma/EtotMeasured, use MacroED2.C
   	 	f<<std::endl<<cbMean.getValV()<<" "<<cbSigma.getValV()*cbSigma.getValV()<<" "<<cbMean.getError()<<" "<<3*2*cbSigma.getValV()*cbSigma.getError();
   	 	f.close();
		u.open ("output4.txt", std::ofstream::out | std::ofstream::app);    //To draw f(EtotMeasured)=sigma/EtotMeasured, use MacroED2.C
   	 	u<<std::endl<<R->ESmearWithoutNeut.GetMean()<<" "<<R->ESmearWithoutNeut.GetStdDev()*R->ESmearWithoutNeut.GetStdDev();
   	 	u.close();
	}	

	if (Case==4) {
		RESO=R->ESmearWithoutNeut.GetStdDev();
		RESO2=cbSigma.getValV();
	}
		
  }
  
 
  if(Case==3) {

	string OutF[3]={"output1.txt","output2.txt","output3.txt"};

 	std::ofstream g(OutF[d]);                          //Use if you need to draw Restot=f(k);
	int number=0;
  	for (double k=0; k<=1; k=k+0.02) {
		number++;
        	cout<<"Itération numéro : "<<number<<" ..."<<endl;
 	 	MesHistogrammes* histos=PlotEnergy(Ecdm, k, Case,g);
         	delete histos;
  	}

 	 g.close();  
  }

}
