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


using namespace ROOT::Math;
using namespace std;


//-------------------------------------------------------Function to calculate the sum with the right resol---------------------------------------------------


template <class F> double calculSommeEnergie(F func,std::vector<ROOT::Math::XYZTVector> *vec)         
{
	double sum=0;
	for (std::vector<ROOT::Math::XYZTVector>::iterator it=vec->begin(); it != vec->end(); ++it){
	    	sum+=func(it->energy());
	}
  	return sum;
}


//-------------------------------------------------------------------Resolutions-definitions----------------------------------------------------------------- 


double identity(double E) {                           //If no resolution
	return E;
}

double ResolCP(double E) {                            //ResolCP
	return gRandom->Gaus(E,1e-4*E*E);
}

class ResolutionStochastiqueEM                        //ResolEM
{
	public:

 	ResolutionStochastiqueEM(double f=1.0) : facteur(f) {}
  	double operator()(double E) {return gRandom->Gaus(E,sqrt((facteur*sqrt(E))*(facteur*sqrt(E))+((1.1/100)*E)*((1.1/100)*E)/*+(15.5/100)*(15.5/100)*/));}

	private:

  	double facteur;
};


class ResolutionStochastiqueNH                         //ResolNH
{
	public:

 	ResolutionStochastiqueNH(double f=1.0) : facteur(f) {}
  	double operator()(double E) {return gRandom->Gaus(E,sqrt((facteur*sqrt(E))*(facteur*sqrt(E))+((1.8/100)*E)*((1.8/100)*E)+(0.18/100)*(0.18/100)));}

	private:

  	double facteur;
};


//---------------------------------------------------------Class-definition-for-energies-histogramms----------------------------------------------------------


struct MesHistogrammes
{
  TH1F EnergyTotaleVraie;
  TH1F EnergyTotaleSmeared;
  TH1F HNEnergySmeasured;
  TH1F CPEnergySmeasured;
  TH1F EMEnergySmeasured;
  MesHistogrammes(float eCDM);
};

MesHistogrammes::MesHistogrammes(float eCDM) : EnergyTotaleVraie("Evraie","Energie Totale Vraie",200,0,2*eCDM),
				     EnergyTotaleSmeared("Esmear","Energie Totale avec resolution",200,0,2*eCDM),
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

      double EtNH=calculSommeEnergie(ResolutionStochastiqueNH(0.443),pHNTracks);                  //SUM NH
      if(EtNH>0) {
      	histos->HNEnergySmeasured.Fill(EtNH);
      }

      double EtEM=calculSommeEnergie(ResolutionStochastiqueEM(ResolVar),pGammaTracks);            //SUM EM
      if(EtEM>0) {
      	histos->EMEnergySmeasured.Fill(EtEM);
      }

      double EtCP=calculSommeEnergie(ResolCP,pChargedTracks);                                     //SUM CP
      if(EtCP>0) {
      	histos->CPEnergySmeasured.Fill(EtCP);
      }

      double Etot=calculSommeEnergie(identity,pHNTracks)+calculSommeEnergie(identity,pGammaTracks)+calculSommeEnergie(identity,pChargedTracks);

      double Esmear=EtNH+EtEM+EtCP;

	
      if (Etot<Ecdm+1e-3&& Etot>Ecdm-1e-3) {                            //We consider only event with Etot=Ecdm (to have a gaussian)
      	histos->EnergyTotaleVraie.Fill(Etot);
      	histos->EnergyTotaleSmeared.Fill(Esmear);
	}
    }      

 
     

  if (Case==3) {
  	out<<ResolVar<<" "<<(histos->EnergyTotaleSmeared.GetStdDev())/histos->EnergyTotaleSmeared.GetMean()<<endl;   //Use if you need to draw Restot=f(k);
  }

  return histos;
}


//-----------------------------------------------------------Functions-called-in-main-------------------------------------------------------------------------



void plotEnergy(double Ecdm, int d, int Case, MesHistogrammes *R) {                        
 
  if (Case==1 || Case==2) {

	std::ofstream f;
  	R=PlotEnergy(Ecdm,0, Case,f);

        //TF1 *Fiti= new TF1("Fiti","gaus");                  //Use if you need to fit
  	//R->EnergyTotaleSmeared.Fit("Fiti");
        
  	if (d==0) {                                                //New file Results.root if (Case=1) or if(Case=2 and d=0)
		
                TFile f1("Results.root","RECREATE");              //Sending new histogramms to a new file Results.root
  		R->EnergyTotaleVraie.Write();
  		R->EnergyTotaleSmeared.Write();
  		R->HNEnergySmeasured.Write();
  		R->EMEnergySmeasured.Write();
  		R->CPEnergySmeasured.Write();
  		f1.Close();
		f.open("output.txt");                       //Usefull to clear the file before writing in
 		f.close();
         }

         else {	
	
		TFile f1("Results.root","UPDATE");            //Sending new histogramms by updating the existing file Results.root
  		R->EnergyTotaleVraie.Write();
  		R->EnergyTotaleSmeared.Write();
  		R->HNEnergySmeasured.Write();
  		R->EMEnergySmeasured.Write();
  		R->CPEnergySmeasured.Write();
  		f1.Close();
	 }

   	if (Case==2) {
    		f.open ("output.txt", std::ofstream::out | std::ofstream::app);    //To draw f(EtotMeasured)=sigma/EtotMeasured, use MacroED2.C
   	 	f<<std::endl<<R->EnergyTotaleSmeared.GetMean()<<" "<<(R->EnergyTotaleSmeared.GetStdDev())/R->EnergyTotaleSmeared.GetMean();
   	 	f.close();
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
