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
#include "RooGaussian.h"
#include "RooFitResult.h"


using namespace RooFit;
using namespace ROOT::Math;
using namespace std;

//------------------------------------------------------------Function-to-choose-a-case-----------------------------------------------------------------------

int affiche() {

	int Case;

 	cout<<"Choose a case : "<<endl;
 	cout<<"1. Study of one simulation at one energy in center of mass"<<endl;
 	cout<<"2. Study of several simulations at different energy in center of mass"<<endl;
 	cout<<"3. Study of the effects of the variation of a resolution parameter"<<endl;
	do {
 	cin>>Case;
	} while (Case<1 || Case>3);
	
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

 	ResolutionStochastiqueEM(double Stochas, double Constant,double Noise);
  	double operator()(double E) {return gRandom->Gaus(E,sqrt((TStochas*sqrt(E))*(TStochas*sqrt(E))+(TConstant*E)*(TConstant*E)+(TNoise)*(TNoise)));}

	private:

  	double TStochas;
	double TNoise;
	double TConstant;
};

ResolutionStochastiqueEM::ResolutionStochastiqueEM(double Stochas, double Constant, double Noise) {
		TStochas=Stochas;
		TNoise=Noise;
		TConstant=Constant;
}

class ResolutionStochastiqueNH                         //ResolNH
{
	public:

 	ResolutionStochastiqueNH(double Stochas,double Constant,double Noise);
  	double operator()(double E) {return gRandom->Gaus(E,sqrt((TStochas*sqrt(E))*(TStochas*sqrt(E))+(TConstant*E)*(TConstant*E)+(TNoise)*(TNoise)));}

	private:

  	double TStochas;
	double TNoise;
	double TConstant;
};

ResolutionStochastiqueNH::ResolutionStochastiqueNH(double Stochas, double Constant,double Noise) {
		TStochas=Stochas;
		TNoise=Noise;
		TConstant=Constant;
}


//---------------------------------------------------------Class-definition-for-energies-histogramms----------------------------------------------------------


struct MesHistogrammes
{
  TH1F EnergyTotaleVraie;
  TH1F EnergyTotaleSmeared;
  TH1F ESmearWithoutNeut;
  TH1F HNEnergySmeasured;
  TH1F CPEnergySmeasured;	
  TH1F EMEnergySmeasured;
  TH1F ResDis;
  TH1F ResDisNoNeut;
  MesHistogrammes(float eCDM);
};

MesHistogrammes::MesHistogrammes(float eCDM) : EnergyTotaleVraie("Evraie","Energie Totale Vraie",200,0,2*eCDM),
				     EnergyTotaleSmeared("Esmear","Energie Totale avec resolution",200,0,2*eCDM),
				     ESmearWithoutNeut("EsmearwoutN","Energie totale avec résolution pour les evenements sans neutrinos",200,0,2*eCDM),
				     HNEnergySmeasured("HNmeas","Energy sum distribution of neutral hadrons (measured)",200,0,eCDM),
				     CPEnergySmeasured("CPmeas","Energy sum distribution of charged particles (measured)",200,0,eCDM),
				     EMEnergySmeasured("EMmeas","Energy sum distribution of photons (measured)",200,0,eCDM),
				     ResDis("ResDis", "Distribution of totale resolution (sampled)",200,0,30),
				     ResDisNoNeut("ResDisNoNeut","Distribution of totale resolution (no neutrinos) (sampled)", 200,0,30)		   
{}

//-------------------------------------------------------------------Class-for-fitting------------------------------------------------------------------------

class PlotFitGauss {

 public :

 PlotFitGauss (TH1F *H);
 float Mean;
 float MeanErr;
 float Sigma;
 float SigmaErr;
 int FitStatus;
};

PlotFitGauss::PlotFitGauss(TH1F *H) {

	//Finding gaussian parameters to fit to the measured energy histogramm (with RooFit)
	RooRealVar x("x","x",H->GetBinCenter(H->FindFirstBinAbove(0,1)),H->GetBinCenter(H->FindLastBinAbove(0,1)));
	RooDataHist data("data","data set of x1", x, H);
	RooRealVar gMean("gMean", "gMean", H->GetMean(), H->GetMean()-10, H->GetMean()+10);
	RooRealVar gSigma("gSigma", "gSigma", H->GetStdDev(),0,3*H->GetStdDev());
	RooGaussian Gauss("Gauss", "Gaussian fit", x, gMean, gSigma);
	RooFitResult* K=Gauss.fitTo(data,PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1),Save());
	//Fitting the histogramm with a gaussian (with Root)
	TF1* Fiti=new TF1("Fiti","gaus");
	Fiti->SetParameter(1,gMean.getValV());
	Fiti->SetParameter(2,gSigma.getValV());
	H->Fit("Fiti","QM0");
	gStyle->SetOptFit();					  //Show fit parameters on cd
	if (Fiti->GetParameter(2)>H->GetStdDev()) {               //In this condition, Gaussian is wrong, so we keep histos values
			Mean=H->GetMean();
			Sigma=H->GetStdDev();
			MeanErr=H->GetMeanError();
			SigmaErr=H->GetStdDevError();
		} else {
			H->Fit("Fiti","QM");
			Mean=Fiti->GetParameter(1);
 			Sigma=Fiti->GetParameter(2);
			MeanErr=Fiti->GetParError(1);
			SigmaErr=Fiti->GetParError(2);
		}
		FitStatus=K->status();
	delete K;
	delete Fiti;
 }

class PlotFitCB {

 public :

 PlotFitCB (TH1F *H);
 float Mean;
 float MeanErr;
 float Sigma;
 float SigmaErr;
 float alpha;
 int FitStatus;

};

PlotFitCB::PlotFitCB(TH1F *H) {

	//Finding crystall ball parameters to fit to the measured energy histogramm (with RooFit)
	RooRealVar x("x","x",H->GetBinCenter(H->FindFirstBinAbove(0,1)),H->GetBinCenter(H->FindLastBinAbove(0,1)));
	RooDataHist data("data","data set of x1", x, H);
	RooRealVar cbMean("cbMean", "cbMean", H->GetMean(), H->GetMean()-10, H->GetMean()+10);
	RooRealVar cbSigma("cbSigma", "cbSigma", H->GetStdDev(),0,3*H->GetStdDev());
	RooRealVar cbalpha("cbalpha","",1,0,5);
	RooRealVar n("n","",3,0,5);
	RooCBShape CBall("CBall", "Crystal Ball fit", x, cbMean, cbSigma, cbalpha, n);
	RooFitResult* K=CBall.fitTo(data,PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1),Save());
	//Fitting the histogramm with a crystal ball (with Root)
	TF1* Fiti=new TF1("Fiti","crystalball");
	Fiti->SetParameter(1,cbMean.getValV());
	Fiti->SetParameter(2,cbSigma.getValV());
	Fiti->SetParameter(3,cbalpha.getValV());
	Fiti->SetParameter(4,n.getValV());
	H->Fit("Fiti","QM");
	gStyle->SetOptFit();					  //Show fit parameters on cd
	Mean=Fiti->GetParameter(1);
 	Sigma=Fiti->GetParameter(2);
	MeanErr=Fiti->GetParError(1);
	SigmaErr=Fiti->GetParError(2);
	alpha=Fiti->GetParameter(3);
	FitStatus=K->status();
	delete K;
	delete Fiti;
 }


//--------------------------------------------------------------Magnetic-field-effect-------------------------------------------------------------------------


float RCourbe(float Pt, float B, float q) {                           //Radius of curvature
	float c=299792458;        //Light celerity
	float J=1.60218e-10;	  //Conversion GeV to Joules
	float PtR=Pt*J/c;	  //[PtR]=m*kg/s
	float e=1.602176634e-19; //elementary charge
	return (PtR/(B*q*e));
}

vector<double> Center(float R, ROOT::Math::XYZTVector Q, vector<double> PosProd) {            //Calculating the center of curvature with momentum

	float c=299792458;        //Light celerity
 	float J=1.60218e-10;	  //Conversion GeV to Joules
        vector<double> P={Q.X()*J/c, Q.Y()*J/c, Q.Z()*J/c};
	double xC=PosProd[0]-(R*(-P[1]))/sqrt(P[0]*P[0]+P[1]*P[1]);
	double yC=PosProd[1]-(R*P[0])/sqrt(P[0]*P[0]+P[1]*P[1]);
	vector<double> C;
	C.push_back(xC);
	C.push_back(yC);
	C.push_back(PosProd[2]);
	return C;
}
	
class Cercle {
 public :
 friend XYZTVector intersect  (Cercle C1, Cercle C2,vector<double> Pos, XYZTVector Q);
 Cercle(float Ray, vector<double> Cent);
 private :
 float Rayon;
 std::vector<double> Centre;

};

 Cercle::Cercle(float Ray, vector<double> Cent) {
	Rayon=Ray;
	Centre=Cent;
 }



 double Angle(vector<double> U, vector<double> V) {                             //Return the angle between U and V
	double normeU=sqrt(U[0]*U[0]+U[1]*U[1]+U[2]*U[2]);
	double normeV=sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
	double ProdScal=U[0]*V[0]+U[1]*V[1]+U[2]*V[2];
	double theta=ProdScal/(normeU*normeV);
	return std::acos(theta);
 }


XYZTVector intersect (Cercle C1, Cercle C2, vector <double> Pos, XYZTVector Q) { //Check if Circle C1 and C2 are intersecting and return the position of detection (CHARGED PARTICLES)

 float c=299792458;        //Light celerity
 float J=1.60218e-10;	  //Conversion GeV to Joules
 float N=(pow(C2.Rayon,2)-pow(C1.Rayon,2)-pow(C2.Centre[0],2)+pow(C1.Centre[0],2)-pow(C2.Centre[1],2)+pow(C1.Centre[1],2))/(2*(C1.Centre[1]-C2.Centre[1]));
 float A=pow((C1.Centre[0]-C2.Centre[0])/(C1.Centre[1]-C2.Centre[1]),2)+1;	
 float B=2*C1.Centre[1]*((C1.Centre[0]-C2.Centre[0])/(C1.Centre[1]-C2.Centre[1]))-2*N*((C1.Centre[0]-C2.Centre[0])/(C1.Centre[1]-C2.Centre[1]))-2*C1.Centre[0];
 float C=(pow(C1.Centre[0],2)+pow(C1.Centre[1],2)+pow(N,2)-pow(C1.Rayon,2)-2*C1.Centre[1]*N);
 float Delta=pow(B,2)-4*A*C;
 XYZTVector x1;
 XYZTVector x2;
 XYZTVector null(0,0,0,0);
 vector<double> P={Q.X()*J/c,Q.Y()*J/c,Q.Z()*J/c};
 vector<double> Pz={0,0,Q.Z()*J/c};

 if (Delta>=0) {

	x1.SetPx((-B-sqrt(Delta))/(2*A));
	x1.SetPy((N-x1.X()*((C1.Centre[0]-C2.Centre[0])/(C1.Centre[1]-C2.Centre[1]))));

	x2.SetPx((-B+sqrt(Delta))/(2*A));
	x2.SetPy((N-x2.X()*((C1.Centre[0]-C2.Centre[0])/(C1.Centre[1]-C2.Centre[1]))));

	double o1=sqrt(x1.X()*x1.X()+x1.Y()*x1.Y());
	double o2=sqrt(x2.X()*x2.X()+x2.Y()*x2.Y());
	
	x1.SetPz(o1/tan(Angle(P,Pz)));
	x2.SetPz(o1/tan(Angle(P,Pz)));
	x1.SetE(Q.E());
	x2.SetE(Q.E());
 
	if (sqrt(pow(x1.X()-Pos[0],2)+pow(x1.Y()-Pos[1],2))<=sqrt(pow(x2.X()-Pos[0],2)+pow(x2.Y()-Pos[1],2))) {
		return x1;
	} else { return x2;}
}

else {return null;}

}


XYZTVector NHDetection (XYZTVector Q,double R) {         //Position of the detection of a neutral Hadron

	float c=299792458;        //Light celerity
 	float J=1.60218e-10;	  //Conversion GeV to Joules
	vector<double> P={Q.X()*J/c,Q.Y()*J/c,Q.Z()*J/c};
	vector<double> Pz={0,0,Q.Z()*J/c};
	XYZTVector PosDetec;
	PosDetec.SetPx((R*P[0])/sqrt(P[0]*P[0]+P[1]*P[1]));
	PosDetec.SetPy((R*P[1])/sqrt(P[0]*P[0]+P[1]*P[1]));
	
	double o=sqrt(PosDetec.X()*PosDetec.X()+PosDetec.Y()*PosDetec.Y());
	PosDetec.SetPz(o/tan(Angle(P,Pz)));
	PosDetec.SetE(Q.E());
	return PosDetec;
}


std::vector<ROOT::Math::XYZTVector> *ConfusionNHCP(std::vector<ROOT::Math::XYZTVector> &vecNH,std::vector<ROOT::Math::XYZTVector> &vecCP) {	
	std::vector<ROOT::Math::XYZTVector> NewNH;
	for (int itNH=0; itNH<vecNH.size(); itNH++){
		int NoDetec=0;
		for (int itCP=0; itCP<vecCP.size(); itCP++){
			double distDetectMax=0.05;
			double DIST=sqrt((vecCP[itCP].X()-vecNH[itNH].X())*(vecCP[itCP].X()-vecNH[itNH].X())+(vecCP[itCP].Y()-vecNH[itNH].Y())*(vecCP[itCP].Y()-vecNH[itNH].Y())+(vecCP[itCP].Z()-vecNH[itNH].Z())*(vecCP[itCP].Z()-vecNH[itNH].Z()));
			if (DIST<distDetectMax || vecNH[itNH].E()<(sqrt( vecCP[itCP].E())/3)) {
				NoDetec++;
			}
		}
		if (NoDetec==0) {
			XYZTVector NHDetec(vecNH[itNH].X(),vecNH[itNH].Y(),vecNH[itNH].Z(),vecNH[itNH].E());
			NewNH.push_back(NHDetec);
		}	
	}
	std::vector<ROOT::Math::XYZTVector> *TEST=&NewNH;
	return TEST;
}

//------------------------------------Function-wich-read-the-tree-from-main-then-do-sums-and-fill-histogramms-------------------------------------------------


MesHistogrammes* PlotEnergy(double Ecdm, double ResolVar, int Case, std::ofstream &out, const char* FileName)    
{

  //ResolVar is the parameter that will evolve in resolution (NH/EM/CP) formula


  TFile *f1=new TFile(FileName);
  //get tree
  TTree *t1=NULL;
  t1 = (TTree*)f1->Get("event");
  
  MesHistogrammes *histos=new MesHistogrammes(Ecdm);
  
  gRandom->SetSeed(time(0));

  std::vector<int> * pChargedPDGid = 0;
  t1->SetBranchAddress("CP_PDG",&pChargedPDGid);
  std::vector<ROOT::Math::XYZTVector> * pChargedTracks = 0;
  t1->SetBranchAddress("CP_LV",&pChargedTracks);
  std::vector<ROOT::Math::XYZTVector> CPDetectPos;
  std::vector<ROOT::Math::XYZTVector> *pCPDetectPos=0;
  t1->SetBranchAddress("CP_DetectPos",&pCPDetectPos);
  std::vector<int> * phnPDGid = 0;
  t1->SetBranchAddress("HN_PDG",&phnPDGid);
  std::vector<ROOT::Math::XYZTVector> *pHNTracks = 0;
  t1->SetBranchAddress("HN_LV",&pHNTracks);
  std::vector<ROOT::Math::XYZTVector> HNDetectPos;
  std::vector<ROOT::Math::XYZTVector> *pHNDetectPos=0;
  t1->SetBranchAddress("HN_DetectPos",&pHNDetectPos);
  std::vector<ROOT::Math::XYZTVector> * pGammaTracks = 0;
  t1->SetBranchAddress("GAM_LV",&pGammaTracks);
  std::vector<ROOT::Math::XYZTVector> * pNeutrinoTracks =0;
  t1->SetBranchAddress("NU_LV",&pNeutrinoTracks);

  int n = (int) t1->GetEntries();
  vector<float> EsmearTot;               //Vecto that will contains all values of Esmear to sample
  vector<float> EnoNeutTot;


  for (int i = 0; i < n; ++i)
    {
      t1->GetEntry(i);

      //if no confusion
     // double EtNH=calculSommeEnergie(ResolutionStochastiqueNH(0.443,1.8/100,0.18/100),pHNTracks,1);                  //SUM NH
      //with confusion on nh because of cp (same event)
      std::vector<ROOT::Math::XYZTVector> *ConfuNH=ConfusionNHCP(*pHNDetectPos,*pCPDetectPos);
      double EtNH=calculSommeEnergie(ResolutionStochastiqueNH(0.443,1.8/100,0.18/100),ConfuNH,1);


      if(EtNH>0) {
      	histos->HNEnergySmeasured.Fill(EtNH);
      }

      double EtEM=calculSommeEnergie(ResolutionStochastiqueEM(0.166,1.1/100,0),pGammaTracks,1);            //SUM EM
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
      EsmearTot.push_back(Esmear);
	
      if (Etot<Ecdm+1e-3&& Etot>Ecdm-1e-3) {                            //If you want to consider only event with Etot=Ecdm (to have a gaussian)
     	histos-> ESmearWithoutNeut.Fill(Esmear);
	EnoNeutTot.push_back(Esmear);
	}
    }  
    
    if (Case==1 || Case==2) {
      cout<<"Determination of uncertainties on total resolution..."<<endl;
      

      for (int k=0; k<10000; k++) {         //Determination of uncertainties of the total resolution ( number of entries in resolution distribution)

	if (k%100==0) {
		cout<<"Sample number "<<k<<endl;
	}

	TH1F *SampleEsmear=new TH1F("SampleEsmear","Sample of Esmear measure (random)", 200, 0, 2*Ecdm);
	TH1F *SampleNoNeutE=new TH1F("SampleNoNeutE", "Sample of ENoNeut (random)", 200,0,2*Ecdm);
	int s, p;
	vector <float> Eg=EsmearTot;
	vector <float> En=EnoNeutTot;
	for (int c=0; c<30000; c++) {				//Number of entries in sampled histos
		do {
			s=gRandom->Integer(Eg.size());
		} while (Eg[s]==0);                              //Never fill with the same value
		do { 
			p=gRandom->Integer(En.size());
		} while (En[p]==0);
		SampleEsmear->Fill(Eg[s]);
		Eg[s]=0;
		SampleNoNeutE->Fill(En[p]);
		En[p]=0;
	}
	
	PlotFitCB *ResPara= new PlotFitCB(SampleEsmear);
	if (ResPara->FitStatus==0 && ResPara->alpha>0) {                          //Check if the plot worked
		histos->ResDis.Fill(ResPara->Sigma);
	}
	PlotFitGauss *ResPara2 = new PlotFitGauss(SampleNoNeutE);
	if (ResPara2->FitStatus==0) {                         //Check if the plot worked
		histos->ResDisNoNeut.Fill(ResPara2->Sigma);
	}

	delete ResPara;
        delete ResPara2;
	delete SampleEsmear;
	delete SampleNoNeutE;
    }
  }	
     

  if (Case==3) {
	TH1F *pointer=&histos->EnergyTotaleSmeared;
	PlotFitCB CB(pointer);
  	out<<ResolVar<<" "<<CB.Sigma/CB.Mean<<endl;   //Use if you need to draw Restot=f(k);
  }

  return histos;
 
  delete f1;
}


	
//-----------------------------------------------------------Functions-called-in-main-------------------------------------------------------------------------


void plotEnergy(double Ecdm, int d, int Case, MesHistogrammes *R, const char* FileName) {                        
 
  if (Case==1 || Case==2) {

	std::ofstream f;
	std::ofstream u;
	
  	R=PlotEnergy(Ecdm,0, Case,f,FileName);                             //Cacul of histogramms
	TH1F *ETotSMEAR=&R->EnergyTotaleSmeared;
	PlotFitCB *CB=NULL;
	int Bug=0;
	do {
	     CB=new PlotFitCB(ETotSMEAR);
	     Bug++;
	     if (Bug==100000) {
		break;
	     }
	} while (CB->alpha<=0 || CB->alpha>10);
	
	TH1F* DistribRes=&R->ResDis;
	PlotFitGauss *G1=new PlotFitGauss(DistribRes);		

  	if (d==0) {                                                //New file Results.root if (Case=1 or 4) or if(Case=2 and d=0)
		
                TFile f1("Results.root","RECREATE");              //Sending new histogramms to a new file Results.root
  		R->EnergyTotaleVraie.Write();
  		R->EnergyTotaleSmeared.Write();
		R->ESmearWithoutNeut.Write();
  		R->HNEnergySmeasured.Write();
  		R->EMEnergySmeasured.Write();
  		R->CPEnergySmeasured.Write();
		R->ResDis.Write();
		R->ResDisNoNeut.Write();
  		f1.Close();
		f.open("output12.txt");                       //Usefull to clear the file before writing in
 		f.close();
		u.open("output13.txt");                       //Usefull to clear the file before writing in
 		u.close();
         }

         else {	
	
		TFile f1("Results.root","UPDATE");            //Sending new histogramms by updating the existing file Results.root
  		R->EnergyTotaleVraie.Write();
  		R->EnergyTotaleSmeared.Write();
		R->ESmearWithoutNeut.Write();
  		R->HNEnergySmeasured.Write();
  		R->EMEnergySmeasured.Write();
  		R->CPEnergySmeasured.Write();
		R->ResDis.Write();
		R->ResDisNoNeut.Write();
	 }

   	if (Case==2) {                                     
    		f.open ("output12.txt", std::ofstream::out | std::ofstream::app);    //To draw f(EtotMeasured)=sigma/EtotMeasured, use MacroED2.C
   	 	f<<std::endl<<CB->Mean<<" "<<G1->Mean*G1->Mean<<" "<<2*G1->Mean*G1->Sigma;
   	 	f.close();
		u.open ("output13.txt", std::ofstream::out | std::ofstream::app);    //To draw f(EtotMeasured)=sigma/EtotMeasured, use MacroED2.C
   	 	u<<std::endl<<R->ESmearWithoutNeut.GetMean()<<" "<<R->ResDisNoNeut.GetMean()*R->ResDisNoNeut.GetMean()<<" "<<2*R->ResDisNoNeut.GetMean()*R->ResDisNoNeut.GetStdDev();
   	 	u.close();
	}	
	

	delete CB;
	delete G1;	
  }
  
 
  if(Case==3) {

	string OutF[3]={"output1.txt","output2.txt","output3.txt"};    //One .txt by energy

 	std::ofstream g(OutF[d]);                          //Use if you need to draw Restot=f(k);
	int number=0;
  	for (double k=0; k<=1; k=k+0.04) {
		number++;
        	cout<<"Iteration number : "<<number<<" ..."<<endl;
 	 	MesHistogrammes* histos=PlotEnergy(Ecdm, k, Case,g,FileName);
  	}

 	 g.close();  
  }

}
