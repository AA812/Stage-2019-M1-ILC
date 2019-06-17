////////
///
///  Tire en partie de  tutorials/math/mathcoreVectorCollection.C 
////////

#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "Math/Vector4D.h"
#include "TH1F.h"
#include "PlotEnergy.cc"
#include "TNtuple.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TH1D.h"
#include "RooGlobalFunc.h"
#include "TStyle.h"

using namespace RooFit ;
using namespace Pythia8;
using namespace ROOT::Math;
using namespace std;


int main()
{

 RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ; //Avoid RooFit display on terminal

//--------------------------------------------------------------Initialisation-of-the-choice------------------------------------------------------------------

 int Case=affiche();                  //Choose of the case
 int Case2;
 
 cout<<endl<<endl<<"1. Do simulations"<<endl;
 cout<<"2. Use existing files"<<endl;
 do {
 	cin>>Case2;
 } while (Case2<1 || Case2>2);

 std::vector<string> fileCDM;  //Vector of files used
 std::vector<const char*> SimuFiles; //Vector where Simulation files are stored

 if(Case==1) {

 	fileCDM={"FILES/file.cmnd"}; //File for one simulation
	SimuFiles={"Simulations/Simu.root"};
 }


 if(Case==2) {	

 	fileCDM={"FILES/file30.cmnd","FILES/file50.cmnd","FILES/file80.cmnd","FILES/file100.cmnd","FILES/file150.cmnd","FILES/file300.cmnd","FILES/file600.cmnd","FILES/file1000.cmnd","FILES/file1300.cmnd","FILES/file1750.cmnd","FILES/file2000.cmnd","FILES/file2300.cmnd","FILES/file2750.cmnd","FILES/file3000.cmnd"}; //List of files used to do simulations with different Ecdm
     	SimuFiles={"Simulations/Simu30.root","Simulations/Simu50.root","Simulations/Simu80.root","Simulations/Simu100.root","Simulations/Simu150.root","Simulations/Simu300.root","Simulations/Simu600.root","Simulations/Simu1000.root","Simulations/Simu1300.root","Simulations/Simu1750.root","Simulations/Simu2000.root","Simulations/Simu2300.root","Simulations/Simu2750.root","Simulations/Simu3000.root"};
 }

 if(Case==3) {
	
 	 fileCDM={"FILES/file200.cmnd","FILES/file350.cmnd","FILES/file500.cmnd"}; //List of files used to draw Restot=f(k)
  	 SimuFiles={"Simulations/Simu200.root","Simulations/Simu350.root","Simulations/Simu500.root"};
 }

 MesHistogrammes *Results[fileCDM.size()];      //Histogramms for each Ecdm

 for (int d=0; d<fileCDM.size(); d++) {  //---------------------Beginning of the simulation loop--------------------------------------------------------------
	
  // Generator.
  Pythia pythia;
  
  // Shorthand for the event record in pythia.
  Event& event = pythia.event;


  // Read in commands from external file.                            
    pythia.readFile(fileCDM[d]);                                           
  
  // Extract settings to be used in the main program.
  int nEvent= pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

  Results[d]=new MesHistogrammes(pythia.info.eCM());

if (Case2==1) {                                    //Do simulations or not ?

  TFile *f1=new TFile(SimuFiles[d],"RECREATE"); //ROOT file opening

  // create tree
  TTree *t1=new TTree("event", "Arbre des particules");


//--------------------------------------------Initialisation-of-differents-objects-in-function-of-particles-categories----------------------------------------

  //Detector geometry infos
  vector<double> origin={0,0,0};
  double Rdec=1.8;
  Cercle Tracker(Rdec,origin);
  float MF=4; //Magnetic Field

  //Charged Particle infos
  std::vector<int>  chargedPDGid;  
  std::vector<int> * pChargedPDGid = &chargedPDGid;
  t1->Branch("CP_PDG","std::vector<int>",&pChargedPDGid);
  std::vector<ROOT::Math::XYZTVector>  chargedTracks;
  std::vector<ROOT::Math::XYZTVector> * pChargedTracks = &chargedTracks;
  t1->Branch("CP_LV","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pChargedTracks);
  std::vector<ROOT::Math::XYZTVector>  CPDetectPos;
  std::vector<ROOT::Math::XYZTVector> * pCPDetectPos = &CPDetectPos;
  t1->Branch("CP_DetectPos","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pCPDetectPos);

  //Neutral hadron infos
  std::vector<int>  hnPDGid;
  std::vector<int> * phnPDGid = &hnPDGid;
  t1->Branch("HN_PDG","std::vector<int>",&phnPDGid);
  std::vector<ROOT::Math::XYZTVector>  HNTracks;
  std::vector<ROOT::Math::XYZTVector> * pHNTracks = &HNTracks;
  t1->Branch("HN_LV","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pHNTracks);
  std::vector<ROOT::Math::XYZTVector>  HNDetectPos;
  std::vector<ROOT::Math::XYZTVector> * pHNDetectPos = &HNDetectPos;
  t1->Branch("HN_DetectPos","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pHNDetectPos);


  //photon infos
  std::vector<ROOT::Math::XYZTVector>  gammaTracks;
  std::vector<ROOT::Math::XYZTVector> * pGammaTracks = &gammaTracks;
  t1->Branch("GAM_LV","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pGammaTracks);

  //neutrino infos
  std::vector<int>  NeutPDGid;
  std::vector<int> * pNeutPDGid = &NeutPDGid;
  t1->Branch("NU_PDG","std::vector<int>",&pNeutPDGid);
  std::vector<ROOT::Math::XYZTVector>  neutrinoTracks;
  std::vector<ROOT::Math::XYZTVector> * pNeutrinoTracks = &neutrinoTracks;
  t1->Branch("NU_LV","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pNeutrinoTracks);

  
//---------------------------------------------------------------------Event-loop-----------------------------------------------------------------------------


  // Begin event loop..
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
    {
      // Generate events. Quit if many failures.
      if (!pythia.next())
	{
	  if (++iAbort < nAbort) continue;
	  cout << " Event generation aborted prematurely, owing to error!\n";
	  break;
	}

      //clear vector
      chargedPDGid.clear();
      chargedTracks.clear();
      CPDetectPos.clear();
      hnPDGid.clear();
      HNTracks.clear();
      HNDetectPos.clear();
      gammaTracks.clear();
      neutrinoTracks.clear();

      //fill vectors
      for (int ipart=0; ipart<event.size(); ++ipart)
	{
	  Particle& pt=event[ipart];
	
	  if (pt.status()<=0) continue; //to sort only final particles

	  XYZTVector q(pt.px(),pt.py(),pt.pz(),pt.e()); //Care, its in GeV
          vector<double> PosProd={pt.xProd()*(1e-3),pt.yProd()*(1e-3),pt.zProd()*(1e-3)};  //mm to meter

	  if (pt.charge()!=0 && pt.m()!=0)            //Charged Particles
	    {
	      XYZTVector null(0,0,0,0);
	      chargedPDGid.push_back(pt.id());
	      chargedTracks.push_back(q);
	      float Rcourbure=RCourbe(pt.pT(),MF,pt.charge());          //Charged Particles interaction with magnetic field
	      Cercle PartCurve(Rcourbure,Center(Rcourbure,q,PosProd));	
	      XYZTVector DetectPos=intersect(Tracker,PartCurve,PosProd,q);
              if (DetectPos!=null) {
	      CPDetectPos.push_back(DetectPos);  //Where the particle is detected OR NOT
		}
	    }

	  else if (pt.charge()==0 && pt.m()!=0 && abs(pt.id())>18)  //Neutral Hadrons
	    {
	      hnPDGid.push_back(pt.id());
	      HNTracks.push_back(q);
	      XYZTVector DetectPos=NHDetection(q,Rdec);
	      HNDetectPos.push_back(DetectPos);
	    }
	  else if (abs(pt.id())<=18 && pt.charge()==0){ //Neutrinos
		
	      NeutPDGid.push_back(pt.id());					
	      neutrinoTracks.push_back(q);
	    }
	  else if (pt.id()==22)  {                        //Photons
	    gammaTracks.push_back(q);
		}
	  else
	    {
	      std::cout << " This particle isn't in any case " << pt.nameWithStatus() << " PDG code=" << pt.id() << std::endl;
	    }
	}
      
      t1->Fill();	//filling of the tree at each event
      
    }// End of event loop

  	f1->Write();  //Sending the tree to Simulation.root (with info of each particle of each event)
  	f1->Close();

	delete f1;
	f1=NULL;
	

   } //	The tree is full
//--------------------------------------------------------------Studying-energy-results-of-simulations--------------------------------------------------------

	cout<<"Calculating for simulation at E= "<<pythia.info.eCM()<<" GeV in CM..."<<endl;
  	plotEnergy(pythia.info.eCM(),d,Case, Results[d],SimuFiles[d]);                    //Allow to study energy results of the simulation
	
  } //End of simulation loop

for (int d=0; d<fileCDM.size(); d++) {
	delete Results[d];
}

  // Final statistics.
 // pythia.stat();//

  return 0;
}

  
