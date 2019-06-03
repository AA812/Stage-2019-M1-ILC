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

using namespace RooFit ;
using namespace Pythia8;
using namespace ROOT::Math;
using namespace std;


int main()
{

//--------------------------------------------------------------Initialisation-of-the-choice------------------------------------------------------------------

 int Case=affiche();                  //Choose of the case
 int Dr=1;
 if (Case==4) {
	cout<<"How many Res ?"<<endl;
	cin>>Dr;
 }

 std::vector<string>  fileCDM;  //Vector of files used


 if(Case==1) {

        
 	fileCDM={"FILES/file.cmnd"}; //File for one simulation

 }


 if(Case==2 || Case==4) {	

 	fileCDM={"FILES/file30.cmnd","FILES/file50.cmnd","FILES/file80.cmnd","FILES/file100.cmnd","FILES/file150.cmnd","FILES/file300.cmnd","FILES/file600.cmnd","FILES/file1000.cmnd","FILES/file1300.cmnd","FILES/file1750.cmnd","FILES/file2000.cmnd","FILES/file2300.cmnd","FILES/file2750.cmnd","FILES/file3000.cmnd"}; //List of files used to do simulations with different Ecdm
     
 }

 if(Case==3) {
	
 	 fileCDM={"FILES/file200.cmnd","FILES/file350.cmnd","FILES/file500.cmnd"}; //List of files used to draw Restot=f(k)
  
 }

 MesHistogrammes *Results[fileCDM.size()]; 

 float  Reso[fileCDM.size()][Dr];          
 float  Reso2[fileCDM.size()][Dr];      



 for (int d=0; d<fileCDM.size(); d++) {  //Beginning of the simulation loop
	
	for (int k=0; k<Dr; k++) {    //Resolution distribution Loop (optionnal)
  // Generator.
  Pythia pythia;
  
  // Shorthand for the event record in pythia.
  Event& event = pythia.event;


  // Read in commands from external file.                            
    pythia.readFile(fileCDM[d]);    
                                             
  
  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");
  

  // Initialize.
  pythia.init();
  Results[d]=new MesHistogrammes(pythia.info.eCM());

  TFile f1("Simulation.root","RECREATE"); //ROOT file opening

  // create tree
  TTree t1("event", "Arbre des particules");



//--------------------------------------------Initialisation-of-differents-objects-in-function-of-particles-categories----------------------------------------

  //Charged Particle info
  std::vector<int>  chargedPDGid;  
  std::vector<int> * pChargedPDGid = &chargedPDGid;
  t1.Branch("CP_PDG","std::vector<int>",&pChargedPDGid);
  std::vector<ROOT::Math::XYZTVector>  chargedTracks;
  std::vector<ROOT::Math::XYZTVector> * pChargedTracks = &chargedTracks;
  t1.Branch("CP_LV","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pChargedTracks);

  //Neutral hadron info
  std::vector<int>  hnPDGid;
  std::vector<int> * phnPDGid = &hnPDGid;
  t1.Branch("HN_PDG","std::vector<int>",&phnPDGid);
  std::vector<ROOT::Math::XYZTVector>  HNTracks;
  std::vector<ROOT::Math::XYZTVector> * pHNTracks = &HNTracks;
  t1.Branch("HN_LV","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pHNTracks);

  //photon info
  std::vector<ROOT::Math::XYZTVector>  gammaTracks;
  std::vector<ROOT::Math::XYZTVector> * pGammaTracks = &gammaTracks;
  t1.Branch("GAM_LV","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pGammaTracks);

  //neutrino info
  std::vector<int>  NeutPDGid;
  std::vector<int> * pNeutPDGid = &NeutPDGid;
  t1.Branch("NU_PDG","std::vector<int>",&pNeutPDGid);
  std::vector<ROOT::Math::XYZTVector>  neutrinoTracks;
  std::vector<ROOT::Math::XYZTVector> * pNeutrinoTracks = &neutrinoTracks;
  t1.Branch("NU_LV","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pNeutrinoTracks);

  
//---------------------------------------------------------------------Event-loop-----------------------------------------------------------------------------

  // Begin event loop..
  int iAbort = 0;
  cout<<"Nb événement : "<<nEvent<<endl;
  
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
      hnPDGid.clear();
      HNTracks.clear();
      gammaTracks.clear();
      neutrinoTracks.clear();

      //fill vectors
      for (int ipart=0; ipart<event.size(); ++ipart)
	{
	  Particle& pt=event[ipart];

	  if (pt.status()<=0) continue; //to sort only final particles

	  XYZTVector q(pt.px(),pt.py(),pt.pz(),pt.e());

	  if (pt.charge()!=0 && pt.m()!=0)            //Charged Particles
	    {
	      chargedPDGid.push_back(pt.id());
	      chargedTracks.push_back(q);
	    }
	  else if (pt.charge()==0 && pt.m()!=0 && abs(pt.id())>18)  //Neutral Hadrons
	    {
	      hnPDGid.push_back(pt.id());
	      HNTracks.push_back(q);
	    }
	  else if (abs(pt.id())<=18 && pt.charge()==0){ //Neutrinos
		
	      NeutPDGid.push_back(pt.id());					
	      neutrinoTracks.push_back(q);
	    }
	  else if (pt.id()==22)                          //Photons
	    gammaTracks.push_back(q);
	  else
	    {
	      std::cout << " This particle isn't in any case " << pt.nameWithStatus() << " PDG code=" << pt.id() << std::endl;
	    }
	}
      
      t1.Fill();	//filling of the tree at each event
      
    }// End of event loop


//------------------------------------------------------------------------------------------------------------------------------------------------------------
  
		  				
  	
	f1.Write();  //Sending the tree to Simulation.root (with info of each particle of each event)
	f1.Close();

	cout<<"Calculating for events at E= "<<pythia.info.eCM()<<" GeV in CDM..."<<endl;
  	plotEnergy(pythia.info.eCM(),d,Case, Results[d], Reso[d][k], Reso2[d][k]);                        //Allow to study energy results of the simulation
                                             

  // Final statistics.
 // pythia.stat();//

  } //End of case 4 loop
	
 } //End of simulation loop

if (Case==4) {

 TFile f2("Resoluncer.root","RECREATE"); //ROOT file opening

 TH1F *RES[fileCDM.size()];
 TH1F *RES2[fileCDM.size()];

 for (int p=0; p<fileCDM.size();p++) {

	RES[p]= new TH1F("DistribResG", "Distribution de la résolution totale, gauss", 100,0,30);
	RES2[p]= new TH1F("DistribResCB", "Distribution de la résolution totale, CB", 100,0,30);
	
	for (int v=0; v<Dr; v++) {

		RES[p]->Fill(Reso[p][v]);
		RES2[p]->Fill(Reso2[p][v]);
	}
	
	RES[p]->Write();
	RES2[p]->Write();
 }
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------
  for(int i=0; i<fileCDM.size(); i++) {
      delete Results[i]; 
  }
 

  return 0;
}

  
