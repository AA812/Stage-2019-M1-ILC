
#include <cmath>
#include <iostream>
#include "TVectorF.h"

{	
 	
	TCanvas *c3 = new TCanvas ("Total energy resolution in fonction of total measured energy (ZEUS)");
	

        TMultiGraph *F=new TMultiGraph();
        F->GetXaxis()->SetTitle("Total energy measured (GeV)");
	F->GetYaxis()->SetTitle("#sigma*#sigma");
        F->SetTitle("Square of the total resolution in fonction of total measured energy (ZEUS)");
	
	TNtuple *Data=new TNtuple("Data", "Data","Et:sigsqr:Eterr:sigerr");
	int n=Data->ReadFile("output12.txt");

	float Et,sigsqr,Eterr,sigerr;

	Data->SetBranchAddress("Et", &Et);
	Data->SetBranchAddress("sigsqr", &sigsqr);
	Data->SetBranchAddress("Eterr", &Eterr);
	Data->SetBranchAddress("sigerr", &sigerr);
	
	TVectorF EtV,sigV,EterrV,sigerrV;

		
	TGraphErrors *gr= new TGraphErrors("output12.txt");
	TGraphErrors *gsq= new TGraphErrors("output13.txt");
	//TGraphErrors *gf= new TGraphErrors("output10.txt");

	TF1 *plot = new TF1("plot", "[0]*x+[1]*x*x+[2]*x*x*x*x+[3]",0,3000);
        TF1 *ploti = new TF1("ploti","[0]*x+[1]*x*x+[2]*x*x*x*x+[3]",0,3000);
	//TF1 *plotu = new TF1("plotu","[0]*x+[1]*x*x+[2]*x*x*x*x+[3]",0,3000);

	gr->SetTitle("Fit with [0]*x+[1]*x*x+[2]*x*x*x*x+[3]");
	plot->SetLineColor(2);
	gr->SetLineColor(2);
	gr->SetMarkerStyle(22);
	gr->Fit("plot","0");

	gsq->SetTitle("Resolution determined for events without neutrinos");
	gsq->SetMarkerStyle(21);
	gsq->SetLineColor(3);
	ploti->SetLineColor(3);
	gsq->Fit("ploti");
	
	TGraphErrors *SigV= new TGraphErrors();
	SigV->SetTitle("Total resolution over E in fonction of total measured energy (ZEUS)");
	SigV->GetXaxis()->SetTitle("Total energy measured (GeV)");
	SigV->GetYaxis()->SetTitle("#sigma/E");

	for (int i=0; i<n; i++) {
		Data->GetEntry(i);
		SigV->SetPoint(i,Et,sqrt(sigsqr)/Et);
		SigV->SetPointError(i,Eterr,sigerr/(2*sqrt(sigsqr)*Et));
	}

	TF1 *SigLaw= new TF1("SigLaw", "(1/x)*sqrt([0]*x+[1]*x*x+[2]*x*x*x*x+[3])",0,3000);
	SigLaw->SetParameters(plot->GetParameter(0),plot->GetParameter(1),plot->GetParameter(2),plot->GetParameter(3));
	SigV->SetMarkerStyle(21);
	SigV->Fit("SigLaw");

	TF1 *SigLawSQR= new TF1("SigLawSQR", "[0]*x+[1]*x*x+[2]*x*x*x*x+[3]",0,3000);
	SigLawSQR->SetParameters(SigLaw->GetParameter(0),SigLaw->GetParameter(1),SigLaw->GetParameter(2),SigLaw->GetParameter(3));
	gr->Fit("SigLawSQR");
		

	cout<<"Total resolution law : #sigma/E= ";
	cout<<sqrt(SigLaw->GetParameter(0))<<"/sqrt(E) ";
	if (SigLaw->GetParameter(3)>0) {
	cout<<"*+* ";
	} else {cout<<"*-* ";}
	cout<<sqrt(abs(SigLaw->GetParameter(3)))<<"/E ";
	if (SigLaw->GetParameter(2)>0) {
	cout<<"*+* ";
	} else {cout<<"*-* ";}
	cout<<sqrt(abs(SigLaw->GetParameter(2)))<<"*E ";
	if (SigLaw->GetParameter(1)>0) {
	cout<<"*+* ";
	} else {cout<<"*-* ";}
	cout<<sqrt(SigLaw->GetParameter(1))<<endl;

	//gf->SetTitle("Resolution determined only with a CB fit");
	//gf->SetLineColor(4);
	//gf->SetMarkerStyle(20);
	//gf->Fit("plotu");
	//plotu->SetLineColor(4);

	
	F->Add(gr,"P"); //Star for each point
        F->Add(gsq,"P");
	//F->Add(gf,"P");
	c3->cd(1);
        F->Draw("AP");
        c3->BuildLegend(0.2, 0.2, .4, .5);
	gStyle->SetOptFit();
	//gr->GetYaxis()->SetRangeUser(0.,0.1);

	TCanvas *c4= new TCanvas ("Total energy resolution in fonction of total measured energy 2 (ZEUS)");

	c4->cd(1);
	SigV->Draw("AP");
	//SigLaw->Draw();
	
	//gf->Fit("ploti");
	//c3->Update();
		


	c3->SaveAs("/home/ilc/armatol/Bureau/TestCMS") 
	
	
}
	
