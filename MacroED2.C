
#include <cmath>
#include <iostream>

{	

	/*TFile f1("Resoluncer.root");
		  
  	TH1F *H1 = (TH1F*)f1.Get("DistribRes;2");

	cout<<H1->GetMean()<<endl;*/

 	
	TCanvas *c3 = new TCanvas ("Total energy resolution in fonction of total measured energy");
	//c3->Divide(2);

        TMultiGraph *F=new TMultiGraph();
        F->GetXaxis()->SetTitle("Total energy measured (GeV)");
	F->GetYaxis()->SetTitle("#sigma²");
        F->SetTitle("Square of the total resolution in fonction of total measured energy");

	TGraphErrors *gr= new TGraphErrors("output.txt");
	TGraph *gsq= new TGraph("output4.txt");
	//TGraph *gf= new TGraph("output5.txt");

	TF1 *plot = new TF1("plot", "[0]*x+[1]*x*x+[2]*x*x*x*x+[3]",0,3000);
        TF1 *ploti = new TF1("ploti","[0]*x+[1]*x*x+[2]*x*x*x*x+[3]",0,3000);

	gr->SetTitle("Total energy measured fit with a crystall-ball");
	plot->SetLineColor(2);
	gr->SetLineColor(2);
	gr->SetMarkerStyle(22);
	gr->Fit("plot");

	gsq->SetTitle("Total energy measured fit with a gaussian (events without neutrinos)");
	gsq->SetMarkerStyle(21);
	gsq->SetLineColor(3);
	ploti->SetLineColor(3);
	gsq->Fit("ploti");
	
	/*
	gf->SetTitle("No CP resolution");
	gf->SetLineColor(4);

	gsq->Fit("ploti");
	*/

	F->Add(gr,"P"); //Star for each point
        F->Add(gsq,"P");
	//F->Add(gf,"PC*");
        F->Draw("AP");
        c3->BuildLegend(0.6, 0.6, .4, .4);
	//gr->GetYaxis()->SetRangeUser(0.,0.1);

	//gf->Fit("ploti");
	//c3->Update();
		
	gStyle->SetOptFit();

	c3->SaveAs("/home/ilc/armatol/Bureau/TestCMS") 
	
	
}
	
