{
 	
	TCanvas *c3 = new TCanvas ("Total energy resolution in fonction of total measured energy");


        TMultiGraph *F=new TMultiGraph();
        F->GetXaxis()->SetTitle("Stochastic factor of photons resolution");
	F->GetYaxis()->SetTitle("#sigma/Etot");
        F->SetTitle("Total energy resolution in fonction of photons resolution (ZEUS)");

	TGraph *gr=new TGraph("output31.txt");
	TGraph *gf=new TGraph("output32.txt");
        TGraph *gk=new TGraph("output33.txt");

       //TF1 *plot = new TF1("plot", "[0]*sqrt(1/x)",0,3000);
       //TF1 *ploti = new TF1("ploti","[0]*sqrt(1/x)+[1]*(1/x)+[2]",0,2500);

	
	
	

	gr->SetTitle("Ecdm=200 GeV");
	gr->SetLineColor(2);
	//gr->Fit("plot");
	gr->SetMarkerStyle(22);
	


	gf->SetTitle("Ecdm=350 GeV");
	gf->SetMarkerStyle(21);
	gf->SetLineColor(3);
	gf->SetDrawOption("C");
	
	
	gk->SetTitle("Ecdm=500 GeV");
	gk->SetMarkerStyle(23);
	gk->SetLineColor(4);
	gk->SetDrawOption("C");
	

	
	
	F->Add(gr,"PC"); //Star for each point
        F->Add(gf,"C*");
        F->Add(gk,"PC");
        F->Draw("AP");
        c3->BuildLegend(0.5, 0.5, .3, .6);
	//gr->GetYaxis()->SetRangeUser(0.,0.1);

	//gf->Fit("ploti");
	//c3->Update();
		
	//gStyle->SetOptFit();

	c3->SaveAs("/home/ilc/armatol/Bureau/TestCMS")
	
	
}
	
