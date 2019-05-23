{
 	
	TCanvas *c3 = new TCanvas ("Total energy resolution in fonction of total measured energy");


        TMultiGraph *F=new TMultiGraph();
        F->GetXaxis()->SetTitle("Total energy measured (GeV)");
	F->GetYaxis()->SetTitle("#sigma/Etot");
        F->SetTitle("Total energy resolution in fonction of total measured energy");

	TGraph *gr=new TGraph("output.txt");
	//TGraph *gf= new TGraph("output2.txt");

	TF1 *plot = new TF1("plot", "[0]*sqrt(1/x)+[1]*(1/x)+[2]*x*x+[3]",0,3000);
        TF1 *ploti = new TF1("ploti","[0]*sqrt(1/x)",0,3000);

        ploti->SetParameter(0,0.16);

	
	
	

	gr->SetTitle("CALICE resolution fit");
	gr->SetLineColor(2);
	gr->Fit("plot");
	
	gr->SetMarkerStyle(22);
	


	/*gf->SetTitle("CMS resolution, w/ charged particles resolution");
	gf->SetMarkerStyle(21);
	gf->SetLineColor(3);
	gf->SetDrawOption("C");*/
	
	

	
	
	F->Add(gr); //Star for each point
       // F->Add(gf,"PC");
        F->Draw("AP");
        c3->BuildLegend(0.9, 0.7, .4, .4);
	//gr->GetYaxis()->SetRangeUser(0.,0.1);

	//gf->Fit("ploti");
	//c3->Update();
		
	gStyle->SetOptFit();

	c3->SaveAs("/home/ilc/armatol/Bureau/TestCMS")
	
	
}
	
