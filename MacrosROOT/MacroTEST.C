{

TF1 *CaliceRes=new TF1("ResCalice","sqrt((0.132/sqrt(x))*(0.132/sqrt(x))+(0.332/x)*(0.332/x)+0.0049*0.0049)",0,1000);
CaliceRes->GetXaxis()->SetTitle("Total energy measured");
CaliceRes->GetYaxis()->SetTitle("#sigma/Etot");
TF1 *PandoraPFA=new TF1("ResPFA","0.79*sqrt(((0.21*0.21)/x)+(0.007*0.007)+((0.00004*0.00004)*x*x)+((0.0001)*pow(x/100,0.6)))",0,1000);
PandoraPFA->SetLineColor(1);
TF1 *Morgunov=new TF1("MorgRes","0.14/sqrt(x)",0,1000);
Morgunov->SetLineColor(4);
PandoraPFA->SetTitle("PandoraPFA resolution");
Morgunov->SetTitle("Morgunov resolution");
CaliceRes->SetTitle("Simulated Calice resolution");

TF1 *CaliceResSQR=new TF1("ResCaliceSQR","x*x*((0.132/sqrt(x))*(0.132/sqrt(x))+(0.332/x)*(0.332/x)+(0.0049*0.0049))",0,1000);
TF1 *PandoraPFASQR=new TF1("ResPFASQR","(0.79*0.79)*x*x*(((0.21*0.21)/x)+(0.007*0.007)+((0.00004*0.00004)*x*x)+((0.0001)*pow(x/100,0.6)))",0,1000);
PandoraPFASQR->SetLineColor(3);


TF1 *Conf=new TF1("Conf","sqrt(abs(ResPFASQR-ResCaliceSQR))",0,1000);
Conf->SetTitle("Confusion factor in fonction of E");
Conf->GetXaxis()->SetTitle("Etot");
Conf->GetYaxis()->SetTitle("Confusion factor");

TCanvas c3;
CaliceRes->Draw();
PandoraPFA->Draw("SAME");
Morgunov->Draw("SAME");
c3.BuildLegend(0.2, 0.2, .4, .5);
//PandoraPFASQR->Draw();
//CaliceResSQR->Draw("SAME");
TCanvas c4;
Conf->Draw();
c4.BuildLegend(0.2, 0.2, .4, .5);
}
		
