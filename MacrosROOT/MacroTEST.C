{

TF1 *CaliceRes=new TF1("ResCalice","sqrt((0.132/sqrt(x))*(0.132/sqrt(x))+(0.332/x)*(0.332/x)+0.0049*0.0049)",0,1000);
TF1 *PandoraPFA=new TF1("ResPFA","0.79*sqrt(((0.21*0.21)/x)+(0.007*0.007)+((0.00004*0.00004)*x*x)+((0.0001)*pow(x/100,0.6)))",0,1000);
PandoraPFA->SetLineColor(1);
TF1 *Morgunov=new TF1("MorgRes","0.14/sqrt(x)",0,1000);
Morgunov->SetLineColor(4);



TF1 *CaliceResSQR=new TF1("ResCaliceSQR","x*x*((0.132/sqrt(x))*(0.132/sqrt(x))+(0.332/x)*(0.332/x)+(0.0049*0.0049))",0,1000);
TF1 *PandoraPFASQR=new TF1("ResPFASQR","(0.79*0.79)*x*x*(((0.21*0.21)/x)+(0.007*0.007)+((0.00004*0.00004)*x*x)+((0.0001)*pow(x/100,0.6)))",0,1000);
PandoraPFASQR->SetLineColor(3);


TF1 *Conf=new TF1("Conf","sqrt(abs(ResPFASQR-ResCaliceSQR))",0,1000);

TCanvas c3;
c3.Divide(2);
c3.cd(1);
CaliceRes->Draw();
PandoraPFA->Draw("SAME");
Morgunov->Draw("SAME");
//PandoraPFASQR->Draw();
//CaliceResSQR->Draw("SAME");
c3.cd(2);
Conf->Draw();
}
		
