#include "rootheader.h"
double myfun(double *xx, double *par)
{
    double E = xx[0];
    double c = par[0];
    double alpha = par[1];
    double mu = par[2];
    double sigma = par[3];
    double beta = par[4];
    double lambda = par[5];
    double pi = 3.141592654;

    double gaus = alpha/sigma/sqrt(2*pi)*exp(-0.5*(E-mu)*(E-mu)/sigma/sigma);
    double constant = beta/mu*(TMath::Erf((mu-E)/sqrt(2)/sigma)-TMath::Erf(-E/sqrt(2)/sigma));
    double exp_c1 = (beta)*lambda*exp((sigma*sigma*lambda*lambda+2*lambda*E)/2)/(exp(lambda*mu)-1);
    double exp1 = TMath::Erf((mu-E-sigma*sigma*lambda)/(TMath::Sqrt(2)*sigma))-TMath::Erf((-E-sigma*sigma*lambda)/sqrt(2)/sigma);

    return c*(gaus+exp_c1*exp1+constant);
}

void fit()
{
    gStyle->SetOptFit(1);
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.14);
    gStyle->SetStatH(0.2);
    TFile *file;
    TF1 *f = new TF1("f",myfun,0,20000,6);
	double *pars;
	pars = new double[6];

	TChain *t = new TChain("evt");
	//t->Add("/junofs/production/public/users/zhangfy/non-uniform/offline_J17v1r1-Pre1/Examples/Tutorial/share/cls/Co60/0_0/evt_8*root");
	t->Add("/junofs/production/public/users/zhangfy/non-uniform/offline_J17v1r1-Pre1/Examples/Tutorial/share/cls/Co60/10000_10000/evt_*root");
	
	t->Draw("totalPE>>h(200,0,0)","","",40000);
	TH1F *h = (TH1F*)gDirectory->Get("h");

	int maxbin = h->GetMaximumBin();
	double maxbincenter = h->GetBinCenter(maxbin);

	t->Draw(Form("totalPE>>h(100,%i,%i)",int(maxbincenter-100),int(maxbincenter+100)),"","",40000);	
	h = (TH1F*)gDirectory->Get("h");
	//h->Fit("gaus");
	
	
//	h->Fit("f");

/*  
    file = new TFile("1.022_AfterC.root","read");
  	t = (TTree*)file->Get("tt");

	TF1 *f = new TF1("f",myfun,0,20000,6);
	t->Draw("newPE>>h(200,1000,1600)","edep>1.0218");	
	TH1F *h = (TH1F*)gDirectory->Get("h");
*/	h->Fit("gaus","Q");
    TF1* fun = h->GetFunction("gaus");

    double C = fun->GetParameter(0);
    double mean = fun->GetParameter(1);
    double emean = fun->GetParError(2);
    double sigma = fun->GetParameter(2);

	//f->SetParameters(C*sqrt(2*3.14159)*sigma,0.99,mean,sigma,0.001,0.01,0.02);
	pars[0] = C*sqrt(2*3.14159)*sigma;
	pars[1] = 0.99;
	pars[2] = mean;
	pars[3] = sigma;
	pars[4] = 0.001;
	pars[5] = 0.01;
	f->SetParameters(pars);

	f->SetParNames("C","#alpha","#mu","#sigma","#beta","#lambda");
	t->Draw(Form("totalPE>>h(%i,%i,%i)",int(5*sigma),int(mean-6*sigma),int(mean-6*sigma)+int(5*sigma)*2),"","",40000);	
	h = (TH1F*)gDirectory->Get("h");
	h->Fit(f);

}
