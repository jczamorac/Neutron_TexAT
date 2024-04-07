#include </home/jcz/proyectos/rootscope_develop/ROOTSCOPE/myROOTSCOPE.C>


void testfit(){


	//TFile* file1 = TFile::Open("HistoTeflon_TTL10_afterR35.root");
	//TFile* file1 = TFile::Open("suma.root");
	TFile* file1 = TFile::Open("tamu_sim_elas_Xsec.root");
	//TFile* file1 = TFile::Open("HistoC12_TTL15_afterR35.root");
	//{..some of your code which creates TH2 or TH1 objects..}
	//TH1F* h_1d = (TH1F*) file1->Get("Calib/Z1Cal");
	TH1F* h_1d = (TH1F*) file1->Get("Det3/det3_th_122.5");
	//TH2F* h_2d = (TH2F*) file1->Get("Calib/T2");

	//TH1F* h_1d = (TH1F*) file1->Get("Calib/E1T");
	//TH1F* h_1d = (TH1F*) file1->Get("Ex/T1EX");
	//TH1F* h_1d = (TH1F*) file1->Get("Ex/Z8EX");


	ROOTSCOPE* app = new ROOTSCOPE( gClient->GetRoot() );
	app->AddHisto( h_1d);
	//app->AddHisto( h_2d);


}

/*
develop version
* type d
* write -> set bgopt (bg option)
* type B for estimate Background
* type d
* write -> set fitN (numerofpeaks)
* select range
* select centroids
* type G for fitting
*/
