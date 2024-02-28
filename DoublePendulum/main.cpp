#include "EquazioniDifferenziali.h"
#include "FunzioniVettoriali.h"
#include "OperazioniVector.h"

#include <iostream>
#include <iomanip>
#include <vector>

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

using namespace std;

int main(int argc, char** argv){
	
	if(argc != 3){
		cout << "Funzionamento programma: " << argv[0] << " <tempo massimo> <numero step>" << endl;
		exit(-1);
	}

	RungeKutta4 lavora;
	DoppioPendolo test(1, 1, 5, 4);
	vector<double> stato {1.6, 0, 0, 1};
	double tmax = atof(argv[1]);
	int nstep = atoi(argv[2]);
	double h = tmax/nstep;
	double t = 0;
	double x1 = test.GetL1()*sin(stato[0]);
	double y1 = -test.GetL1()*cos(stato[0]);
	double x2 = x1 + test.GetL2()*sin(stato[1]);
	double y2 = y1 - test.GetL2()*cos(stato[1]);

	//Cose per root
	TApplication myApp("myApp", 0, 0);
	TGraph traj1; TGraph traj2;
	TGraph pendolo;
	TCanvas* c1 = new TCanvas("c1", "Moto in evoluzione");
	c1 -> Divide(1,1);
	c1 -> cd(1);
	//Setto punti necessari
	traj1.SetPoint(0, x1, y1); traj2.SetPoint(0, x2, y2);
	pendolo.SetPoint(0, 0, 0); pendolo.SetPoint(1, x1, y1); pendolo.SetPoint(2, x2, y2);

	traj2.SetTitle("Traiettoria; X; Y");
	traj2.GetXaxis()->SetLimits(-10, 10);

	//Disegno le traiettorie
	traj2.Draw("ACP"); traj1.SetMarkerColor(7);
	traj1.SetLineColor(7); traj1.Draw("sameCP");
	//Disegno il pendolo
	pendolo.SetMarkerColor(6); pendolo.SetLineColor(6);
		pendolo.Draw("*sameL");
	c1 -> Modified();
	c1 -> Update();	


	for(int i=1; i<=nstep; i++){
		//Evolvo il sistema
		stato = lavora.Passo(t, stato, h, test);
		//Calcolo le coordinate
		x1 = test.GetL1()*sin(stato[0]);
		y1 = -test.GetL1()*cos(stato[0]);
		x2 = x1 + test.GetL2()*sin(stato[1]);
		y2 = y1 - test.GetL2()*cos(stato[1]);
		//Setto i vari punti necessari
		traj1.SetPoint(i, x1, y1); traj2.SetPoint(i, x2, y2);
		pendolo.SetPoint(0, 0, 0); pendolo.SetPoint(1, x1, y1); pendolo.SetPoint(2, x2, y2);
		//Disegno le traiettorie
		traj2.Draw("ACP"); traj1.SetMarkerColor(7);
		traj1.SetLineColor(7); traj1.Draw("sameCP");
		//Disegno il pendolo
		pendolo.SetMarkerColor(6); pendolo.SetLineColor(6);
		pendolo.Draw("*sameL");
		c1 -> Modified();
		c1 -> Update();
	}

	
	
	myApp.Run();
return 0;
}