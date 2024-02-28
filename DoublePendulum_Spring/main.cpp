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
	DoppioPendolo test(1, 1, 5, 4, 10);
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
	TGraph traj; TGraph pendolo; TGraph zoomg; TGraph molla;
	TCanvas* c1 = new TCanvas("c1", "Moto in evoluzione");
	c1 -> Divide(1,1);
	c1 -> cd(1);
	//Setto punti necessari
	traj.SetTitle("Traiettoria; X; Y");
	pendolo.SetPoint(0, 0, 0); pendolo.SetPoint(1, x1, y1); pendolo.SetPoint(2, x2, y2);
	zoomg.SetPoint(0, 3, 0); zoomg.SetPoint(1, -8, -5); zoomg.SetPoint(2, 8, -9);	
	molla.SetPoint(0, 0, -9); molla.SetPoint(1, x1, y1);
	traj.SetPoint(0, x2, y2);

	//Disegno traiettoria - pendolo - zoomg - molla
	zoomg.SetMarkerColor(10); zoomg.SetLineColor(10); zoomg.Draw("ALP");
	traj.SetMarkerColor(50); traj.SetLineColor(50);	traj.Draw("sameCP");
	pendolo.SetMarkerColor(6); pendolo.SetLineColor(6); pendolo.Draw("*sameL");
	molla.Draw("*sameL");
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

		//Setto punti necessari

		pendolo.SetPoint(0, 0, 0); pendolo.SetPoint(1, x1, y1); pendolo.SetPoint(2, x2, y2);
		zoomg.SetPoint(0, 0, 0); zoomg.SetPoint(1, -8, -5); zoomg.SetPoint(2, 8, -9);	
		molla.SetPoint(1, x1, y1); 	traj.SetPoint(i, x2, y2);
		
		//Disegno traiettoria - pendolo - zoomg - molla
		zoomg.Draw("ALP"); traj.Draw("sameCP"); pendolo.Draw("*sameL"); molla.Draw("*sameL");
		c1 -> Modified();
		c1 -> Update();	
	
	}

	
	
	myApp.Run();
return 0;
}