#include "EquazioniDifferenziali.h"
#include "FunzioniVettoriali.h"
#include "OperazioniVector.h"

#include <iostream>
#include <iomanip>
#include <cmath>

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH2F.h"

int main(int argc, char** argv) {

	if(argc != 3) {
		cout << "Funzionamento programma: " << argv[0] << "<tempo massimo> <numero di step> " << endl;
		exit(-1);
	}
	RungeKutta4 lavora;
  DoppioPendolo pendolo1(1, 1, 5, 4);
	DoppioPendolo pendolo2(1, 1, 5, 4);
	vector<double> stato1 {1.001, 0, 0, -5};
	vector<double> stato2 {1, 0, 0, -5};

	TApplication myApp("myApp", 0 , 0);
	TGraph pos1_1; TGraph pos1_2; TGraph pos2_1; TGraph pos2_2;
	TGraph pend1; TGraph pend2; TGraph zoomg;

	double tmax = atof(argv[1]);
	int nstep = atoi(argv[2]);
	double h = tmax/nstep;
	double t = 0;
	//Posizioni delle due masse del primo pendolo
	double x1_1 = pendolo1.GetL1()*sin(stato1[0]);
	double y1_1 = -pendolo1.GetL1()*cos(stato1[0]);
	double x1_2 = x1_1 + pendolo1.GetL2()*sin(stato1[1]);
	double y1_2 = y1_1 - pendolo1.GetL2()*cos(stato1[1]);
	//Posizioni delle due masse del secondo pendolo
	double x2_1 = pendolo2.GetL1()*sin(stato2[0]);
	double y2_1 = -pendolo2.GetL1()*cos(stato2[0]);
	double x2_2 = x2_1 + pendolo2.GetL2()*sin(stato2[1]);
	double y2_2 = y2_1 - pendolo2.GetL2()*cos(stato2[1]);

	TCanvas* c1 = new TCanvas("c1", "Moto in evoluzione");
	c1 -> Divide(1,1);
	c1 -> cd(1);
	//Setto punti necessari
	pos1_1.SetPoint(0, x1_1, y1_1); pos1_2.SetPoint(0, x1_2, y1_2);
	pos2_1.SetPoint(0, x2_1, y2_1); pos2_2.SetPoint(0, x2_2, y2_2);	
	pend1.SetPoint(0, 0, 0); pend1.SetPoint(1, x1_1, y1_1); pend1.SetPoint(2, x1_2, y1_2);
	pend2.SetPoint(0, 0, 0); pend2.SetPoint(1, x2_1, y2_1); pend2.SetPoint(2, x2_2, y2_2);
	zoomg.SetPoint(0, 0, 8); zoomg.SetPoint(1, 8, -5 ); zoomg.SetPoint(2, -8, -9);

	
	//Disegno i pendoli - triettorie
	zoomg.SetMarkerColor(10); zoomg.SetLineColor(10);
	pend1.SetMarkerColor(6); pend1.SetLineColor(6);
	pend2.SetMarkerColor(4); pend2.SetLineColor(4);
	pos1_2.SetMarkerColor(30); pos2_2.SetMarkerColor(50);
	pos1_2.SetLineColor(30); pos2_2.SetLineColor(50);	
	zoomg.Draw("ALP");
	pos1_2.Draw("sameCP"); pos2_2.Draw("sameCP");
	//pos1_1.Draw("sameCP"); pos2_1.Draw("sameCP");
	pend1.Draw("*sameL"); pend2.Draw("*sameL");

	c1 -> Modified();
	c1 -> Update();	


	for(int i=1; i<=nstep; i++){
		//Evolvo il sistema
		stato1 = lavora.Passo(t, stato1, h, pendolo1);
		stato2 = lavora.Passo(t, stato2, h, pendolo2);
		//Calcolo le coordinate delle masse del primo pendolo
		x1_1 = pendolo1.GetL1()*sin(stato1[0]);
		y1_1 = -pendolo1.GetL1()*cos(stato1[0]);
		x1_2 = x1_1 + pendolo1.GetL2()*sin(stato1[1]);
		y1_2 = y1_1 - pendolo1.GetL2()*cos(stato1[1]);
		//Calcolo le coordinate delle masse del secondo pendolo
		x2_1 = pendolo2.GetL1()*sin(stato2[0]);
		y2_1 = -pendolo2.GetL1()*cos(stato2[0]);
		x2_2 = x2_1 + pendolo2.GetL2()*sin(stato2[1]);
		y2_2 = y2_1 - pendolo2.GetL2()*cos(stato2[1]);
		
		//Setto i vari punti necessari
		pos1_1.SetPoint(i, x1_1, y1_1); pos1_2.SetPoint(i, x1_2, y1_2);
		pos2_1.SetPoint(i, x2_1, y2_1); pos2_2.SetPoint(i, x2_2, y2_2);	
		pend1.SetPoint(0, 0, 0); pend1.SetPoint(1, x1_1, y1_1); pend1.SetPoint(2, x1_2, y1_2);
		pend2.SetPoint(0, 0, 0); pend2.SetPoint(1, x2_1, y2_1); pend2.SetPoint(2, x2_2, y2_2);

		//Disegno i pendoli - traiettorie
		zoomg.Draw("ALP");
		pos1_2.Draw("sameCP"); pos2_2.Draw("sameCP");
		//pos1_1.Draw("sameCP"); pos2_1.Draw("sameCP");
		pend1.Draw("*sameL"); pend2.Draw("*sameL");
		c1 -> Modified();
		c1 -> Update();	
		
	}
	myApp.Run();
	
return 0;
} 