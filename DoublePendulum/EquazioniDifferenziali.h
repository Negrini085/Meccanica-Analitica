#ifndef __EquazioniDifferenziali_h__
#define __EquazioniDifferenziali_h__

#include "FunzioniVettoriali.h"
#include "OperazioniVector.h"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class EquazioneDifferenzialeBase{

	public:
	virtual vector<double> Passo(double t, const vector<double>& x, double h, FunzioneVettorialeBase& f) const =  0;

};

class RungeKutta4 : public EquazioneDifferenzialeBase{

	public:
	virtual vector<double> Passo(double t, const vector<double>& x, double h, FunzioneVettorialeBase& f) const override {
		vector<double> K1;
		vector<double> K2;
		vector<double> K3;
		vector<double> K4;

		K1 = f.Eval(t, x);
		K2 = f.Eval(t + h/2, x + h/2*K1);
		K3 = f.Eval(t + h/2, x + h/2*K2);
		K4 = f.Eval(t + h, x + h*K3);


		return x + (double)h/6 * (K1 + K2 * 2. + K3 * 2. + K4);
	}
};
#endif //__EquazioniDifferenziali_h__