#ifndef __FunzioniVettoriali_h__
#define __FunzioniVettoriali_h__

#include <iostream>
#include <vector>
#include <cmath>

#define USE_MATH_NDEFINES
#define acc_g 9.81

using namespace std;

class FunzioneVettorialeBase{

	public:
	virtual vector<double> Eval(double t, const vector<double>& stato) const = 0;

};

class DoppioPendolo : public FunzioneVettorialeBase{

	public:
	//Cotruttore di default
	DoppioPendolo() {m_m1=0; m_m2=0; m_l1=0; m_l2=0; m_k=0;}
	//Cotruttore che accetta in ingresso masse e lunghezze
	DoppioPendolo(double m1, double m2, double l1, double l2, double k) {
		m_m1 = m1; m_m2 = m2; m_l1 = l1; m_l2 = l2; m_k = k;
	}
	//Metodo per conoscere m_l1
	double GetL1() { return m_l1; }
	//Metodo per modificare m_l1
	void SetL1(double tmp) { m_l1 = tmp; }
	//Metodo per conoscere m_l2
	double GetL2() { return m_l2; }
	//Metodo per modificare m_l2
	void SetL2(double tmp) { m_l2 = tmp; }

	vector<double> Eval(double t, const vector<double>& stato) const override {
		vector<double> derivate;
		double a, b, c, d, e, f;
		double rad;
		double x, y;

		//Calcolo la radice
		rad = sqrt(2*pow(m_l1, 2) + pow(m_l2, 2) + 2*m_l1*m_l2 - 2*m_l1*cos(stato[0])*(m_l1 + m_l2));

		//Calcolo i valori dei sei coefficienti
		a = (m_m1 + m_m2)*pow(m_l1, 2);
		b = m_m2*m_l1*m_l2*cos(stato[0] - stato[1]);
		c = m_m2*m_l1*m_l2*stato[3]*sin((stato[0] - stato[1]))*(stato[2] - stato[3]) - m_m2*m_l1*m_l2*stato[2]*stato[3]*sin(stato[0] - stato[1]) - (m_m1 + m_m2)*acc_g*m_l1*sin(stato[0]) - m_k*(rad - 4)*m_l1*sin(stato[0])*(m_l1 + m_l2)/rad;
		d = m_m2*m_l1*m_l2*cos(stato[0] - stato[1]);
		e = m_m2 * pow(m_l2, 2);
		f = m_m2*m_l1*m_l2*stato[2]*(sin(stato[0] - stato[1]))*(stato[2] - stato[3]) + m_m2*m_l1*m_l2*stato[2]*stato[3]*sin(stato[0] - stato[1]) - m_m2*acc_g*m_l2*sin(stato[1]);

		y = (a*f - d*c)/(e*a - d*b);
		x = c/a - (b/a)*y;

		derivate.push_back(stato[2]);
		derivate.push_back(stato[3]);
		derivate.push_back(x);
		derivate.push_back(y);

		/*double energia = 0.5*(m_m1 + m_m2)*pow(m_l1*stato[2], 2) + 0.5*m_m2*pow(m_l2*stato[3], 2) + m_m2*m_l1*m_l2*cos(stato[0] - stato[1])*stato[2]*stato[3] - (m_m1 + m_m2)*acc_g*m_l1*cos(stato[0]) - m_m2*acc_g*m_l2*cos(stato[1]) + 0.5*m_k*pow(rad - 4, 2);
		cout << "Energia sistema: " << energia << endl;*/

	return derivate;
	}

	protected:
	double m_m1, m_m2, m_l1, m_l2, m_k;
};
#endif //__FunzioniVettoriali_h__