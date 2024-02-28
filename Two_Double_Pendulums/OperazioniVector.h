#ifndef __OperazioniVector_h__
#define __OperazioniVector_h__


#include <algorithm>
#include <iostream>
#include <numeric>
#include <cassert>
#include <vector>

using namespace std;


/******************************************
           Somma fra due vector
******************************************/
template <typename T> vector<T> operator+(const vector<T> &a, const vector<T> &b){

	if(a.size() != b.size()){
		cout << "Somma fra due vettori di dimensioni differenti: esco!!" << endl;
		exit(-2);
	}

	vector<T> somma;
	for(int i=0; i<a.size(); i++){
		somma.push_back(a[i] + b[i]);
	}

return somma;
}


/******************************************
        Differenza fra due vector
******************************************/
template <typename T> vector<T> operator-(const vector<T>& a, const vector<T>& b){

	if(a.size() != b.size()){
		cout << "Differenza fra due vettori di dimensioni diverse: esco!!" << endl;
		exit(-2);
	}

	vector<T> diff;
	for(int i=0; i<a.size(); i++){
		diff.push_back(a[i] - b[i]);
	}

return diff;
}


/******************************************
        Overloading operator +=
******************************************/
template <typename T> vector<T>& operator+=(vector<T>& a, vector<T>& b){
	
	if(a.size() != b.size()){
		cout << "Somma fra due vettori di dimensioni differenti: esco!!" << endl;
		exit(-2);
	}

	for(int i=0; i<a.size(); i++){
		a[i] += b[i];
	}

return a;
}


/******************************************
        Overloading operator -=
******************************************/
template <typename T> vector<T>& operator-=(vector<T>& a, vector<T>& b){

	if(a.size() != b.size()){
		cout << "Differenza fra due vettori di dimensioni differenti: esco!!" << endl;
		exit(-2);
	}

	for(int i=0; i<a.size(); i++){
		a[i] -= b[i];
	}

return a;
}


/******************************************
        Prodotto per scalare (prima)
******************************************/
template <typename T> vector<T> operator*(T scal, const vector<T>& a){
	
	vector<T> risultato;
	for(int i=0; i<a.size(); i++){
		risultato.push_back(scal*a[i]);
	}

return risultato;
}


/******************************************
        Prodotto per scalare (dopo)
******************************************/
template <typename T> vector<T> operator*(const vector<T>& a ,T scal){
	
	vector<T> risultato;
	for(int i=0; i<a.size(); i++){
		risultato.push_back(scal*a[i]);
	}

return risultato;
}


/******************************************
        Divisione per scalare (dopo)
******************************************/
template <typename T> vector<T> operator/(const vector<T>& a, T scal){
	
	vector<T> risultato;
	for(int i=0; i<a.size(); i++){
		risultato.push_back(a[i]/scal);
	}

return risultato;
}


/******************************************
           Prodotto scalare
******************************************/
template <typename T> T operator*(const vector<T>& a, const vector<T>& b){

	if(a.size() != b.size()){
		cout << "Prodotto scalare fra vettori con dimensioni diverse: esco!!" << endl;
		exit(-2);
	}

	T appo=0;
	for(int i=0; i<a.size(); i++){
		appo += a[i] * b[i];
	}
	

return appo;	
}


/******************************************
        Prodotto fra int e vector
******************************************/
template <typename T> vector<T> operator*(int k, const vector<T>& a){

	vector<T> ris;
	for(int i=0; i<a.size(); i++) {ris.push_back(k*a[i]);}

return ris;
}

/******************************************
        Prodotto fra vector e int
******************************************/
template <typename T> vector<T> operator*(const vector<T>& a, int k){

	vector<T> ris;
	for(int i=0; i<a.size(); i++) {ris.push_back(k*a[i]);}

return ris;
}


/******************************************
           Stampa di un vector
******************************************/
template <typename T> void Stampa(vector<T>& a){

	cout << "Stampa vettore: ";
  for(auto i = a.begin(); i != a.end(); i++){
		cout << *i << " ";
	}
	cout << endl;
	
}


#endif //__OperazioniVector_h__