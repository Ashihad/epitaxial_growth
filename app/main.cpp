/**********************************************************************************************************
 * 
 * program symuluje wzrost kryształu 1+1
 * oddziaływanie pomiedzy atomami opisywane jest w przybliżeniu parabolicznym
 * po przesunieciu atomu relaksowana jest lokalnie sieć w celu zmniejszenia napręzeń
 * 
 * atomy początkowo są rozmieszczone w węzłach siatki, w wyniku naprężeń mogą się nieco przesunąć poza węzeł
 **********************************************************************************************************/

#include "sim.hpp"

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<vector>
#include<random>
#include<ctime>
#include<chrono>

// #include"mkl.h"
// #include<omp.h>

// #include <Eigen/Sparse>

using namespace std;

int main(){
	
	
	int nx,ny,nsurf;
	double as,ag,al,skl,skd,D0,D,E;
	double gamma,E0; 
	double ab, ha;
	double mu;
	double time_fluence;
	double time_simulation;
	int k_max_step;
	double f_deposition;
	double temperature;
	int idiffusion;
	int iterations;
	double tolerance;
	int irange_min;
	int irange_max;
	double tol_local_max;
	
	int n_write_files;
	int n_global_relaxation;
	
	
	// omp_set_num_threads(1); //ograniczamy liczbe watkow
	//test_gen_discrete_sign();
	
	/***************************************************************************************************
	 * jednostki:   
	 * 			dlugosc - [Angstrem], 
	 * 			energia - [eV]
	 * 			czas -    [sekunda]
	 * 			temperatura - [K]	
	 ***************************************************************************************************/

	as=5.43/2; //Si [angstrem]  - podloze
	ag=5.66/2; //Ge [angstrem]  - atomy deponowane
	
	mu=(ag-as)/as;
 
	skl=15.85/as/as;  //[eV/angstrem^2] - stala sprzezystosci
	skd=skl/2.;       //[Ha/ab^2]
	
	al=ag+as*mu*skd/(skl+skd); //vertical lattice spacing
	
	D0=3.83E+13;    //[angstrem^2/sekunda] - fitting factor
	D=D0/as/as;     //[1/sekunda]
	
	gamma=0.4;      // [eV] - bond energy
	E=0.53;         //[eV] - fitting factor for experimental data
	  
	nx=500;
	ny=80;
	
	
	nsurf=32; //wysokosc podloza Si
	
	irange_min=10; // w zakresie (xi,yj)+/- irange liczne sa lokalne zmiany polozen atomow
	irange_max=10; //maksymalny rozmiar otoczenia w relaksacji lokalnej
	tol_local_max=1.0E-2; //tolerancja bledu w relaksacji lokalnej 
	
	k_max_step=8; //maksymalny zasieg dyfuzji w jednym kroku
	
	f_deposition=0.08; // number of layers/second 
	temperature=600.; //temperatura [K]
	
	idiffusion=2; //1-dyfuzja atomow podloza i deponowanych, 2-dyfuzja tylko atomow deponowanych (podloze stabilne)
	
	
	time_fluence=8.0;
	time_simulation=100.0;
	
	iterations=1000; //maksymalna liczba iteracji w trakcie rozwiazywania ukladu rownan
	tolerance=1.0E-3; //tolerancja bledu w iteracyjnej metodzie rozwiazywania ukaldu rownan
	
	n_write_files=1000;        // zapis danych do pliku co tyle iteracji
	n_global_relaxation=1000; // wykonanie relaksacji globalnej co tyle iteracji
		
	//symulacja
	growth_simulation(as,ag,al,skl,skd,mu,D,E,gamma,nx,ny,nsurf,irange_min,irange_max,tol_local_max,
				time_fluence,time_simulation,k_max_step,
				f_deposition,temperature,idiffusion,iterations,tolerance,
				n_write_files, n_global_relaxation	);
}// main
