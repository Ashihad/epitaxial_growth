/**********************************************************************************************************
 * 
 * program symuluje wzrost kryształu 1+1
 * oddziaływanie pomiedzy atomami opisywane jest w przybliżeniu parabolicznym
 * po przesunieciu atomu relaksowana jest lokalnie sieć w celu zmniejszenia napręzeń
 * 
 * atomy początkowo są rozmieszczone w węzłach siatki, w wyniku naprężeń mogą się nieco przesunąć poza węzeł
 **********************************************************************************************************/

#include "sim.hpp"

using namespace std;

int main()
{
	/***************************************************************************************************
	 * jednostki:   
	 * 			dlugosc - [Angstrem], 
	 * 			energia - [eV]
	 * 			czas -    [sekunda]
	 * 			temperatura - [K]	
	 ***************************************************************************************************/

	double as {5.43/2}; //Si [angstrem]  - podloze
	double ag {5.66/2}; //Ge [angstrem]  - atomy deponowane
	
	double mu {(ag-as)/as};
 
	double skl {15.85/as/as};  //[eV/angstrem^2] - stala sprzezystosci
	double skd {skl/2.};       //[Ha/ab^2]
	
	double al {ag + as*mu*skd/(skl+skd)}; //vertical lattice spacing
	
	double D0 {3.83E+13};    //[angstrem^2/sekunda] - fitting factor
	double D {D0/as/as};     //[1/sekunda]
	
	double gamma {0.4};      // [eV] - bond energy
	double E {0.53};         //[eV] - fitting factor for experimental data
	  
	std::size_t nx {500};
	std::size_t ny {80};
	
	
	std::size_t nsurf {32}; //wysokosc podloza Si
	
	std::size_t irange_min {10}; // w zakresie (xi,yj)+/- irange liczne sa lokalne zmiany polozen atomow
	std::size_t irange_max {10}; //maksymalny rozmiar otoczenia w relaksacji lokalnej
	double tol_local_max {1.0E-2}; //tolerancja bledu w relaksacji lokalnej 
	
	std::size_t k_max_step {8}; //maksymalny zasieg dyfuzji w jednym kroku
	
	double f_deposition {0.08}; // number of layers/second 
	double temperature {600.}; //temperatura [K]
	
	int idiffusion {2}; //1-dyfuzja atomow podloza i deponowanych, 2-dyfuzja tylko atomow deponowanych (podloze stabilne)
	
	
	double time_fluence {8.0};
	double time_simulation {100.0};
	
	long unsigned iterations {1000}; //maksymalna liczba iteracji w trakcie rozwiazywania ukladu rownan
	double tolerance {1.0E-3}; //tolerancja bledu w iteracyjnej metodzie rozwiazywania ukaldu rownan
	
	int n_write_files {1000};        // zapis danych do pliku co tyle iteracji
	int n_global_relaxation {1000}; // wykonanie relaksacji globalnej co tyle iteracji
		
	//symulacja
	growth_simulation(as,ag,al,skl,skd,mu,D,E,gamma,nx,ny,nsurf,irange_min,irange_max,tol_local_max,
				time_fluence,time_simulation,k_max_step,
				f_deposition,temperature,idiffusion,iterations,tolerance,
				n_write_files, n_global_relaxation	);
}// main
