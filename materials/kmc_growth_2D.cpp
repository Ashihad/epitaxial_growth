/**********************************************************************************************************
 * 
 * program symuluje wzrost krysztalu 1+1
 * oddzia�ywanie pmiedzy atomami opisywane jest w przyblizeniu parabolicznym
 * po przesunieciu atomu relaksowana jest lokalnie siec w celu zmniejszenia naprezen
 * 
 * atomy poczatkowo sa rozmieszczone w wezlach siatki, w wyniku naprezen moga sie nieco przesunac poza wezel
 **********************************************************************************************************/

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

void solve_linear_system_u_v(vector<vector<vector<double>>> & crystal,int imin, int i_nodes,int jmin,int jmax,
				     double as, double ag, double al, double skl, double skd,int *iterations, double *tolerance, 
				     int ierr, double *bmax);

void compute_sparse_Ax_y(const int & n, double * acsr, int * icsr, int * jcsr, double * x, double * y);
double scalar_product_x_y(const int & n, double * x, double * y);

void solve_linear_system_CG_standard(const int & n, double * acsr, int * icsr, int  * jcsr, 
									  double * b, double * x, int * itmax, double *tol);


void conduct_relaxation(const int & ipos, const int & jpos, const int & irange_min, const int & irange_max, 
				int * iterations, double * tolerance,
				vector<vector<vector<double>>> & crystal, double tol_local_max, int iloc,
				double as,double ag, double al, double skl, double skd);



void compute_u_v_from_wxx(const int & number, const int & k,const int & i, const int & j, const int & nx, 
				  const double & skl, const double & skd, 
				  const vector<vector<int>> & ip,
				  const vector<vector<int>> & iboundary,
				  const vector<vector<double>> & d,
				  const vector<vector<vector<double>>> & crystal,
				  vector<double> & acol,	
				  vector<int> & jcol,
				  double * ff);


void compute_u_v_from_wxy(const int & number, const int & k, const int & i, const int & j, const int & nx, const double & skd,
				  const vector<vector<int>> & ip,
				  const vector<vector<int>> & iboundary,
				  const vector<vector<double>> & d,
				  const vector<vector<vector<double>>> & crystal,
				  vector<double> & acol,	
				  vector<int> & jcol,
				  double * ff);

void sort_and_add_matrix_elements(const int & nrow, const int & k, vector<int> & jcol, vector<double> & acol, 
					    double * acsr, int * icsr, int * jcsr);

void growth_simulation(double as, double ag, double al, double skl, double skd, double mu, double D, double E, double gamma,
				int nx, int ny, int nsurf, int irange_min, int irange_max, double tol_local_max, 
			     double time_fluence, double time_simulation, int k_max_step,
			     double f_deposition, double temperature, int idiffusion, const int iterations, const double tolerance,
			     int n_write_files, int n_global_relaxation	);

double compute_elastic_energy_wij(const int & i, const int & j, const vector<vector<vector<double>>> & crystal, 
				   const double & as,const double & ag,const double & al,const double & skl,const double & skd);


double compute_gradient(int nx, int ny,  vector<vector<vector<double>>> & crystal, 
				   const double & as,const double & ag,const double & al,const double & skl,const double & skd);


double gen_uniform();
int gen_discrete_1_K(const int & K);
int gen_sign();
int gen_discrete_1_K_multiply_sign(const int & K);
void test_gen_discrete_sign();

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







/********************************************************************************************************************************
 ******************************************************************************************************************************** 
 * 
 * 				procedura symulacji wzrostu krzystalu: 1+1 
 * 
 * 
 ******************************************************************************************************************************** 
 ********************************************************************************************************************************/

void growth_simulation(double as, double ag, double al, double skl, double skd, double mu, double D, double E, double gamma,
				int nx, int ny, int nsurf, int irange_min, int irange_max, double tol_local_max, 
			     double time_fluence, double time_simulation, int k_max_step,
			     double f_deposition, double temperature, int idiffusion,const int iterations, const double tolerance,
			     int n_write_files, int n_global_relaxation)
{
	
	printf("=============================================================== \n\n");
	printf("as [A] =   %15.5f \n",as);
	printf("ag [A] =   %15.5f \n",ag);
	printf("al [A] =   %15.5f \n",al);
	printf("skl [eV/A/A] =  %15.5f \n",skl);
	printf("skd [eV/A/A] =  %15.5f \n",skd);
	printf("mu [A] =   %15.5f \n",mu);
	printf("D [1/s]=    %15.5E \n",D);
	printf("E [eV] =    %15.5f \n",E);
	printf("gamma [eV] = %15.5f \n",gamma);
	printf("nx=    %d \n",nx);
	printf("ny=    %d \n",ny);
	printf("nsurf= %d \n",nsurf);
	printf("irange_min = %d \n",irange_min);
	printf("irange_max = %d \n",irange_max);
	printf("time_fluence [s] = %15.3f \n",time_fluence);
	printf("time_simulation [s] = %15.3f \n",time_simulation);
	printf("k_max_step= %d \n",k_max_step);
	printf("f_deposition [ML/s]= %15.3f \n",f_deposition);
	printf("temperature [K] = %15.0f \n",temperature);
	printf("================================================================ \n\n");
	
	
	/***************************************************************************************************
	 *  tablice z polozeniami atomow - polozenie (x,y)
	 *  odleglosci atomowe zrenormalizowane - oddzialywania/naprezenia skalowane wzgledem typu atomow
	 * 
	 ***************************************************************************************************/
	vector<vector<vector<double>>> crystal;
	crystal.resize(nx,vector<vector<double>>(ny,vector<double>(10,0.)));
	
	vector<vector<vector<double>>> crystal_copy;
	crystal_copy.resize(nx,vector<vector<double>>(ny,vector<double>(10,0.)));
	
	
	// ustawiamy atomy podloza
	for(int i=0;i<nx;i++){
		for(int j=0;j<nsurf;j++){
			crystal[i][j][0]=1.0; //0-empty, 1-Si, 2-Ge
			crystal[i][j][1]=0.0; //u
			crystal[i][j][2]=0.0; //v
		}
	}	
	
	
	int istart=0;//(nx/2)-10;
	int ikoniec=nx-1;//(nx/2+10);
	
	// kladziemy 1 monowarstwe atomow typu-2
	for(int i=istart;i<=ikoniec;i++){
		int j=nsurf;
	
	
		
		crystal[i][j][0]=2.0; //0-empty, 1-Si, 2-Ge
		crystal[i][j][1]=0.0; //u
		crystal[i][j][2]=0.0; //v
		
	
		
		j=nsurf+1;
		crystal[i][j][0]=2.0; //0-empty, 1-Si, 2-Ge
		crystal[i][j][1]=0.0; //u
		crystal[i][j][2]=0.0; //v
	
			
		
		j=nsurf+2;
		crystal[i][j][0]=2.0; //0-empty, 1-Si, 2-Ge
		crystal[i][j][1]=0.0; //u
		crystal[i][j][2]=0.0; //v
		
		
		j=nsurf+3;
		crystal[i][j][0]=2.0; //0-empty, 1-Si, 2-Ge
		crystal[i][j][1]=0.0; //u
		crystal[i][j][2]=0.0; //v
	
		j=nsurf+4;
		crystal[i][j][0]=2.0; //0-empty, 1-Si, 2-Ge
		crystal[i][j][1]=0.0; //u
		crystal[i][j][2]=0.0; //v
	/*
	*/
		
	}	
	
	
/*	
	
		FILE *fpp;
		fpp=fopen("atoms_zle_old.dat","r");
		for(int j=0;j<ny;j++){
			for(int i=0;i<nx;i++){
				double a0,a1,a2,a3,a4,a5,a6;
				fscanf(fpp,"%lf %lf  %lf  %lf  %lf  %lf  %lf ",&a0,&a1,&a2,&a3,&a4,&a5,&a6);
				crystal[i][j][0]=a0;
				crystal[i][j][1]=a1;
				crystal[i][j][2]=a2;
				crystal[i][j][3]=a3;
				crystal[i][j][4]=a4;
				crystal[i][j][5]=a5;
				crystal[i][j][6]=a6;
			}
		}
		fclose(fpp);
		printf("WCZYTANE\n");
		
*/		
	
	
	
	//ewolucja
	int iter;
	double time;
	double dt;
	double ab,ha;
	double kbt;
	double ev;
	double k_boltzmann;
	double fluence;
	double rdep;
	int n_at_diff;
	double ep,ri,r0;
	int neigh;
	int number_atoms=0;
	double n_proposed=1.;  //liczba wylosowanych procesow dyfuzji
	double n_accepted=1.; // liczba zaakceptowanych procesow dyfuzji
	
	ev=1.602E-19; //1-elektronowolt
	k_boltzmann=1.38E-23;
	kbt=k_boltzmann*temperature/ev; //energia termiczna w [eV]
	
	
	vector<vector<double>> atoms_diff; //tu trzymac bedziemy informacje o atomach podlegajacych dyfuzji
	atoms_diff.resize(nx+1,vector<double>(10,0));
	
	
	auto czas_start = std::chrono::high_resolution_clock::now();
	
	iter=0;	
	time=0.;
	while(time<time_simulation){
		
		iter++;	
		
		
		
	/******************************************************************************
	 * zapis do pliku
	 * 
	 ******************************************************************************/
		if(iter%n_write_files==0){	
			
			// liczymy energie sperezystosci dla kazdego atomu
			for(int i=0;i<nx;i++){
				for(int j=1;j<ny-1;j++){
						crystal[i][j][8]=compute_elastic_energy_wij(i,j,crystal,as,ag,al,skl,skd);
				}
			}
			
			FILE *fp;	
			fp=fopen("atoms.dat","w");	
			for(int i=0;i<nx;i++){
				for(int j=nsurf-1;j<ny;j++){
				if(crystal[i][j][0]>0.5)fprintf(fp,"%8d %8d  %8d    %15.5E\n",i,j,lround(crystal[i][j][0]),crystal[i][j][8]);
				}
			}	
			fclose(fp);
			
			
			double gradient_norm=compute_gradient(nx,ny,crystal,as,ag,al,skl,skd);
			
			
			fp=fopen("elastic_en.dat","w");	
			for(int i=0;i<nx;i++){
				for(int j=1;j<ny-1;j++){
					double grad_loc_2=pow(crystal[i][j][6],2)+pow(crystal[i][j][7],2);
				 fprintf(fp,"%8d %8d  %8d    %15.5E   %15.5E\n",i,j,lround(crystal[i][j][0]),crystal[i][j][8],grad_loc_2);
				}
				fprintf(fp,"\n");
			}	
			fclose(fp);
			
			
			//dane tymczasowe 
			FILE *fpp;
			fpp=fopen("atoms_data.dat","w");
			for(int j=0;j<ny;j++){
				for(int i=0;i<nx;i++){
					fprintf(fpp,"%15.5E  ",crystal[i][j][0]);
					fprintf(fpp,"%15.5E  ",crystal[i][j][1]);
					fprintf(fpp,"%15.5E  ",crystal[i][j][2]);
					fprintf(fpp,"%15.5E  ",crystal[i][j][3]);
					fprintf(fpp,"%15.5E  ",crystal[i][j][4]);
					fprintf(fpp,"%15.5E  ",crystal[i][j][5]);
					fprintf(fpp,"%15.5E  ",crystal[i][j][6]);
					fprintf(fpp,"\n");
				}
				
			}
			fclose(fpp);
			
		//********* INFO ************************************************	
			auto czas_teraz = std::chrono::high_resolution_clock::now();
			auto czas = std::chrono::duration_cast<std::chrono::milliseconds>(czas_teraz - czas_start);
			
			double duv_max=0.;
			for(int j=0;j<ny;j++){
				for(int i=0;i<nx;i++){
					double duv=pow(crystal[i][j][1],2)+pow(crystal[i][j][2],2);
					duv=sqrt(duv);
					if(duv>duv_max && crystal[i][j][0]>0.5)duv_max=duv;
				}
				
			}
			
			printf("iter, time, new_atoms, pr_acc, (real time [s]), duv_max =    %6d    %12.6f  %6d  %12.4f   (%12.1f)    %10.3f\n", 
				 iter,time,number_atoms,n_accepted/n_proposed,czas.count()*1.0E-3, duv_max);
			
			
		}
		
		
	/******************************************************************************
	 * globalna relaksacja naprezen w sieci
	 ******************************************************************************/		
		if(iter%n_global_relaxation==0 || iter==1){	
			int iloc=0; //relaksacja globalna
			int iter0=iterations;
			double tol0=tolerance;
			conduct_relaxation(0,0,irange_min,irange_max,&iter0,&tol0,crystal,tol_local_max,iloc,as,ag,al,skl,skd);			
		}
				
		
	/*********************************************************************************************************************
	 * okreslamy tempo depozycji atomow: deposition rate
	 * rate= fluence*nx -> tempo (prawdopodobienstwo/sekunde) depozycji jednego atomu w calym ukladzie (gdziekolwiek)
	 * 
	 *********************************************************************************************************************/
		if(time<time_fluence)fluence=f_deposition*nx; 
		else fluence=0.;
		rdep=(k_max_step+1.0)*(2*k_max_step+1.0)*fluence/6.; 
		atoms_diff[0][5]=rdep; //tempo depozycji atomow z wiazki
		atoms_diff[0][6]=rdep;
				
	/**********************************************************************************************************************
	 * szukamy atomow powierzchniowych podlegajacych dyfuzji -> tworzymy liste, z ktorej wybierzemy jeden lub depozycje
	 * idiffusion=1,2:   1-atomy podloza i zdeponowane, 2-tylko zdeponowane
	 * 
	 **********************************************************************************************************************/	
		n_at_diff=0; //liczba atomow powierzchniowych podlegajacych dyfuzji
		for(int i=0;i<nx;i++){
			for(int j=ny-1;j>=1;j--){
				int ll=lround(crystal[i][j][0]); //ll=0(brak), 1(podloze), 2(atom wiazki)
				if(ll>0){					
					if(ll>=idiffusion){
						n_at_diff++;
						atoms_diff[n_at_diff][0]=lround(crystal[i][j][0]); //typ atomu
						atoms_diff[n_at_diff][1]=i; //polozenie x
						atoms_diff[n_at_diff][2]=j; //polozenie y
						
						neigh=0; //sasiedzi: NN i NNN (NN=poziomo+pionowo, NNN=diagonala+antydiagonala)
						for(int ii=-1;ii<=1;ii++){
							for(int jj=-1;jj<=1;jj++){
								if((abs(ii)+abs(jj))>0){
									int i2=(i+ii+nx)%nx;
									int j2=j+jj;
									if(crystal[i2][j2][0]>0.5)neigh++; //jest sasiad - dodajemy	
								}
							}
						}
						atoms_diff[n_at_diff][3]=neigh; 
						ep=compute_elastic_energy_wij(i,j,crystal,as,ag,al,skl,skd); //local elastic energy
						atoms_diff[n_at_diff][4]=ep;
						
						r0=12.*D/(k_max_step+1.0)/(2*k_max_step+2.);
						double delta_w=ep; //neigh<=2
						if(neigh==3){
							delta_w=ep*1.5;
						}else if(neigh==4){
							delta_w=ep*2.0;
						}else if(neigh>=5){
							delta_w=ep*3.5;
						}
						ri=r0*exp((-gamma*neigh+delta_w+E)/kbt); //szacowane tempo dyfuzji atomu
						atoms_diff[n_at_diff][5]=ri; 
						atoms_diff[n_at_diff][6]=ri; //kopia dla porownania dokladnego prawdopodobienstwa
						atoms_diff[n_at_diff][7]=neigh;
						atoms_diff[n_at_diff][8]=delta_w;
					}
					break; //przerywamy sprawdzanie - natrafilismy na atom idac od gory
				}
			}
		}

		
		
		
	//krok czasowy -> zmiana czasu [odwrotnosc aktualnej sumy czestosci procesow dt=1/sum_{i}(Gamma_i)]
		double sum_ri=0;
		for(int i=0;i<=n_at_diff;i++){
			sum_ri+=atoms_diff[i][5];
		}
		dt=1.0/sum_ri;
		time+=dt;
	
		
	// losujemy proces: depozycja lub dyfuzja atomu	
		for(int i=1;i<=n_at_diff;i++){
			atoms_diff[i][5]+=atoms_diff[i-1][5]; //suma Gamma_i do generatora dyskretnego
		}
		
		double u1=gen_uniform()*sum_ri;
		int ktory=n_at_diff; //zabezpieczenie na wypadek gdyby petla nie zadziala
		for(int i=0;i<=n_at_diff;i++){
			if(u1<=atoms_diff[i][5]){
				ktory=i;	
				break;
			}
		}
		
			
		
		if(ktory==0){ 
		/***************************************************************************************************
		 * 			   ktory=0:  	losowa depozycja atomu
		 * 
		 ***************************************************************************************************/ 
			number_atoms++; //zwiekszamy liczbe atomow w ukladzie	

			int ipos=gen_discrete_1_K(nx)-1; //polozenie atomu w kierunku x: 0-(nx-1)
			int jpos=0;
			for(int j=ny-1;j>=1;j--){
				int ll=lround(crystal[ipos][j][0]);
				if(ll>0){
					jpos=j+1;
					if(jpos<ny){
						crystal[ipos][jpos][0]=2.;//dodajemy atom typu 2
						crystal[ipos][jpos][1]=0.;//u
						crystal[ipos][jpos][2]=0.;//v			
						break;
					}else{
						printf("za duzo atomow w tablicy: jpos=ny\n");
						printf("omijamy punkt i=%d\n",ipos);
						break;
					}
				}
			}
			
		
		int iter0=iterations;
		double tol0=tolerance;
		int iloc=1; //relaksacja lokalna
		conduct_relaxation(ipos,jpos,irange_min,irange_max,&iter0,&tol0,crystal,tol_local_max,iloc,as,ag,al,skl,skd);
			
			
		}//ktory==0: depozycja atomu
		else{  
			
		/****************************************************************************************************
		 * 				ktory>0: dyfuzja losowego atomu   
		 * 
		 ****************************************************************************************************/
		
			int ii_shift=gen_discrete_1_K_multiply_sign(k_max_step); //losowe przesuniecie lewo-prawo
			int i_old=lround(atoms_diff[ktory][1]); //aktualna pozycja atomu dyfundujacego
			int j_old=lround(atoms_diff[ktory][2]); //aktualna pozycja
			int i_new=(i_old+ii_shift+nx)%(nx);     //nowa pozycja atomu
			double typ;				
			int j_new;
			
			int imin,imax,jmin,jmax,i_nodes;
			int iter0;
			double tol0;
			int ierr=0;
			
			int irange=25;
			
			
			auto s1 = std::chrono::high_resolution_clock::now();
		//**** liczymy stara energie [-irange,irange]:  atom-on  *****************************
			double en_old=0.;
			for(int i=-irange;i<=irange; i++){
				for(int j=-irange;j<=irange; j++){
					int ii=(i_old+i+nx)%nx;
					int jj=j_old+j;
					if(jj>1 && jj<(ny-1) && crystal[ii][jj][0]>0.5){
						en_old+=compute_elastic_energy_wij(ii,jj,crystal,as,ag,al,skl,skd); //elastic energy
					}
				}
			}
			
			
			auto s2 = std::chrono::high_resolution_clock::now();
		//**** liczymy energie po usunieciu atomu (i_old,j_old):   atom-off  *************				
			crystal_copy=crystal; //kopia: atom-on
			crystal_copy[i_old][j_old][0]=0.; // usuwamy atom
			
		//***** relaksacja naprezen w sieci atom-off *************************************
			iter0=iterations;
			tol0=tolerance;
			int iloc=1; //relaksacja lokalna
			conduct_relaxation(i_old,j_old,irange_min,irange_max,&iter0,&tol0,crystal_copy,tol_local_max,iloc,as,ag,al,skl,skd);
			
			
			auto s3 = std::chrono::high_resolution_clock::now();
		//**** liczymy nowa energie - atom-off [-range,irange]*****************************
			double en_new=0.;
			for(int i=-irange;i<=irange; i++){
				for(int j=-irange;j<=irange; j++){
					int ii=(i_old+i+nx)%nx;
					int jj=j_old+j;
					if(jj>1 && jj<(ny-1) && crystal_copy[ii][jj][0]>0.5){
						en_new+=compute_elastic_energy_wij(ii,jj,crystal_copy,as,ag,al,skl,skd); //elastic energy
					}
				}
			}
			
			auto s4 = std::chrono::high_resolution_clock::now();
			auto s21 = std::chrono::duration_cast<std::chrono::microseconds>(s2-s1);
			auto s32 = std::chrono::duration_cast<std::chrono::microseconds>(s3-s2);
			auto s43 = std::chrono::duration_cast<std::chrono::microseconds>(s4-s3);
			//printf("%15.5E   %15.5E   %15.5E   %d\n",s21*1.0E-6,s32*1.0E-6,s43*1.0E-6,iter0);
	
			
			
			
			
		//******* sprawdzamy czy prawdopodobienstwo r_new<r_old=r_approx - jesli tak to atom dyfunduje   ***********************
				int neigh=lround(atoms_diff[ktory][7]); //aktualna liczba najblizszych sasiadow
				r0=12.*D/(k_max_step+1.0)/(2*k_max_step+2.);
				double ri_new=r0*exp((-gamma*neigh+(en_old-en_new)/2+E)/kbt); //dokladne tempo dyfuzji atomu
				double ri_old=atoms_diff[ktory][6]; //szacowane prawodpodobienstwo dyfuzji
			
			
			
				n_proposed++;
	
			if(gen_uniform()<(ri_new/ri_old)){ 
				n_accepted++;
				typ=crystal[i_old][j_old][0]; // zachowujemy typ atomu
				crystal[i_old][j_old][0]=0;   // kasujemy atom w starej pozycji		
				
				crystal=crystal_copy;
				//osadzamy atom w j_new
				for(int j=ny-1;j>=1;j--){
					int ll=lround(crystal[i_new][j][0]);
					if(ll>0){
						if(j<(ny-1)){
							j_new=j+1; //tu przesuwamy atom na wolne miejsce
							crystal[i_new][j_new][0]=typ; 
							break;
						}else{
							printf("nie mozemy przesunac atomu do nowej pozycji - brak miejsca w tabilicy \n\n");
							break;
						}
							
					}
				}
				//lokalna relaksacja polozen atomow w starym i w nowym polozeniu
				iter0=iterations;
				tol0=tolerance;
				int iloc=1; //relaksacja lokalna
				conduct_relaxation(i_new,j_new,irange_min,irange_max,&iter0,&tol0,crystal,tol_local_max,iloc,as,ag,al,skl,skd);
				
			}

		}// ktory>0: dyfuzja
		
		
	}//time<time_simulation
	
}//growth_simulation





/*************************************************************************************************************************
 *************************************************************************************************************************
 *************************************************************************************************************************
			procedury pomocnicze 
 
**************************************************************************************************************************
**************************************************************************************************************************
**************************************************************************************************************************/



/*************************************************************************************************************************
 *  relakasacja lokalna/globalna:
 *  iloc=0,1: 
 * 			0 - globalna
 * 			1 - lokalna
 * irange_min - minimalny promien otoczenia
 * irange_max - maksymalny promien otoczenia
 * tolerance - tolerancja bledu podczas rozwiazywania ukladu rownan
 * iterations - maksymalna liczba iteracji w metodzie CG
 * tol_local_max - maksymalna akceptowalna wartosc bledu lokalnego
 * 
 *************************************************************************************************************************/		
inline void conduct_relaxation(const int & ipos, const int & jpos, const int & irange_min, const int & irange_max, 
				int * iterations, double * tolerance,
				vector<vector<vector<double>>> & crystal, double tol_local_max, int iloc,
				double as,double ag, double al, double skl,double skd)
{
	
	int nx=crystal.size();    //liczba komorek w x
	int ny=crystal[0].size(); //liczba komorek w y
	int imin;
	int jmin;
	int jmax;
	int i_nodes;
	int iter0;
	double tol0;
	int ierr;
	double tol_local=0.;
	double mu;
	int irange=irange_min;
	double bmax=0;
	
	
	/*************************************************
	 * warunek na relaksacje lokalna
	 * 
	 *************************************************/
	
	
		if(iloc==1){
			//do{	
				
								
				/*************************************************************************************
				 * lokalna zmiana polozen atomow w otoczeniu irange_min - minimalizacja naprezen
				 * nie ma sensu zwiekszac rozmiaru otoczenia i liczenia bledu bmax wielokrotnie
				 * bo wyznaczenie Acsr trwa zawsze 2.5 ms
				 * 
				 *************************************************************************************/
				imin=(ipos-irange+nx)%nx;
				jmin=max(jpos-irange,1);
				jmax=min(jpos+irange,ny-2);
				i_nodes=2*irange; //roznica: imin->imax
				iter0=*iterations;
				tol0=*tolerance;
				ierr=0;
				solve_linear_system_u_v(crystal,imin,i_nodes,jmin,jmax,as,ag,al,skl,skd,&iter0,&tol0,ierr,&bmax);
				
				//liczymy blad lokalny dla zasiegu: irange_max
				/*
					imin=(ipos-irange_max+nx)%nx;
					jmin=max(jpos-irange_max,1);
					jmax=min(jpos+irange_max,ny-2);
					i_nodes=2*irange_max; 
					ierr=1; //liczymy: max|R|=max|Au-f| -> wynik zwracany w tol0
				
					solve_linear_system_u_v(crystal,imin,i_nodes,jmin,jmax,as,ag,al,skl,skd,&iter0,&tol0,ierr,&bmax);
					irange+=10; // np.: 10,20,30,40,50 - zwiekszamy rozmiar otoczenia atomu
				
				*/
				//wartosc bledu lokalnego 
				mu=(ag-as)/as;
				tol_local=bmax/mu/as/skl;
								
				
			//}while(tol_local>tol_local_max && irange_min <= irange_max);
		}	
		
		
				
	/******************************************************
	 * warunek na wykonanie relaksacji globalnej:
	 * 
	 ******************************************************/
	
		if( tol_local>tol_local_max || iloc!=1 ){	
			imin=0;
			jmin=1;
			jmax=ny-2;
			i_nodes=nx-1;
			iter0=*iterations;
			tol0=*tolerance;
			ierr=0;
			solve_linear_system_u_v(crystal,imin,i_nodes,jmin,jmax,as,ag,al,skl,skd,&iter0,&tol0,ierr,&bmax);
		}
	
	
		*iterations=iter0;
		*tolerance=tol0;
	
	return;
		
}//conduct_relaxation






/**************************************************************************************************************************
 * liczymy wektor gradientu dla wszystkich atomow poza dolnym brzegiem
 * 
 *  zwracamy norme wektora gradientu oraz gradient:
 * 
 *  crystal[][][6]=dW/du_ij
 *  crystal[][][7]=dW/dv_ij
 * 
 **************************************************************************************************************************/
double compute_gradient(int nx, int ny,  vector<vector<vector<double>>> & crystal, 
				   const double & as,const double & ag,const double & al,const double & skl,const double & skd)
{
	
	double gradient_norm=0;
	
	for(int i=0;i<nx;i++){
		for(int j=1;j<(ny-1);j++){
			
			if(crystal[i][j][0]>0.5){
			
				double u0=crystal[i][j][1];
				double v0=crystal[i][j][2];
				double delta=0.01; //krok do liczenia pochodnych
				
				double epx,emx,epy,emy;
				
				crystal[i][j][1]=u0+delta;
				epx=compute_elastic_energy_wij(i,j,crystal,as,ag,al,skl,skd);
				
				crystal[i][j][1]=u0-delta;
				emx=compute_elastic_energy_wij(i,j,crystal,as,ag,al,skl,skd);
				
				crystal[i][j][1]=u0;
				
				
				crystal[i][j][2]=v0+delta;
				epy=compute_elastic_energy_wij(i,j,crystal,as,ag,al,skl,skd);
				
				crystal[i][j][2]=v0-delta;
				emy=compute_elastic_energy_wij(i,j,crystal,as,ag,al,skl,skd);
				
				crystal[i][j][2]=v0;
				
				
				crystal[i][j][6]=(epx-emx)/2/delta; // pochodna w x
				crystal[i][j][7]=(epy-emy)/2/delta; // pochodna w y
				
				gradient_norm+=pow(crystal[i][j][6],2)+pow(crystal[i][j][7],2);
			}else{
				crystal[i][j][6]=0.;
				crystal[i][j][7]=0.;
			}
		}
	}	

	gradient_norm=sqrt(gradient_norm);

	return gradient_norm;
}



/**************************************************************************************************************************
 * liczymy energie sprezystosci dla atomu (i,j)
 * 
 **************************************************************************************************************************/

double compute_elastic_energy_wij(const int & i, const int & j, const vector<vector<vector<double>>> & crystal, 
				   const double & as,const double & ag,const double & al,const double & skl,const double & skd){
					   
	
	if(crystal[i][j][0]<0.5){
		return 0.0;
	}
	
	int nx=crystal.size();    //liczba komorek w x
	int ny=crystal[0].size(); //liczba komorek w y

	static int ip[3][3];
	static double d1[3][3];
	static double d2[3][3];
	static double u[3][3];
	static double v[3][3];
	
	//----------wypelniamy lokalne macierze pomocnicze------------------------------------------
	for(int ii=0;ii<3;ii++){
		for(int jj=0;jj<3;jj++){
			int i3=(i+ii-1+nx)%(nx);
			int j3=j+jj-1;
			if(lround(crystal[i3][j3][0])>0)ip[ii][jj]=1; //jest atom					
			else ip[ii][jj]=0; //brak atomu
			
			u[ii][jj]=crystal[i3][j3][1];
			v[ii][jj]=crystal[i3][j3][2];
			
			int id=lround(crystal[i][j][0]*crystal[i3][j3][0]); // typ oddzialywania
			if(id==1){ //s-s
				d1[ii][jj]=0;
				d2[ii][jj]=0;			
			}
			else if(id==2 || id==4){ //g-g (4) lub g-s (2)
				d1[ii][jj]=ag-as;
				d2[ii][jj]=ag-al;			
			}		
		}
	}
	
		
	int ii,jj;
	double wxx,wyy,wxy,ep;
	
	ii=1; //atom centralny
	jj=1;
	
	wxx=
	 skl/2.*ip[ii][jj]*ip[ii+1][jj]*pow(u[ii+1][jj]-u[ii][jj]-d1[ii+1][jj],2)
	+skl/2.*ip[ii][jj]*ip[ii-1][jj]*pow(u[ii-1][jj]-u[ii][jj]+d1[ii-1][jj],2)
	+skd/4.*ip[ii][jj]*ip[ii+1][jj+1]*pow(u[ii+1][jj+1]-u[ii][jj]-d1[ii+1][jj+1],2)
	+skd/4.*ip[ii][jj]*ip[ii-1][jj-1]*pow(u[ii-1][jj-1]-u[ii][jj]+d1[ii-1][jj-1],2)
	+skd/4.*ip[ii][jj]*ip[ii+1][jj-1]*pow(u[ii+1][jj-1]-u[ii][jj]-d1[ii+1][jj-1],2)
	+skd/4.*ip[ii][jj]*ip[ii-1][jj+1]*pow(u[ii-1][jj+1]-u[ii][jj]+d1[ii-1][jj+1],2);
	
	
	wyy=
	 skl/2.*ip[ii][jj]*ip[ii][jj+1]*pow(v[ii][jj+1]-v[ii][jj]-d2[ii][jj+1],2)
	+skl/2.*ip[ii][jj]*ip[ii][jj-1]*pow(v[ii][jj-1]-v[ii][jj]+d2[ii][jj-1],2)
	+skd/4.*ip[ii][jj]*ip[ii+1][jj+1]*pow(v[ii+1][jj+1]-v[ii][jj]-d2[ii+1][jj+1],2)
	+skd/4.*ip[ii][jj]*ip[ii-1][jj-1]*pow(v[ii-1][jj-1]-v[ii][jj]+d2[ii-1][jj-1],2)
	+skd/4.*ip[ii][jj]*ip[ii+1][jj-1]*pow(v[ii+1][jj-1]-v[ii][jj]+d2[ii+1][jj-1],2)
	+skd/4.*ip[ii][jj]*ip[ii-1][jj+1]*pow(v[ii-1][jj+1]-v[ii][jj]-d2[ii-1][jj+1],2);
	
	wxy=
	 skd/4.*ip[ii][jj]*ip[ii-1][jj-1]*(u[ii-1][jj-1]-u[ii][jj]+d1[ii-1][jj-1])*(v[ii-1][jj-1]-v[ii][jj]+d2[ii-1][jj-1])
	+skd/4.*ip[ii][jj]*ip[ii+1][jj+1]*(u[ii+1][jj+1]-u[ii][jj]-d1[ii+1][jj+1])*(v[ii+1][jj+1]-v[ii][jj]-d2[ii+1][jj+1])
	-skd/4.*ip[ii][jj]*ip[ii+1][jj-1]*(u[ii+1][jj-1]-u[ii][jj]-d1[ii+1][jj-1])*(v[ii+1][jj-1]-v[ii][jj]+d2[ii+1][jj-1])
	-skd/4.*ip[ii][jj]*ip[ii-1][jj+1]*(u[ii-1][jj+1]-u[ii][jj]+d1[ii-1][jj+1])*(v[ii-1][jj+1]-v[ii][jj]-d2[ii-1][jj+1]);
	
	ep=wxx+wyy+2*wxy;
	
	return ep;
	
}





/***************************************************************************
 *   generatory liczb pseudolosowych
 * 
 *    UWAGA: generator rand() z biblioteki kompilatora C jest restartowany 
 *    gdy wywolywana jest procedura PARDISO do rozwiazania ukladu 
 *    (teraz uzywamy metody iteracyjnej CG)
 * 
 ***************************************************************************/

double gen_uniform(){
	return ((double)rand()/RAND_MAX);
}


int gen_discrete_1_K(const int & K){
	int k1=lround(floor(gen_uniform()*K)+1);	
	return k1;
}


int gen_sign(){
	if(gen_uniform()<0.5) return -1;
	else return 1;
}

int gen_discrete_1_K_multiply_sign( const int & K){
	return gen_discrete_1_K(K)*gen_sign();
}


void test_gen_discrete_sign(){
		
		int k=8;
		int n=100000;
		vector<double> hist(2*k+1,0.);
		
		
		for(int i=0;i<n;i++){
			int m=gen_discrete_1_K_multiply_sign(k);
			m=m+k;
			hist[m]++;
		}
	
		FILE *fp;
		fp=fopen("hist_d_s.dat","w");
		
		for(int i=0;i<(2*k+1);i++){
			fprintf(fp,"%6d   %15.5E\n",i-k,hist[i]);
		}
		
		
		fclose(fp);
		return;
	}



/***********************************************************************************************
 * rozwiazujemy uklad rownan (IERR=0) LUB liczymy blad lokalny rozwiazania (IERR=1)
 * 
 * tablica: crystal[0:nx-1][0:ny-1][k]
 *          crystal[0:nx-1][0:ny-1][0]=0,1,2 - empty/substrate,deposited atom
 *          crystal[0:nx-1][0:ny-1][1]=uij - shift in x
 *          crystal[0:nx-1][0:ny-1][2]=vij - shift in y
 *          crystal[0:nx-1][0:ny-1][3]=l - global index for uij 
 *          crystal[0:nx-1][0:ny-1][4]=l+1 - global index for vij
 *          crystal[0:nx-1][0:ny-1][3]=-1 - boundary or outside of linear system 
 * 
 * zakres ukladu rownan obejmuje: [imin:imax,jmin:jmax]
 * 
 * imin         - poczatek zakresu atomow do relaksacji
 * imax=imin+i_nodes - koniec zakresu atomow do relaksacji
 * 
 * macierz ukladu rownan w formacie CSR: acsr,icsr,jcsr
 * rozwiazanie ukladu iteracyjne: Conjugate Gradients
 * 
 * 
 * iterations - maksymalna liczba iteracji CG
 * tolerance - dopuszczalna tolerancja bledu rozwiazania CG
 * ierr=0,1:  0-rozwiazujemy uklad rownan, 1-liczymy norme max aktualnego rozwiazania
 * 
 ************************************************************************************************/
void solve_linear_system_u_v(vector<vector<vector<double>>> & crystal,int imin, int i_nodes,int jmin,int jmax, 
				     double as, double ag, double al, double skl, double skd,int *iterations, double *tolerance, 
				     int ierr, double *bmax){

	auto s1 = std::chrono::high_resolution_clock::now();
	
	int nx=crystal.size();    //liczba komorek w x
	int ny=crystal[0].size(); //liczba komorek w y
	
	//sprawdzamy dolny i gorny zakres aby nie wyjsc za tablice
	if(jmin<1 || jmax>ny-2){
		jmin=max(jmin,1);
		jmax=min(jmax,ny-2);
	}
	
	int imax=imin+abs(i_nodes); //imax moze byc wieksze od (nx-1) - indeks jest renormalizowany

	/*
	 * usuwamy stare indeksy globalne, wpisujemy blokade (wb Dirichleta): -1
	 * jesli pozniej zmienimy numer (0,1,2,3,...) to bedzie warunek Neumanna
	 */
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
		      crystal[i][j][3]=-1;
			crystal[i][j][4]=-1;
		}
	}
	
	/*
	 * numeracja globalna w ukladzie rownan: nrow -liczba zmiennych/wierszy
	 * indeksacja wierszy: 0:(nrow-1) 
	 * 
	 */
	
	int nrow_max=(abs(i_nodes)+1)*(jmax-jmin+1)*2; //maksymalna liczba wierszy (dwa kierunki (x,y))
	int nrow=0; //faktyczna liczba wierszy - do ustalenia
	vector<vector<int>> indx; //tablica indeksow
	indx.resize(nrow_max,vector<int>(3,-1));
	
	nrow=0;
	for(int ii=imin;ii<=imax;ii++){
		for(int j=jmin;j<=jmax;j++){
			int i=(ii+nx)%(nx); //fizyczny numer komorki	 
			if(crystal[i][j][0]>0.5){ //tylko komorki zajete przez atomy
				crystal[i][j][3]=nrow; //blokada zniesiona
				indx[nrow][0]=i;
				indx[nrow][1]=j;
				indx[nrow][2]=3; //indeks przesuniecie w 'x'
				nrow++;
				
				crystal[i][j][4]=nrow; //blokada zniesiona
				indx[nrow][0]=i;
				indx[nrow][1]=j;
				indx[nrow][2]=4; //indeks przesuniecie w 'y'
				nrow++;
			}
		}
	}
	
	
	/*******************************************************************************************************
	 * tablice w postaci CSR -  wyznaczamy elementy w wierszu i wpisujemy posortowane do tablicy glownej 
	 * tablica dla wartosci w pojedynczym wierszu - po wypelnieniu sortujemy
	 *******************************************************************************************************/
	
	int ncol=9*2; //maksymalna liczba elementow w wierszu - liczba sasiadow * liczba kierunkow
	vector<double> acol;
	vector<int> jcol;
	
	acol.resize(ncol+10,0.);
	jcol.resize(ncol+10,0); //w zerowym indeksie zapisujemy liczbe elementow w wierszu
	
	/*
	 * tablice globalne do rozwiazywania ukladu rownan - allokacja jak w C 
	 * 
	 */
	
	int nmax=nrow*9*2; // maksymalna liczba niezerowych elementow w wierszu * liczba wierszy
	double * acsr;
	int * icsr;
	int * jcsr;
	double * ff;
	double * xx;
	double * bb;
	acsr=(double *)malloc(nmax*sizeof(double));
	jcsr=(int *)malloc(nmax*sizeof(int));
	icsr=(int *)malloc((nrow+1)*sizeof(int));
	ff=(double *)malloc(nrow*sizeof(double));
	xx=(double *)malloc(nrow*sizeof(double));
	bb=(double *)malloc(nrow*sizeof(double));
	icsr[nrow]=0; //aktualna liczba NNZ	w macierzy ukladu 
		
	//tworzymy tablice lokalnego otoczenia punktu 3x3
	/*
	 *   00  01  02    - numeracja wezlow w otoczeniu wezla (i,j) centralnego (11)
	 *   10 (11) 12
	 *   20  21  22 
	 * 
	 */
	
	// obsadzenie sasiadow - pij
	vector<vector<int>> ip;   
	ip.resize(3,vector<int>(3,0));
	
	// rodzaj brzegu: 0-Dirichlet, 1-Neumann 
	vector<vector<int>> iboundary;
	iboundary.resize(3,vector<int>(3,0));
	
	// tablica oddzialywania - d1
	vector<vector<double>> d1;   
	d1.resize(3,vector<double>(3,0));
	
	// tablica oddzialywania - d2
	vector<vector<double>> d2;   
	d2.resize(3,vector<double>(3,0));
	
	/*================================================================================================
	 * generujemy elementy macierzowe i wektor wyrazow wolnych 
	 *================================================================================================*/
	for(int k=0;k<nrow;k++){ //numer wiersza globalnego
		
			int i=indx[k][0]; //atom centralny dla wiersza
			int j=indx[k][1];

			
			for(int ii=0;ii<3;ii++){
				for(int jj=0;jj<3;jj++){
					ip[ii][jj]=0;
					d1[ii][jj]=0;
					d2[ii][jj]=0;
					iboundary[ii][jj]=0;
				}
			}
			
			
			//----------wypelniamy lokalne macierze pomocnicze------------------------------------------
			for(int ii=0;ii<3;ii++){
				for(int jj=0;jj<3;jj++){
					int i3=(i+ii-1+nx)%(nx);
					int j3=j+jj-1;
					if(lround(crystal[i3][j3][0])>0)ip[ii][jj]=1; //jest atom					
					else ip[ii][jj]=0; //brak atomu
								
					if(crystal[i3][j3][3]<0){ 
						iboundary[ii][jj]=0; // brzeg: Dirichlet (wyraz przenosimy do wyrazow wolnych)
					}else{
						iboundary[ii][jj]=1; // brzeg: Neumann (wyrazy zostawiamy w macierzy A)
					}
					
					int id=lround(crystal[i][j][0]*crystal[i3][j3][0]); // typ oddzialywania
					if(id==1){ //s-s
						d1[ii][jj]=0;
						d2[ii][jj]=0;
					}
					else if(id==2 || id==4){ //g-g (4) lub g-s (2)
						d1[ii][jj]=ag-as;
						d2[ii][jj]=ag-al;
					}
					
				}
			}
		
		/*================================================================================
		* *********** liczymy elementy: A, F ******************************
		* A: format CSR - macierz rzadka (acsr,icsr,jcsr)
		* F=ff[nrow] - wektor wyrazow wolnych
		* 
		*================================================================================*/
			jcol[0]=0; // 0-brak elementow: liczbe elementow trzymamy w elemencie  jcol[0]
			
			int number;
			number=indx[k][2];  //number:  3-uij, 4-vij
			
			ff[k]=0.; //zerujemy element wektora wyrazow wolnych - usuwamy smieci z poprzednich iteracji
			fill(acol.begin(), acol.end(), 0.0);
			fill(jcol.begin(), jcol.end(), 0.0);
			
			if(number==3){
				compute_u_v_from_wxx(number,k,i,j,nx,skl,skd,ip,iboundary,d1,crystal,acol,jcol,ff);
				compute_u_v_from_wxy(number,k,i,j,nx,skd,ip,iboundary,d2,crystal,acol,jcol,ff);
			}else if(number==4){
				compute_u_v_from_wxx(number,k,i,j,nx,skl,skd,ip,iboundary,d2,crystal,acol,jcol,ff);
				compute_u_v_from_wxy(number,k,i,j,nx,skd,ip,iboundary,d1,crystal,acol,jcol,ff); 
			}
			sort_and_add_matrix_elements(nrow,k,jcol,acol,acsr,icsr,jcsr);
	}//k=row index
	
	/******************************************************
	 * rozwiazujemy uklad rownan A*(uv)=ff
	 *         Conjugate Gradients
	 * 
	 ******************************************************/
		int itmax0=*iterations;
		double tol0=*tolerance;
		
		for(int i=0;i<nrow;i++)xx[i]=0.0; 
		
		// wektor startowy to poprzednie rozwiazanie
		for(int k=0;k<nrow;k++){ //numer wiersza globalnego
			int i=indx[k][0];
			int j=indx[k][1];
			int number=indx[k][2]; //3-uij, 4-vij
			xx[k]=crystal[i][j][number-2]; //number-2: 1-uij, 2-vij
		}
		
	auto s2 = std::chrono::high_resolution_clock::now();	
	/****************************************************************************************
	 * ierr=0,1:    
	 * 			0 - rozwiazujemy uklad rownan
	 * 			1 - liczymy blad lokalny jak w publikacji
	 * 
	 ****************************************************************************************/
		
		
		
		if(ierr==0){	//rozwiazujemy uklad rownan
			
			solve_linear_system_CG_standard(nrow,acsr,icsr,jcsr,ff,xx,iterations,tolerance);	
			
			if(*tolerance>=1.0E-3 || *iterations>=itmax0){
				printf("solution:  iterations,  tolerance  =   %6d   %15.5E  \n\n",*iterations,*tolerance);
			}
		
			//zachowujemy nowe polozenia/przesuniecia atomow
			for(int k=0;k<nrow;k++){ //numer wiersza globalnego
				int i=indx[k][0];
				int j=indx[k][1];
				int number=indx[k][2]; //3-uij, 4-vij
				crystal[i][j][number-2]=xx[k]; //number-2: 1-uij, 2-vij
			}
		}	
		
	
		//norma max z wektora reszt - liczymy zawsze: ierr-dowolne
			compute_sparse_Ax_y(nrow,acsr,icsr,jcsr,xx,bb);  // bb = Acsr*xx	
			*bmax=0.;
			for(int i=0;i<nrow;i++){
				bb[i]=bb[i]-ff[i];
				if(abs(bb[i])>*bmax)*bmax=abs(bb[i]);
			}
	
	auto s3 = std::chrono::high_resolution_clock::now();
	auto s21 = std::chrono::duration_cast<std::chrono::microseconds>(s2-s1);
	auto s32 = std::chrono::duration_cast<std::chrono::microseconds>(s3-s2);
	//printf("+++>   %15.5E   %15.5E    %d   %d   %15.5E\n",s21*1.0E-6,s32*1.0E-6,*iterations,itmax0,*tolerance);
	
	
	free(acsr);
	free(icsr);
	free(jcsr);
	free(ff);
	free(xx);
	free(bb);
	
}//solve Au=F:end


/***************************************************************************************************************
*	mnozenie macierz-wektor: Ax=y
* 	macierz rzadka
* 
***************************************************************************************************************/
inline void compute_sparse_Ax_y(const int & n, double * acsr, int * icsr, int  * jcsr, double * x, double * y){
		
		for(int i=0;i<n;i++){
			double sum=0.;
			int col;
			for(int j=icsr[i];j<=icsr[i+1]-1;j++){
				col=jcsr[j];
				sum+=acsr[j]*x[col];
			}
			y[i]=sum;
		}
		return;
	}

	
		
/***************************************************************************************************************
*	iloczyn skalarny dwoch wektorow
* 
***************************************************************************************************************/
inline double scalar_product_x_y(const int & n, double * x, double * y){
		double res=0.;
		for(int i=0;i<n;i++)res+=x[i]*y[i];
		return res;
	}	


	
	
/***************************************************************************************************************
****************************************************************************************************************
*	solve linear equations system: CG - standard algorithm - Saad
* 
**************************************************************************************************************** 
***************************************************************************************************************/
	void solve_linear_system_CG_standard(const int & n, double * acsr, int * icsr, int  * jcsr, 
									  double * b, double * x, int * itmax, double *tol)
	{
		
		
		// sprawdzamy wektor wyrazow wolnych: jesli jest pusty to zwracamy rozwiazanie trywialne vec{x}=0
		double b_2=scalar_product_x_y(n,b,b);
		if( b_2 < 1.0E-10 ){
			for(int i=0;i<n;i++) x[i]=0.;
			*itmax=0;
			*tol=0.;
			return;
		}
		
		//algorytm CG
		
		double * rj;
		double * rjp1;
		double * xj;
		double * xjp1;
		double * tmp1;
		double * Apj;
		double * pj;
		double * pjp1;
		
		rj=  (double *)malloc(n*sizeof(double));
		rjp1=(double *)malloc(n*sizeof(double));
		xj=  (double *)malloc(n*sizeof(double));
		xjp1=(double *)malloc(n*sizeof(double));
		tmp1=(double *)malloc(n*sizeof(double));
		Apj= (double *)malloc(n*sizeof(double));
		pj= (double *)malloc(n*sizeof(double));
		pjp1= (double *)malloc(n*sizeof(double));
		
		
		compute_sparse_Ax_y(n,acsr,icsr,jcsr,x,tmp1); //A*x0
		
		for(int i=0;i<n;i++){
			xj[i]=x[i]; //proponowane rozwiazanie 
			rj[i]=b[i]-tmp1[i]; //b-A*x0
			pj[i]=rj[i];
		}
		
		double Apj_2;
		double rj_2;
		double rjp1_2;
		double xj_2;
		double err;
		double alfa;
		double beta;
		int ile;
			

		
		for(int j=0;j<*itmax;j++){
			ile=j;	
			compute_sparse_Ax_y(n,acsr,icsr,jcsr,pj,Apj); //A*pj
			rj_2=scalar_product_x_y(n,rj,rj);
			Apj_2=scalar_product_x_y(n,Apj,pj);
			alfa=rj_2/Apj_2;
			if(fabs(alfa)<1.0E-5)printf("BLAD CG:  alfa= %15.5E \n",alfa);
			
			for(int i=0;i<n;i++)xjp1[i]=xj[i]+alfa*pj[i];
			for(int i=0;i<n;i++)rjp1[i]=rj[i]-alfa*Apj[i];
			rjp1_2=scalar_product_x_y(n,rjp1,rjp1);
			beta=rjp1_2/rj_2;
			for(int i=0;i<n;i++)pjp1[i]=rjp1[i]+beta*pj[i];
			
			for(int i=0;i<n;i++){
				xj[i]=xjp1[i];
				rj[i]=rjp1[i];
				pj[i]=pjp1[i];
			}
			
			xj_2=scalar_product_x_y(n,xj,xj);
			err=sqrt(rj_2)/sqrt(b_2);
			
			if(err< *tol && j>0){
				*itmax=j;
				break;
			}
		}
			*tol=err;
		
		
		//zapisujemy rozwiazanie
		for(int i=0;i<n;i++)x[i]=xj[i];
		
		
		//zwalniamy pamiec
		free(rj);
		free(rjp1);
		free(xj);
		free(xjp1);
		free(tmp1);
		free(Apj);
		free(pj);
		free(pjp1);
		
		return;
	}//CG-standard
	
	

/***************************************************************************************************************
 *************************************************************************************************************** 
 *  sortujemy wektor elementow macierzowych wzgledem kolumn (format CSR) i wkladamy do macierzy
 *  k - numer wiersza w macierzy ukladu   
 *  l=jcol[0] - liczba elementow niezerowych
 *  acol[1]-acol[l] - elementy niezerowe
 * 
 *  indeksowanie elementow od 0 - ostatni element w acsr lezy na pozycji acsr [nnz-1]
 * 
 ***************************************************************************************************************
 ***************************************************************************************************************/
inline void sort_and_add_matrix_elements(const int & nrow, const int & k, vector<int> & jcol, vector<double> & acol, 
					    double * acsr, int * icsr, int * jcsr){
	
	
	//sortowanie
	int l=jcol[0]; //liczba niezerowych elementow w wierszu
	int l1,l2;
	double a1,a2;
	
	for(int i=1;i<l;i++){
		for(int j=i;j>=1;j--){
			l1=jcol[j]; //numery kolumn
			l2=jcol[j+1];
			if(l1>l2){  //zamieniamy miejscami 
				a1=acol[j];
				a2=acol[j+1];
				acol[j]=a2;
				acol[j+1]=a1;
				jcol[j]=l2;
				jcol[j+1]=l1;
			}
		}	
	}
	
	if(l<1){
		printf("brak elementow w wierszu macierzy\n");
		exit(0);
	}
	
	//sprawdzamy numery kolumn
	for(int i=1;i<l;i++){
		if(jcol[i]>=jcol[i+1]){
			printf("blad w numerach kolumn:  %d   %d\n ",jcol[i],jcol[i+1]);
			exit(0);
		}
	}
	
	
	
	//dodajemy elementy do macierzy: acsr, jcsr, icsr
	int nnz;
	nnz=icsr[nrow]; //aktualna liczba elementow niezerowych - indeksowane od 0, 
	icsr[k]=nnz; //pozycja nnz jest pusta - od niej zaczynamy wypelnianie wiersza k-tego
	for(int i=1;i<=l;i++){
		acsr[nnz]=acol[i];
		jcsr[nnz]=jcol[i];
		nnz++;
	}
	icsr[nrow]=nnz; //zachowujemy aktualna wartosc nnz
	
	return;	
}//sort_and_add






/***************************************************************************************************************
 *************************************************************************************************************** 
 *    liczymy wkladu do wiersza dla wyrazu wxx/wyy - identycznie
 *    number=3,4:  
 * 			3-dwxx/duij, d=d1
 * 			4-dwyy/dvij, d=d2
 * 
 *    ii=1, jj=1: to punkt centralny
 * 
 *************************************************************************************************************** 
 ***************************************************************************************************************/
inline void compute_u_v_from_wxx(const int & number, const int & k,const int & i, const int & j, const int & nx, 
				  const double & skl, const double & skd, 
				  const vector<vector<int>> & ip,
				  const vector<vector<int>> & iboundary,
				  const vector<vector<double>> & d,
				  const vector<vector<vector<double>>> & crystal,
				  vector<double> & acol,	
				  vector<int> & jcol,
				  double * ff){
	
	int ii=1;
	int jj=1;
	int i3,lu;
	double val;
	
//element diagonalny 
	
	
	
	
	if(number==3){ //uij
	
		val=  -skl*ip[ii][jj]*ip[ii+1][jj]
			-skl*ip[ii][jj]*ip[ii-1][jj]
			-skd/2.*ip[ii][jj]*ip[ii+1][jj+1]
			-skd/2.*ip[ii][jj]*ip[ii-1][jj-1]
			-skd/2.*ip[ii][jj]*ip[ii+1][jj-1]
			-skd/2.*ip[ii][jj]*ip[ii-1][jj+1];	
		val=val*(-1);//pochodna wewnetrzna	
		
		lu=jcol[0]+1;	
		jcol[0]=lu;
		jcol[lu]=lround(crystal[i][j][number]);
		acol[lu]=val;
				
		// element wolny - wxx
		val=   skl*ip[ii][jj]*ip[ii+1][jj]*d[ii+1][jj]
			-skl*ip[ii][jj]*ip[ii-1][jj]*d[ii-1][jj]
			+skd/2.*ip[ii][jj]*ip[ii+1][jj+1]*d[ii+1][jj+1]
			-skd/2.*ip[ii][jj]*ip[ii-1][jj-1]*d[ii-1][jj-1]
			+skd/2.*ip[ii][jj]*ip[ii+1][jj-1]*d[ii+1][jj-1]
			-skd/2.*ip[ii][jj]*ip[ii-1][jj+1]*d[ii-1][jj+1];
		val=val*(-1);//pochodna wewnetrzna		
		ff[k]+=val;
		
		
	}else if(number==4){//vij
		
		val=  -skl*ip[ii][jj]*ip[ii][jj+1] 
			-skl*ip[ii][jj]*ip[ii][jj-1] 
			-skd/2.*ip[ii][jj]*ip[ii+1][jj+1]
			-skd/2.*ip[ii][jj]*ip[ii-1][jj-1]
			-skd/2.*ip[ii][jj]*ip[ii+1][jj-1]
			-skd/2.*ip[ii][jj]*ip[ii-1][jj+1];	
		val=val*(-1);//pochodna wewnetrzna				
		lu=jcol[0]+1;	
		jcol[0]=lu;
		jcol[lu]=lround(crystal[i][j][number]);
		acol[lu]=val;
		
		// element wolny - wyy
		val=   skl*ip[ii][jj]*ip[ii][jj+1]*d[ii][jj+1]
			-skl*ip[ii][jj]*ip[ii][jj-1]*d[ii][jj-1]
			+skd/2.*ip[ii][jj]*ip[ii+1][jj+1]*d[ii+1][jj+1]
			-skd/2.*ip[ii][jj]*ip[ii-1][jj-1]*d[ii-1][jj-1]
			-skd/2.*ip[ii][jj]*ip[ii+1][jj-1]*d[ii+1][jj-1]
			+skd/2.*ip[ii][jj]*ip[ii-1][jj+1]*d[ii-1][jj+1];
		val=val*(-1);//pochodna wewnetrzna				
		ff[k]+=val;
		
	}
	
//elementy pozadiagonalne: horyzontalne (number=3) i wertykalne (number=4) 
	if(number==3){
		
		for(int im=-1;im<=1;im+=2){
			int jm=0;
			i3=(i+im+nx)%(nx);
			val=skl*ip[ii][jj]*ip[ii+im][jj+jm];
			val=val*(-1);//pochodna wewnetrzna	
			if(iboundary[ii+im][jj+jm]==1){
				lu=jcol[0]+1;	
				jcol[0]=lu;
				jcol[lu]=lround(crystal[i3][j+jm][number]);
				acol[lu]=val;
			}
			else if(iboundary[ii+im][jj+jm]==0)
				ff[k]-=val*crystal[i3][j+jm][number-2];	
		}
		
	}else if(number==4){
		
		for(int jm=-1;jm<=1;jm+=2){
			int im=0;
			i3=(i+im+nx)%(nx);
			val=skl*ip[ii][jj]*ip[ii+im][jj+jm];
			val=val*(-1);//pochodna wewnetrzna	
			if(iboundary[ii+im][jj+jm]==1){
				lu=jcol[0]+1;	
				jcol[0]=lu;
				jcol[lu]=lround(crystal[i3][j+jm][number]);
				acol[lu]=val;
			}
			else if(iboundary[ii+im][jj+jm]==0)
				ff[k]-=val*crystal[i3][j+jm][number-2];	
		}
		
	}
	
//next nearest neighbours: pozostale diagonalne i antydiagonalne liczone identycznie	dla uij i vij
		for(int im=-1;im<=1;im+=2){
			for(int jm=-1;jm<=1;jm+=2){
				i3=(i+im+nx)%(nx);
				val=skd/2.*ip[ii][jj]*ip[ii+im][jj+jm];
				val=val*(-1);//pochodna wewnetrzna	
				if(iboundary[ii+im][jj+jm]==1){ //Neumann
					lu=jcol[0]+1;	
					jcol[0]=lu;	
					jcol[lu]=lround(crystal[i3][j+jm][number]);
					acol[lu]=val;
				}
				else if(iboundary[ii+im][jj+jm]==0) //Dirichlet
					ff[k]-=val*crystal[i3][j+jm][number-2];
				
			}
		}
	
		
	return;
}//compute_u_v_from_wxx



/***************************************************************************************************************
 *************************************************************************************************************** 
 *  liczymy wkladu do wiersza od wxy
 * 
 *  number=3,4:  
 * 			3-dW/duij
 * 			4-dW/dvij
 * 
 *  ii=1, jj=1: to punkt centralny
 *************************************************************************************************************** 
 ***************************************************************************************************************/
inline void compute_u_v_from_wxy(const int & number, const int & k, const int & i, const int & j, 
				  const int & nx, const double & skd, 
				  const vector<vector<int>> & ip,
				  const vector<vector<int>> & iboundary,
				  const vector<vector<double>> & d,
				  const vector<vector<vector<double>>> & crystal,
				  vector<double> & acol,	
				  vector<int> & jcol,
				  double * ff){
	int ii=1;
	int jj=1;
	double wsp=2.0; //mnoznik dla wxy w wij
	int i3,lu; 
	double val;
	
	int number_2;
	
	if(number==3){
		number_2=4; //indeks dla elementu v
	}else if(number==4){
		number_2=3; // indeks dla elementu u
	}
	
	for(int im=-1;im<=1;im+=2){
		for(int jm=-1;jm<=1;jm+=2){
			int sign=im*jm*(-1);
			i3=(i+im+nx)%(nx);
			double val=sign*skd/4.*ip[ii][jj]*ip[ii+im][jj+jm]*wsp;
			if(iboundary[ii+im][jj+jm]==1){
				lu=jcol[0]+1;	
				jcol[0]=lu;
				jcol[lu]=lround(crystal[i3][j+jm][ number_2 ]); //oddzialywanie: u->v, v->u
				acol[lu]=val;
			}
			else if(iboundary[ii+im][jj+jm]==0)
			ff[k]-=val*crystal[i3][j+jm][number_2-2];
			
		}
	}
	
	//element: vij*uij  - do diagonali w Acsr
		val=(
		skd/4.*ip[ii][jj]*ip[ii-1][jj-1]
		+skd/4.*ip[ii][jj]*ip[ii+1][jj+1]
		-skd/4.*ip[ii][jj]*ip[ii+1][jj-1]
		-skd/4.*ip[ii][jj]*ip[ii-1][jj+1])*wsp;
		
		lu=jcol[0]+1;	
		jcol[0]=lu;
		i3=(i+nx)%(nx);
		jcol[lu]=lround(crystal[i3][j][number_2]); //v
		acol[lu]=val;
		
	//element wolny: f(k)  - wxy
		
	if(number == 3){	
		val=(
		skd/4.*ip[ii][jj]*ip[ii-1][jj-1]*d[ii-1][jj-1]
		-skd/4.*ip[ii][jj]*ip[ii+1][jj+1]*d[ii+1][jj+1]
		-skd/4.*ip[ii][jj]*ip[ii+1][jj-1]*d[ii+1][jj-1]
		+skd/4.*ip[ii][jj]*ip[ii-1][jj+1]*d[ii-1][jj+1])*wsp;
		ff[k]=ff[k]+val;	
	}else if(number == 4){
		val=(
		skd/4.*ip[ii][jj]*ip[ii-1][jj-1]*d[ii-1][jj-1]
		-skd/4.*ip[ii][jj]*ip[ii+1][jj+1]*d[ii+1][jj+1]
		+skd/4.*ip[ii][jj]*ip[ii+1][jj-1]*d[ii+1][jj-1]
		-skd/4.*ip[ii][jj]*ip[ii-1][jj+1]*d[ii-1][jj+1])*wsp;
		ff[k]=ff[k]+val;	
	}
			
			return;
}//compute_u_v_from_wxy
