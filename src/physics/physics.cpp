#include "physics.hpp"

#include "linalg.hpp"
#include <cmath>

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
void conduct_relaxation(const int & ipos, const int & jpos, const int & irange_min, const int & irange_max, 
				int * iterations, double * tolerance,
				std::vector<std::vector<std::vector<double>>> & crystal, double tol_local_max, int iloc,
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
				jmin=std::max(jpos-irange,1);
				jmax=std::min(jpos+irange,ny-2);
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
 * liczymy energie sprezystosci dla atomu (i,j)
 * 
 **************************************************************************************************************************/

double compute_elastic_energy_wij(const int & i, const int & j, const std::vector<std::vector<std::vector<double>>> & crystal, 
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
			if(std::lround(crystal[i3][j3][0])>0)ip[ii][jj]=1; //jest atom					
			else ip[ii][jj]=0; //brak atomu
			
			u[ii][jj]=crystal[i3][j3][1];
			v[ii][jj]=crystal[i3][j3][2];
			
			int id=std::lround(crystal[i][j][0]*crystal[i3][j3][0]); // typ oddzialywania
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

/**************************************************************************************************************************
 * liczymy wektor gradientu dla wszystkich atomow poza dolnym brzegiem
 * 
 *  zwracamy norme wektora gradientu oraz gradient:
 * 
 *  crystal[][][6]=dW/du_ij
 *  crystal[][][7]=dW/dv_ij
 * 
 **************************************************************************************************************************/
double compute_gradient(int nx, int ny,  std::vector<std::vector<std::vector<double>>> & crystal, 
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
