#include "linalg.hpp"

#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdio>

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
void solve_linear_system_u_v(std::vector<std::vector<std::vector<double>>> & crystal,int imin, int i_nodes,int jmin,int jmax, 
				     double as, double ag, double al, double skl, double skd,int *iterations, double *tolerance, 
				     int ierr, double *bmax)
{

	auto s1 = std::chrono::high_resolution_clock::now();
	
	int nx=crystal.size();    //liczba komorek w x
	int ny=crystal[0].size(); //liczba komorek w y
	
	//sprawdzamy dolny i gorny zakres aby nie wyjsc za tablice
	if(jmin<1 || jmax>ny-2){
		jmin=std::max(jmin,1);
		jmax=std::min(jmax,ny-2);
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
	std::vector<std::vector<int>> indx; //tablica indeksow
	indx.resize(nrow_max,std::vector<int>(3,-1));
	
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
	std::vector<double> acol;
	std::vector<int> jcol;
	
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
	std::vector<std::vector<int>> ip;   
	ip.resize(3,std::vector<int>(3,0));
	
	// rodzaj brzegu: 0-Dirichlet, 1-Neumann 
	std::vector<std::vector<int>> iboundary;
	iboundary.resize(3,std::vector<int>(3,0));
	
	// tablica oddzialywania - d1
	std::vector<std::vector<double>> d1;   
	d1.resize(3,std::vector<double>(3,0));
	
	// tablica oddzialywania - d2
	std::vector<std::vector<double>> d2;   
	d2.resize(3,std::vector<double>(3,0));
	
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
				  const std::vector<std::vector<int>> & ip,
				  const std::vector<std::vector<int>> & iboundary,
				  const std::vector<std::vector<double>> & d,
				  const std::vector<std::vector<std::vector<double>>> & crystal,
				  std::vector<double> & acol,	
				  std::vector<int> & jcol,
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
				  const std::vector<std::vector<int>> & ip,
				  const std::vector<std::vector<int>> & iboundary,
				  const std::vector<std::vector<double>> & d,
				  const std::vector<std::vector<std::vector<double>>> & crystal,
				  std::vector<double> & acol,	
				  std::vector<int> & jcol,
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
inline void sort_and_add_matrix_elements(const int & nrow, const int & k, std::vector<int> & jcol, std::vector<double> & acol, 
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
