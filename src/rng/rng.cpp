#include "rng.hpp"

#include <cstdlib>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/***************************************************************************
 *   generatory liczb pseudolosowych
 * 
 *    UWAGA: generator rand() z biblioteki kompilatora C jest restartowany 
 *    gdy wywolywana jest procedura PARDISO do rozwiazania ukladu 
 *    (teraz uzywamy metody iteracyjnej CG)
 * 
 ***************************************************************************/

double gen_uniform(){
	return (static_cast<double>(rand())/RAND_MAX);
}


std::size_t gen_discrete_1_K(const std::size_t & K){
	return static_cast<std::size_t>(lround(floor(gen_uniform()*static_cast<double>(K))+1));
}


int gen_sign(){
	if(gen_uniform()<0.5) return -1;
	else return 1;
}

std::size_t gen_discrete_1_K_multiply_sign( const std::size_t & K){
	return static_cast<std::size_t>(static_cast<int>(gen_discrete_1_K(K))*gen_sign());
}


// void test_gen_discrete_sign(){
		
// 		int k=8;
// 		int n=100000;
// 		vector<double> hist(2*k+1,0.);
		
		
// 		for(int i=0;i<n;i++){
// 			int m=gen_discrete_1_K_multiply_sign(k);
// 			m=m+k;
// 			hist[m]++;
// 		}
	
// 		FILE *fp;
// 		fp=fopen("hist_d_s.dat","w");
		
// 		for(int i=0;i<(2*k+1);i++){
// 			fprintf(fp,"%6d   %15.5E\n",i-k,hist[i]);
// 		}
		
		
// 		fclose(fp);
// 		return;
// 	}
	