#include <vector>

void conduct_relaxation(const int & ipos, const int & jpos, const int & irange_min, const int & irange_max, 
				int * iterations, double * tolerance,
				std::vector<std::vector<std::vector<double>>> & crystal, double tol_local_max, int iloc,
				double as,double ag, double al, double skl, double skd);

double compute_elastic_energy_wij(const int & i, const int & j, const std::vector<std::vector<std::vector<double>>> & crystal, 
				   const double & as,const double & ag,const double & al,const double & skl,const double & skd);

double compute_gradient(int nx, int ny,  std::vector<std::vector<std::vector<double>>> & crystal, 
				   const double & as,const double & ag,const double & al,const double & skl,const double & skd);