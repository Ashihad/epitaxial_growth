#include <vector>

void solve_linear_system_u_v(std::vector<std::vector<std::vector<double>>>& crystal, int imin, int i_nodes,int jmin,int jmax,
				     double as, double ag, double al, double skl, double skd,int *iterations, double *tolerance, 
				     int ierr, double *bmax);

void compute_sparse_Ax_y(const int & n, double * acsr, int * icsr, int * jcsr, double * x, double * y);
double scalar_product_x_y(const int & n, double * x, double * y);

void solve_linear_system_CG_standard(const int & n, double * acsr, int * icsr, int  * jcsr, 
									  double * b, double * x, int * itmax, double *tol);

void compute_u_v_from_wxx(const int & number, const int & k,const int & i, const int & j, const int & nx, 
				  const double & skl, const double & skd, 
				  const std::vector<std::vector<int>> & ip,
				  const std::vector<std::vector<int>> & iboundary,
				  const std::vector<std::vector<double>> & d,
				  const std::vector<std::vector<std::vector<double>>> & crystal,
				  std::vector<double> & acol,	
				  std::vector<int> & jcol,
				  double * ff);

void compute_u_v_from_wxy(const int & number, const int & k, const int & i, const int & j, const int & nx, const double & skd,
				  const std::vector<std::vector<int>> & ip,
				  const std::vector<std::vector<int>> & iboundary,
				  const std::vector<std::vector<double>> & d,
				  const std::vector<std::vector<std::vector<double>>> & crystal,
				  std::vector<double> & acol,	
				  std::vector<int> & jcol,
				  double * ff);

void sort_and_add_matrix_elements(const int & nrow, const int & k, std::vector<int> & jcol, std::vector<double> & acol, 
					    double * acsr, int * icsr, int * jcsr);