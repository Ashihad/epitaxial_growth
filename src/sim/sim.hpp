#include <cstddef>

void growth_simulation(double as, double ag, double al, double skl, double skd, double mu, double D, double E, double gamma,
				int nx, int ny, int nsurf, int irange_min, int irange_max, double tol_local_max, 
			     double time_fluence, double time_simulation, int k_max_step,
			     double f_deposition, double temperature, int idiffusion, const int iterations, const double tolerance,
			     int n_write_files, int n_global_relaxation	);