#include <cstddef>

/**
 * procedura symulacji wzrostu krzystalu: 1+1
 */
void growth_simulation(double as,
                       double ag,
                       double al,
                       double skl,
                       double skd,
                       double mu,
                       double D,
                       double E,
                       double gamma,
                       std::size_t nx,
                       std::size_t ny,
                       std::size_t nsurf,
                       std::size_t irange_min,
                       std::size_t irange_max,
                       double tol_local_max,
                       double time_fluence,
                       double time_simulation,
                       std::size_t k_max_step,
                       double f_deposition,
                       double temperature,
                       int idiffusion,
                       long unsigned iterations,
                       const double tolerance,
                       int n_write_files,
                       int n_global_relaxation);