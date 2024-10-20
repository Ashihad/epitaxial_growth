#include <vector>

void conduct_relaxation(const std::size_t& ipos,
                        const std::size_t& jpos,
                        const std::size_t& irange_min,
                        const std::size_t& irange_max,
                        std::size_t* iterations,
                        double* tolerance,
                        std::vector<std::vector<std::vector<double>>>& crystal,
                        double tol_local_max,
                        int iloc,
                        double as,
                        double ag,
                        double al,
                        double skl,
                        double skd);

double compute_elastic_energy_wij(
    const std::size_t& i,
    const std::size_t& j,
    const std::vector<std::vector<std::vector<double>>>& crystal,
    const double& as,
    const double& ag,
    const double& al,
    const double& skl,
    const double& skd);

double compute_gradient(std::size_t nx,
                        std::size_t ny,
                        std::vector<std::vector<std::vector<double>>>& crystal,
                        const double& as,
                        const double& ag,
                        const double& al,
                        const double& skl,
                        const double& skd);