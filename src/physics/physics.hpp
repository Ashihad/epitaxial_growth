#include <vector>

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

/**************************************************************************************************************************
 * liczymy energie sprezystosci dla atomu (i,j)
 *
 **************************************************************************************************************************/
double compute_elastic_energy_wij(
    const std::size_t& i,
    const std::size_t& j,
    const std::vector<std::vector<std::vector<double>>>& crystal,
    const double& as,
    const double& ag,
    const double& al,
    const double& skl,
    const double& skd);

/**************************************************************************************************************************
 * liczymy wektor gradientu dla wszystkich atomow poza dolnym brzegiem
 *
 *  zwracamy norme wektora gradientu oraz gradient:
 *
 *  crystal[][][6]=dW/du_ij
 *  crystal[][][7]=dW/dv_ij
 *
 **************************************************************************************************************************/
double compute_gradient(std::size_t nx,
                        std::size_t ny,
                        std::vector<std::vector<std::vector<double>>>& crystal,
                        const double& as,
                        const double& ag,
                        const double& al,
                        const double& skl,
                        const double& skd);