#include <vector>

void solve_linear_system_u_v(
    std::vector<std::vector<std::vector<double>>>& crystal,
    std::size_t imin,
    std::size_t i_nodes,
    std::size_t jmin,
    std::size_t jmax,
    double as,
    double ag,
    double al,
    double skl,
    double skd,
    std::size_t* iterations,
    double* tolerance,
    int ierr,
    double* bmax);

void compute_sparse_Ax_y(const std::size_t& n,
                         double* acsr,
                         int* icsr,
                         int* jcsr,
                         double* x,
                         double* y);
double scalar_product_x_y(const std::size_t& n, double* x, double* y);

void solve_linear_system_CG_standard(const std::size_t& n,
                                     double* acsr,
                                     int* icsr,
                                     int* jcsr,
                                     double* b,
                                     double* x,
                                     std::size_t* itmax,
                                     double* tol);

void compute_u_v_from_wxx(
    const std::size_t& number,
    const std::size_t& k,
    const std::size_t& i,
    const std::size_t& j,
    const std::size_t& nx,
    const double& skl,
    const double& skd,
    const std::vector<std::vector<int>>& ip,
    const std::vector<std::vector<int>>& iboundary,
    const std::vector<std::vector<double>>& d,
    const std::vector<std::vector<std::vector<double>>>& crystal,
    std::vector<double>& acol,
    std::vector<int>& jcol,
    double* ff);

void compute_u_v_from_wxy(
    const std::size_t& number,
    const std::size_t& k,
    const std::size_t& i,
    const std::size_t& j,
    const std::size_t& nx,
    const double& skd,
    const std::vector<std::vector<int>>& ip,
    const std::vector<std::vector<int>>& iboundary,
    const std::vector<std::vector<double>>& d,
    const std::vector<std::vector<std::vector<double>>>& crystal,
    std::vector<double>& acol,
    std::vector<int>& jcol,
    double* ff);

void sort_and_add_matrix_elements(const std::size_t& nrow,
                                  const std::size_t& k,
                                  std::vector<int>& jcol,
                                  std::vector<double>& acol,
                                  double* acsr,
                                  int* icsr,
                                  int* jcsr);