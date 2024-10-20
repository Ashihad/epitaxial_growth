#include <vector>

/**
 * rozwiazujemy uklad rownan (IERR=0) LUB liczymy blad lokalny rozwiazania
 * (IERR=1)
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
 * ierr=0,1:  0-rozwiazujemy uklad rownan, 1-liczymy norme max aktualnego
 *rozwiazania
 *
 */
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

/**
 * Compute Ax=y where A is given in CSR format
 *
 * not comfortable with CSR?
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)
 *
 * @param n_rows - number of rows in matrix A
 * @param csr_val - CSR VAL array
 * @param csr_row - CSR ROW array
 * @param csr_column - CSR COLUMN array
 * @param input_vector - vector that we multiply matrix A by (x vector)
 * @param output_vector - result vector (y vector)
 * @return nothing, output_vector contains result
 */
void compute_sparse_Ax_y(const std::size_t& n,
                         double* acsr,
                         int* icsr,
                         int* jcsr,
                         double* x,
                         double* y);

/**
 * Compute inner product of two vectors x and y, each of them of length n
 * @param x - 1st vector
 * @param y - 2nd vector
 * @return scalar product
 */
double scalar_product_x_y(const std::size_t& n, double* x, double* y);

/**
 *	solve linear equations system: CG - standard algorithm - Saad
 *
 * @param n
 * @param acsr
 * @param icsr
 * @param jcsr
 * @param b
 * @param x
 * @param itmax
 * @param tol
 */
void solve_linear_system_CG_standard(const std::size_t& n,
                                     double* csr_val,
                                     int* csr_row,
                                     int* csr_col,
                                     double* b,
                                     double* x,
                                     std::size_t* itmax,
                                     double* tolerance);

/**
 * liczymy wkladu do wiersza dla wyrazu wxx/wyy - identycznie
 * number=3,4:
 * 3-dwxx/duij, d=d1
 * 4-dwyy/dvij, d=d2
 *
 * ii=1, jj=1: to punkt centralny
 */
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

/**
 *  liczymy wkladu do wiersza od wxy
 *
 *  number=3,4:
 * 			3-dW/duij
 * 			4-dW/dvij
 *
 *  ii=1, jj=1: to punkt centralny
 */
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

/**
 *
 * sortujemy wektor elementow macierzowych wzgledem kolumn (format CSR) i
 * wkladamy do macierzy k - numer wiersza w macierzy ukladu l=jcol[0] - liczba
 * elementow niezerowych acol[1]..acol[l] - elementy niezerowe
 *
 *  indeksowanie elementow od 0 - ostatni element w acsr lezy na pozycji acsr
 * [nnz-1]
 *
 *
 */
void sort_and_add_matrix_elements(const std::size_t& nrow,
                                  const std::size_t& k,
                                  std::vector<int>& jcol,
                                  std::vector<double>& acol,
                                  double* acsr,
                                  int* icsr,
                                  int* jcsr);