#pragma once

#include "AtomContainers.hpp"
#include "ConfigCG.hpp"
#include "ConfigPhysics.hpp"
#include "MathDriver.hpp"
#include "MatrixTypes.hpp"

class CGDriver : public MathDriver {
 public:
  CGDriver(Simulator*, const ConfigCG&, const ConfigPhysics&);
  virtual ~CGDriver() = default;
  CGDriver(const CGDriver&) = delete;
  CGDriver(CGDriver&&) = delete;
  CGDriver& operator=(const CGDriver&) = delete;
  CGDriver&& operator=(CGDriver&&) = delete;

  virtual void solve_linear_system(const std::size_t imin,
                                   const std::size_t i_nodes,
                                   std::size_t jmin,
                                   std::size_t jmax,
                                   int ierr,
                                   double* bmax,
                                   bool performOnCopy = false) override;

  const long unsigned max_iterations;
  const double tolerance;

  const double m_substrate_lattice_constant;
  const double m_adatom_lattice_constant;
  const double m_vertical_lat_spacing;
  const double m_spring_const_neighbors;
  const double m_spring_const_next_neighbors;
  const double m_misfit_coeff;
  const double m_D;
  const double m_E;
  const double m_bond_energy;

  // these params hold parameters related to last apply_CG_standard run
  std::size_t last_iterations_no;
  double last_tolerance;

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
  double scalar_product(const std::size_t& n, double* x, double* y);

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
  void apply_CG_standard(const std::size_t& row_count,
                         double* csr_val,
                         int* csr_row,
                         int* csr_col,
                         double* b,
                         double* x);

  /**
   * liczymy wkladu do wiersza dla wyrazu wxx/wyy - identycznie
   * number=3,4:
   * 3-dwxx/duij, d=d1
   * 4-dwyy/dvij, d=d2
   *
   * ii=1, jj=1: to punkt centralny
   */
  void compute_u_v_from_wxx(const std::size_t number,
                            const std::size_t k,
                            const std::size_t i_central,
                            const std::size_t j_central,
                            const std::size_t nx,
                            const Matrix3x3I& ip,
                            const Matrix3x3I& iboundary,
                            const Matrix3x3D& d,
                            const Grid& crystal,
                            std::array<double, column_count + 10>& acol,
                            std::array<int, column_count + 10>& jcol,
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
  void compute_u_v_from_wxy(const std::size_t number,
                            const std::size_t k,
                            const std::size_t i_central,
                            const std::size_t j_central,
                            const std::size_t nx,
                            const Matrix3x3I& ip,
                            const Matrix3x3I& iboundary,
                            const Matrix3x3D& d,
                            const Grid& crystal,
                            std::array<double, column_count + 10>& acol,
                            std::array<int, column_count + 10>& jcol,
                            double* ff);

  /**
   *
   * sortujemy wektor elementow macierzowych wzgledem kolumn (format CSR) i
   * wkladamy do macierzy k - numer wiersza w macierzy ukladu l=jcol[0] -
   * liczba elementow niezerowych acol[1]..acol[l] - elementy niezerowe
   *
   *  indeksowanie elementow od 0 - ostatni element w acsr lezy na pozycji
   * acsr [nnz-1]
   *
   *
   */
  void sort_and_add_matrix_elements(const std::size_t&,
                                    const std::size_t&,
                                    std::array<int, column_count + 10>&,
                                    std::array<double, column_count + 10>&,
                                    double*,
                                    int*,
                                    int*);
};