#pragma once

#include <array>
#include <cstdlib>

#include "AtomContainers.hpp"
#include "CSRContainers.hpp"
#include "ConfigMathDriver.hpp"
#include "ConfigPhysics.hpp"
#include "MatrixTypes.hpp"

// maksymalna liczba elementow w wierszu - liczba sasiadow * liczba kierunkow
constexpr std::size_t column_count = 9 * 2;

class MathDriver {
 public:
  MathDriver(const ConfigMathDriver& conf_md, const ConfigPhysics& conf_ph);
  virtual ~MathDriver() = default;
  MathDriver(const MathDriver&) = delete;
  MathDriver(MathDriver&&) = delete;
  MathDriver& operator=(const MathDriver&) = delete;
  MathDriver&& operator=(MathDriver&&) = delete;

  void compute_displacements(Grid& grid,
                             const std::size_t imin,
                             const std::size_t i_nodes,
                             std::size_t jmin,
                             std::size_t jmax,
                             int ierr,
                             double* bmax);

 protected:
  const double m_substrate_lattice_constant;
  const double m_adatom_lattice_constant;
  const double m_vertical_lat_spacing;
  const double m_spring_const_neighbors;
  const double m_spring_const_next_neighbors;
  const double m_misfit_coeff;
  const double m_D;
  const double m_E;
  const double m_bond_energy;

  const long unsigned m_max_iterations;
  const double m_tolerance;

  // these params hold parameters related to last iterative solver run
  std::size_t m_last_iterations_no;
  double m_last_tolerance;

  // must be defined in child
  virtual void solve_linear_system(CSRMatrix& A,
                                   FastVector<double>& b,
                                   FastVector<double>& x) = 0;

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
  void matrix_times_vector(const CSRMatrix& A,
                           const FastVector<double>& x,
                           FastVector<double>& y);

  /**
   * Compute inner product of two vectors x and y, each of them of length n
   * @param x - 1st vector
   * @param y - 2nd vector
   * @return scalar product
   */
  double dot(const FastVector<double>& x, const FastVector<double>& y);

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
                            FastVector<double>& ff);

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
                            FastVector<double>& ff);

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
  void sort_and_add_matrix_elements(const std::size_t,
                                    std::array<int, column_count + 10>&,
                                    std::array<double, column_count + 10>&,
                                    CSRMatrix& A);
};
