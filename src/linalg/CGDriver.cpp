#include "CGDriver.hpp"

#include <iostream>
#include <vector>

#include "Simulator.hpp"

CGDriver::CGDriver(Simulator* simulator,
                   const ConfigCG& conf_cg,
                   const ConfigPhysics& conf_ph)
    : MathDriver{simulator},
      max_iterations{conf_cg.max_iterations},
      tolerance{conf_cg.tolerance},
      m_substrate_lattice_constant{conf_ph.substrate_lattice_constant},
      m_adatom_lattice_constant{conf_ph.adatom_lattice_constant},
      m_vertical_lat_spacing{conf_ph.vertical_lat_spacing},
      m_spring_const_neighbors{conf_ph.spring_const_neighbors},
      m_spring_const_next_neighbors{conf_ph.spring_const_next_neighbors},
      m_misfit_coeff{conf_ph.misfit_coeff},
      m_D{conf_ph.D},
      m_E{conf_ph.E},
      m_bond_energy{conf_ph.bond_energy},
      last_iterations_no{},
      last_tolerance{} {}

void CGDriver::solve_linear_system(const std::size_t imin,
                                   const std::size_t i_nodes,
                                   std::size_t jmin,
                                   std::size_t jmax,
                                   int ierr,
                                   double* bmax,
                                   bool performOnCopy) {
  const std::size_t grid_x = sim->m_grid_x;
  const std::size_t grid_y = sim->m_grid_y;
  Grid& grid = performOnCopy ? sim->m_grid_copy : sim->m_grid;
  // check boundaries to prevent out-of-memory access
  if (jmin < 1 || jmax > grid_y - 2) {
    jmin = std::max(jmin, 1ul);
    jmax = std::min(jmax, grid_y - 2ul);
  }

  // NOTE: imax can be bigger than (grid_x-1) - index is renormalized
  const std::size_t imax = imin + i_nodes;

  // delete old global indexes, insert blockade (-1: Dirichlet boundary
  // condition), numbers (0,1,2,3,...) dictate Neumann boundary condition
  for (auto& row : grid) {
    for (auto& atom : row) {
      atom.boundary1 = -1;
      atom.boundary2 = -1;
    }
  }

  /*
   * numeracja globalna w ukladzie rownan: nrow -liczba zmiennych/wierszy
   * indeksacja wierszy: 0:(nrow-1)
   */

  // max no of rows (two directions (x,y))
  const std::size_t nrow_max = (i_nodes + 1) * (jmax - jmin + 1) * 2;

  // index table, initialized with -1
  std::vector<std::vector<int>> indx;
  indx.resize(nrow_max, std::vector<int>(3, -1));

  // calculate real no of rows
  std::size_t row_count = 0;
  for (std::size_t i = imin; i <= imax; i++) {
    for (std::size_t j = jmin; j <= jmax; j++) {
      const std::size_t i_wrapped =
          (i + grid_x) %
          (grid_x);  // TODO: fizyczny numer komorki, musimy to wrappować?
      // tylko komorki zajete przez atomy
      if (grid[i_wrapped][j].type != ATOM_TYPE::NO_ATOM) {
        // blokada zniesiona
        grid[i_wrapped][j].boundary1 = static_cast<int>(row_count);
        indx[row_count][0] = static_cast<int>(i_wrapped);
        indx[row_count][1] = static_cast<int>(j);
        // indeks przesuniecie w 'x'
        indx[row_count][2] = 3;
        row_count++;

        // blokada zniesiona
        grid[i_wrapped][j].boundary2 = static_cast<int>(row_count);
        indx[row_count][0] = static_cast<int>(i_wrapped);
        indx[row_count][1] = static_cast<int>(j);
        // indeks przesuniecie w 'y'
        indx[row_count][2] = 4;
        row_count++;
      }
    }
  }
  /*******************************************************************************************************
   * tablice w postaci CSR -  wyznaczamy elementy w wierszu i wpisujemy
   *posortowane do tablicy glownej tablica dla wartosci w pojedynczym wierszu -
   *po wypelnieniu sortujemy
   *******************************************************************************************************/
  // w zerowym indeksie acol i jcol zapisujemy liczbe elementow w wierszu
  std::array<double, column_count + 10> acol{};
  std::array<int, column_count + 10> jcol{};

  // tablice globalne do rozwiazywania ukladu rownan - allokacja jak w C

  // maksymalna liczba niezerowych elementow w wierszu * liczba wierszy
  const std::size_t nmax = row_count * 9 * 2;
  std::unique_ptr<double[]> csr_val{new double[nmax]};
  std::unique_ptr<int[]> icsr{new int[row_count + 1]};
  // aktualna liczba elementów niezerowych w macierzy ukladu
  icsr[row_count] = 0;
  std::unique_ptr<int[]> jcsr{new int[nmax]};

  std::unique_ptr<double[]> ff{new double[row_count]};
  std::unique_ptr<double[]> xx{new double[row_count]};
  std::unique_ptr<double[]> bb{new double[row_count]};

  /**
   * tworzymy tablice lokalnego otoczenia punktu 3x3
   *   00  01  02    - numeracja wezlow w otoczeniu wezla (i,j) centralnego (11)
   *   10 (11) 12
   *   20  21  22
   *
   */

  // obsadzenie sasiadow - pij
  Matrix3x3I ip;
  // rodzaj brzegu: 0-Dirichlet, 1-Neumann
  Matrix3x3I iboundary;
  // tablica oddzialywania d1
  Matrix3x3D d1_matrix;
  // tablica oddzialywania d2
  Matrix3x3D d2_matrix;

  /*================================================================================================
   * generujemy elementy macierzowe i wektor wyrazow wolnych
   *================================================================================================*/
  for (std::size_t k = 0; k < row_count; k++) {  // numer wiersza globalnego

    // experimental, wcześniej:
    // int i=indx[k][0]; //atom centralny dla wiersza
    // int j=indx[k][1];

    // atom centralny dla wiersza
    std::size_t i_central = static_cast<std::size_t>(indx[k][0]);
    std::size_t j_central = static_cast<std::size_t>(indx[k][1]);

    // fill helper matrices with 0s
    std::fill(ip.begin()->begin(), ip.back().end(), 0);
    std::fill(iboundary.begin()->begin(), iboundary.back().end(), 0);
    std::fill(d1_matrix.begin()->begin(), d1_matrix.back().end(), 0);
    std::fill(d2_matrix.begin()->begin(), d2_matrix.back().end(), 0);

    // wypelniamy lokalne macierze pomocnicze
    for (std::size_t i = 0; i < 3; i++) {
      for (std::size_t j = 0; j < 3; j++) {
        std::size_t i3 = (i_central + i - 1 + grid_x) % (grid_x);
        std::size_t j3 = j_central + j - 1;
        if (grid[i3][j3].type != ATOM_TYPE::NO_ATOM)
          ip[i][j] = 1;  // jest atom
        else
          ip[i][j] = 0;  // brak atomu

        if (grid[i3][j3].boundary1 < 0) {
          // brzeg: Dirichlet (wyraz przenosimy do wyrazow wolnych)
          iboundary[i][j] = 0;
        } else {
          // brzeg: Neumann (wyrazy zostawiamy w macierzy A)
          iboundary[i][j] = 1;
        }

        // d1 and d2 depend on bond type
        int id = static_cast<int>(grid[i_central][j_central].type) *
                 static_cast<int>(grid[i3][j3].type);
        if (id == 1) {
          // s-s
          d1_matrix[i][j] = 0;
          d2_matrix[i][j] = 0;
        } else if (id == 2 || id == 4) {
          // g-g (4) or g-s (2)
          d1_matrix[i][j] =
              m_adatom_lattice_constant - m_substrate_lattice_constant;
          d2_matrix[i][j] = m_adatom_lattice_constant - m_vertical_lat_spacing;
        }
      }
    }

    // caulculate A and F elements of system o eq
    // A: format CSR (- macierz rzadka (csr_val,icsr,jcsr)
    // F=ff[nrow] - wektor wyrazow wolnych

    // 0-brak elementow: liczbe elementow trzymamy w elemencie  jcol[0]
    jcol[0] = 0;

    // number:  3-uij, 4-vij
    std::size_t number = static_cast<std::size_t>(indx[k][2]);

    // zerujemy element wektora wyrazow wolnych - usuwamy smieci z poprzednich
    // iteracji
    ff[k] = 0.;
    std::fill(acol.begin(), acol.end(), 0.0);
    std::fill(jcol.begin(), jcol.end(), 0.0);

    if (number == 3) {
      compute_u_v_from_wxx(number, k, i_central, j_central, grid_x, ip,
                           iboundary, d1_matrix, grid, acol, jcol, ff.get());
      compute_u_v_from_wxy(number, k, i_central, j_central, grid_x, ip,
                           iboundary, d2_matrix, grid, acol, jcol, ff.get());
    } else if (number == 4) {
      compute_u_v_from_wxx(number, k, i_central, j_central, grid_x, ip,
                           iboundary, d2_matrix, grid, acol, jcol, ff.get());
      compute_u_v_from_wxy(number, k, i_central, j_central, grid_x, ip,
                           iboundary, d1_matrix, grid, acol, jcol, ff.get());
    }
    sort_and_add_matrix_elements(row_count, k, jcol, acol, csr_val.get(),
                                 icsr.get(), jcsr.get());
  }  // k=row index

  // rozwiazujemy uklad rownan A*(uv)=ff
  // Conjugate Gradients
  std::size_t itmax0 = max_iterations;

  for (std::size_t i = 0; i < row_count; i++)
    xx[i] = 0.0;

  // wektor startowy to poprzednie rozwiazanie
  for (std::size_t k = 0; k < row_count; k++) {  // numer wiersza globalnego
    std::size_t i = static_cast<std::size_t>(indx[k][0]);
    std::size_t j = static_cast<std::size_t>(indx[k][1]);
    std::size_t number = static_cast<std::size_t>(indx[k][2]);  // 3-uij, 4-vij
    // xx[k] = grid[i][j][number - 2];  // number-2: 1-uij, 2-vij
    if (number == 3)
      xx[k] = grid[i][j].u;
    else if (number == 4)
      xx[k] = grid[i][j].v;
  }

  // ierr=0,1:
  // 0 - rozwiazujemy uklad rownan
  // 1 - liczymy blad lokalny jak w publikacji

  if (ierr == 0) {
    // rozwiazujemy uklad rownan

    apply_CG_standard(row_count, csr_val.get(), icsr.get(), jcsr.get(),
                      ff.get(), xx.get());

    if (last_tolerance >= 1.0E-3 || last_iterations_no >= itmax0) {
      printf("solution:  iterations,  tolerance  =   %6ld   %15.5E  \n\n",
             last_iterations_no, last_tolerance);
    }

    // zachowujemy nowe polozenia/przesuniecia atomow
    for (std::size_t k = 0; k < row_count; k++) {  // numer wiersza globalnego
      std::size_t i = static_cast<std::size_t>(indx[k][0]);
      std::size_t j = static_cast<std::size_t>(indx[k][1]);
      std::size_t number =
          static_cast<std::size_t>(indx[k][2]);  // 3-uij, 4-vij
      if (number == 3)
        grid[i][j].u = xx[k];  // number-2: 1-uij, 2-vij
      else if (number == 4)
        grid[i][j].v = xx[k];  // number-2: 1-uij, 2-vij
    }
  }

  // norma max z wektora reszt - liczymy zawsze: ierr-dowolne
  compute_sparse_Ax_y(row_count, csr_val.get(), icsr.get(), jcsr.get(),
                      xx.get(),
                      bb.get());  // bb = csr_val*xx
  *bmax = 0.;
  for (std::size_t i = 0; i < row_count; i++) {
    bb[i] = bb[i] - ff[i];
    if (std::fabs(bb[i]) > *bmax)
      *bmax = std::fabs(bb[i]);
  }
}  // solve Au=F:end

void CGDriver::compute_sparse_Ax_y(const std::size_t& n_rows,
                                   double* csr_val,
                                   int* csr_row,
                                   int* csr_column,
                                   double* input_vector,
                                   double* output_vector) {
  // iterate over rows
  for (std::size_t i = 0; i < n_rows; i++) {
    double sum = 0;
    int col;
    for (int j = csr_row[i]; j <= csr_row[i + 1] - 1; j++) {
      col = csr_column[j];
      sum += csr_val[j] * input_vector[col];
    }
    output_vector[i] = sum;
  }
  return;
}

double CGDriver::scalar_product(const std::size_t& n, double* x, double* y) {
  double res = 0.;
  for (std::size_t i = 0; i < n; i++) {
    res += x[i] * y[i];
  }
  return res;
}

void CGDriver::apply_CG_standard(const std::size_t& row_count,
                                 double* csr_val,
                                 int* csr_row,
                                 int* csr_col,
                                 double* b,
                                 double* x) {
  // if b is zero vector, return trivial solution x={0,...}
  double f_2 = scalar_product(row_count, b, b);
  if (f_2 < 1.0E-10) {
    for (std::size_t i = 0; i < row_count; i++) {
      x[i] = 0.;
    }
    last_iterations_no = 0;
    last_tolerance = 0.;
    return;
  }

  // residual vectors
  std::unique_ptr<double[]> rj{new double[row_count]};
  std::unique_ptr<double[]> rj_proposed{new double[row_count]};
  // approximate solution vectors
  std::unique_ptr<double[]> xj{new double[row_count]};
  std::unique_ptr<double[]> xj_proposed{new double[row_count]};
  // vector for storing A*x0
  std::unique_ptr<double[]> A_times_x0{new double[row_count]};
  // matrix-vector product result (A*pj)
  std::unique_ptr<double[]> A_times_pj{new double[row_count]};
  // direction vectors
  std::unique_ptr<double[]> pj{new double[row_count]};
  std::unique_ptr<double[]> pj_proposed{new double[row_count]};

  // compute initial guess for A*x0, store in tmp1
  compute_sparse_Ax_y(row_count, csr_val, csr_row, csr_col, x,
                      A_times_x0.get());

  for (std::size_t i = 0; i < row_count; i++) {
    xj[i] = x[i];                  // initial guess
    rj[i] = b[i] - A_times_x0[i];  // b-A*x0
    pj[i] = rj[i];
  }

  double Apj_2;
  double rj_2;
  double rjp1_2;
  double approximation_error;
  double alfa;

  // iterate for maximum of itmax iterations
  for (std::size_t j = 0; j < max_iterations; j++) {
    // compute A*pj
    compute_sparse_Ax_y(row_count, csr_val, csr_row, csr_col, pj.get(),
                        A_times_pj.get());

    // compute step size alpha = r_j * r_j / (p_j * A * p_j)
    rj_2 = scalar_product(row_count, rj.get(), rj.get());
    Apj_2 = scalar_product(row_count, A_times_pj.get(), pj.get());
    alfa = rj_2 / Apj_2;
    if (std::fabs(alfa) < 1.0E-5)
      // step too small
      std::cerr << "Conjugate Gradient method error, step too small:  alfa="
                << alfa << std::endl;

    // update approximate solution vector
    for (std::size_t i = 0; i < row_count; i++) {
      xj_proposed[i] = xj[i] + alfa * pj[i];
    }

    // update residual vector
    for (std::size_t i = 0; i < row_count; i++) {
      rj_proposed[i] = rj[i] - alfa * A_times_pj[i];
    }

    // compute the update factor beta = r_{j+1} * r_{j+1} / (r_j * r_j)
    rjp1_2 = scalar_product(row_count, rj_proposed.get(), rj_proposed.get());
    double beta = rjp1_2 / rj_2;
    for (std::size_t i = 0; i < row_count; i++)
      pj_proposed[i] = rj_proposed[i] + beta * pj[i];

    // insert computed tmp solutions to corresponding vectors
    for (std::size_t i = 0; i < row_count; i++) {
      xj[i] = xj_proposed[i];
      rj[i] = rj_proposed[i];
      pj[i] = pj_proposed[i];
    }

    // compute solution error
    approximation_error = sqrt(rj_2) / sqrt(f_2);
    // if error is satisfying, end procedure
    if (approximation_error < tolerance && j > 0) {
      last_iterations_no = j;
      break;
    }
  }
  last_tolerance = approximation_error;

  // save computed solution
  for (std::size_t i = 0; i < row_count; i++)
    x[i] = xj[i];
}  // CG-standard

void CGDriver::compute_u_v_from_wxx(const std::size_t mode,
                                    const std::size_t k,
                                    const std::size_t i_central,
                                    const std::size_t j_central,
                                    const std::size_t grid_x,
                                    const Matrix3x3I& ip,
                                    const Matrix3x3I& iboundary,
                                    const Matrix3x3D& d_matrix,
                                    const Grid& grid,
                                    std::array<double, column_count + 10>& acol,
                                    std::array<int, column_count + 10>& jcol,
                                    double* ff) {
  std::size_t ii = 1;
  std::size_t jj = 1;
  std::size_t lu;
  std::size_t i3;
  double val;

  // u_{ij}
  if (mode == 3) {
    val = -m_spring_const_neighbors * ip[ii][jj] * ip[ii + 1][jj] -
          m_spring_const_neighbors * ip[ii][jj] * ip[ii - 1][jj] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii + 1][jj + 1] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii - 1][jj - 1] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii + 1][jj - 1] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii - 1][jj + 1];
    val = val * (-1);  // pochodna wewnetrzna

    lu = static_cast<std::size_t>(jcol[0] + 1);
    jcol[0] = static_cast<int>(lu);
    jcol[lu] =
        static_cast<int>(std::lround(grid[i_central][j_central].boundary1));
    acol[lu] = val;

    // element wolny - wxx
    val = m_spring_const_neighbors * ip[ii][jj] * ip[ii + 1][jj] *
              d_matrix[ii + 1][jj] -
          m_spring_const_neighbors * ip[ii][jj] * ip[ii - 1][jj] *
              d_matrix[ii - 1][jj] +
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii + 1][jj + 1] *
              d_matrix[ii + 1][jj + 1] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii - 1][jj - 1] *
              d_matrix[ii - 1][jj - 1] +
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii + 1][jj - 1] *
              d_matrix[ii + 1][jj - 1] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii - 1][jj + 1] *
              d_matrix[ii - 1][jj + 1];
    val = val * (-1);  // pochodna wewnetrzna
    ff[k] += val;

  }

  // v_{ij}
  else if (mode == 4) {
    val = -m_spring_const_neighbors * ip[ii][jj] * ip[ii][jj + 1] -
          m_spring_const_neighbors * ip[ii][jj] * ip[ii][jj - 1] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii + 1][jj + 1] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii - 1][jj - 1] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii + 1][jj - 1] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii - 1][jj + 1];
    val = val * (-1);  // pochodna wewnetrzna
    lu = static_cast<std::size_t>(jcol[0] + 1);
    jcol[0] = static_cast<int>(lu);
    jcol[lu] = static_cast<int>(lround(grid[i_central][j_central].boundary2));
    acol[lu] = val;

    // element wolny - wyy
    val = m_spring_const_neighbors * ip[ii][jj] * ip[ii][jj + 1] *
              d_matrix[ii][jj + 1] -
          m_spring_const_neighbors * ip[ii][jj] * ip[ii][jj - 1] *
              d_matrix[ii][jj - 1] +
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii + 1][jj + 1] *
              d_matrix[ii + 1][jj + 1] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii - 1][jj - 1] *
              d_matrix[ii - 1][jj - 1] -
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii + 1][jj - 1] *
              d_matrix[ii + 1][jj - 1] +
          m_spring_const_next_neighbors / 2. * ip[ii][jj] * ip[ii - 1][jj + 1] *
              d_matrix[ii - 1][jj + 1];
    val = val * (-1);  // pochodna wewnetrzna
    ff[k] += val;
  }

  // horizontal elements
  if (mode == 3) {
    for (int im = -1; im <= 1; im += 2) {
      std::size_t jm = 0;
      std::size_t i_offset =
          static_cast<std::size_t>(static_cast<int>(ii) + im);
      std::size_t j_offset = jj + jm;
      i3 = static_cast<std::size_t>(
          (static_cast<int>(i_central) + im + static_cast<int>(grid_x)) %
          static_cast<int>(grid_x));
      val = m_spring_const_neighbors * ip[ii][jj] * ip[i_offset][j_offset];
      val = val * (-1);  // pochodna wewnetrzna
      if (iboundary[i_offset][j_offset] == 1) {
        lu = static_cast<std::size_t>(jcol[0]) + 1;
        jcol[0] = static_cast<int>(lu);
        // jcol[lu] = static_cast<int>(lround(grid[i3][j_central + jm][mode]));
        jcol[lu] = static_cast<int>(lround(grid[i3][j_central + jm].boundary1));
        acol[lu] = val;
      } else if (iboundary[i_offset][j_offset] == 0)
        ff[k] -= val * grid[i3][j_central + jm].u;
    }

  }

  // vertical elements
  else if (mode == 4) {
    for (int jm = -1; jm <= 1; jm += 2) {
      std::size_t im = 0;
      std::size_t i_offset = ii + im;
      std::size_t j_offset =
          static_cast<std::size_t>(static_cast<int>(jj) + jm);
      i3 = (i_central + im + grid_x) % grid_x;
      val = m_spring_const_neighbors * ip[ii][jj] * ip[i_offset][j_offset];
      val = val * (-1);  // pochodna wewnetrzna
      if (iboundary[i_offset][j_offset] == 1) {
        lu = static_cast<std::size_t>(jcol[0] + 1);
        jcol[0] = static_cast<int>(lu);
        jcol[lu] = static_cast<int>(lround(
            grid[i3][static_cast<std::size_t>(static_cast<int>(j_central) + jm)]
                .boundary2));
        acol[lu] = val;
      } else if (iboundary[i_offset][j_offset] == 0)
        ff[k] -=
            val *
            grid[i3][static_cast<std::size_t>(static_cast<int>(j_central) + jm)]
                .v;
    }
  }

  // next nearest neighbours: pozostale diagonalne i antydiagonalne liczone
  // identycznie	dla uij i vij
  for (int im = -1; im <= 1; im += 2) {
    for (int jm = -1; jm <= 1; jm += 2) {
      std::size_t i_offset =
          static_cast<std::size_t>(static_cast<int>(ii) + im);
      std::size_t j_offset =
          static_cast<std::size_t>(static_cast<int>(jj) + jm);
      i3 = static_cast<std::size_t>(
          (static_cast<int>(i_central) + im + static_cast<int>(grid_x)) %
          static_cast<int>(grid_x));
      val = m_spring_const_next_neighbors / 2. * ip[ii][jj] *
            ip[i_offset][j_offset];
      val = val * (-1);                          // pochodna wewnetrzna
      if (iboundary[i_offset][j_offset] == 1) {  // Neumann
        lu = static_cast<std::size_t>(jcol[0] + 1);
        jcol[0] = static_cast<int>(lu);
        if (mode == 3)
          jcol[lu] = static_cast<int>(lround(
              grid[i3]
                  [static_cast<std::size_t>(static_cast<int>(j_central) + jm)]
                      .boundary1));
        else if (mode == 4)
          jcol[lu] = static_cast<int>(lround(
              grid[i3]
                  [static_cast<std::size_t>(static_cast<int>(j_central) + jm)]
                      .boundary2));
        acol[lu] = val;
      } else if (iboundary[i_offset][j_offset] == 0)  // Dirichlet
        if (mode == 3)
          ff[k] -= val * grid[i3][static_cast<std::size_t>(
                                      static_cast<int>(j_central) + jm)]
                             .u;
        else if (mode == 4)
          ff[k] -= val * grid[i3][static_cast<std::size_t>(
                                      static_cast<int>(j_central) + jm)]
                             .v;
    }
  }

  return;
}  // compute_u_v_from_wxx

void CGDriver::compute_u_v_from_wxy(const std::size_t mode,
                                    const std::size_t k,
                                    const std::size_t i_central,
                                    const std::size_t j_central,
                                    const std::size_t grid_x,
                                    const Matrix3x3I& ip,
                                    const Matrix3x3I& iboundary,
                                    const Matrix3x3D& d,
                                    const Grid& crystal,
                                    std::array<double, column_count + 10>& acol,
                                    std::array<int, column_count + 10>& jcol,
                                    double* ff) {
  std::size_t ii = 1;
  std::size_t jj = 1;
  double wsp = 2.0;  // mnoznik dla wxy w wij
  std::size_t lu;
  std::size_t i3;
  double val;

  std::size_t number;

  if (mode == 3) {
    number = 4;  // indeks dla elementu v
  } else if (mode == 4) {
    number = 3;  // indeks dla elementu u
  }

  for (int im = -1; im <= 1; im += 2) {
    for (int jm = -1; jm <= 1; jm += 2) {
      int sign = im * jm * (-1);
      std::size_t i_offset =
          static_cast<std::size_t>(static_cast<int>(ii) + im);
      std::size_t j_offset =
          static_cast<std::size_t>(static_cast<int>(jj) + jm);
      i3 = static_cast<std::size_t>(
          (static_cast<int>(i_central) + im + static_cast<int>(grid_x)) %
          static_cast<int>(grid_x));
      double val_local = sign * m_spring_const_next_neighbors / 4. *
                         ip[ii][jj] * ip[i_offset][j_offset] * wsp;
      if (iboundary[i_offset][j_offset] == 1) {
        lu = static_cast<std::size_t>(jcol[0] + 1);
        jcol[0] = static_cast<int>(lu);
        if (number == 3)
          jcol[lu] = static_cast<int>(
              lround(crystal[i3][static_cast<std::size_t>(
                                     static_cast<int>(j_central) + jm)]
                         .boundary1));  // oddzialywanie:
                                        // u->v, v->u
        else if (number == 4)
          jcol[lu] = static_cast<int>(
              lround(crystal[i3][static_cast<std::size_t>(
                                     static_cast<int>(j_central) + jm)]
                         .boundary2));  // oddzialywanie:
                                        // u->v, v->u
        acol[lu] = val_local;
      } else if (iboundary[i_offset][j_offset] == 0)
        if (number == 3)
          ff[k] -=
              val_local * crystal[i3][static_cast<std::size_t>(
                                          static_cast<int>(j_central) + jm)]
                              .u;
        else if (number == 4)
          ff[k] -=
              val_local * crystal[i3][static_cast<std::size_t>(
                                          static_cast<int>(j_central) + jm)]
                              .v;
    }
  }

  // element: vij*uij  - do diagonali w csr_val
  val = (m_spring_const_next_neighbors / 4. * ip[ii][jj] * ip[ii - 1][jj - 1] +
         m_spring_const_next_neighbors / 4. * ip[ii][jj] * ip[ii + 1][jj + 1] -
         m_spring_const_next_neighbors / 4. * ip[ii][jj] * ip[ii + 1][jj - 1] -
         m_spring_const_next_neighbors / 4. * ip[ii][jj] * ip[ii - 1][jj + 1]) *
        wsp;

  lu = static_cast<std::size_t>(jcol[0] + 1);
  jcol[0] = static_cast<int>(lu);
  i3 = (i_central + grid_x) % (grid_x);
  if (number == 3)
    jcol[lu] = static_cast<int>(lround(crystal[i3][j_central].boundary1));  // v
  else if (number == 4)
    jcol[lu] = static_cast<int>(lround(crystal[i3][j_central].boundary2));  // v
  acol[lu] = val;

  // element wolny: f(k)  - wxy

  if (mode == 3) {
    val = (m_spring_const_next_neighbors / 4. * ip[ii][jj] *
               ip[ii - 1][jj - 1] * d[ii - 1][jj - 1] -
           m_spring_const_next_neighbors / 4. * ip[ii][jj] *
               ip[ii + 1][jj + 1] * d[ii + 1][jj + 1] -
           m_spring_const_next_neighbors / 4. * ip[ii][jj] *
               ip[ii + 1][jj - 1] * d[ii + 1][jj - 1] +
           m_spring_const_next_neighbors / 4. * ip[ii][jj] *
               ip[ii - 1][jj + 1] * d[ii - 1][jj + 1]) *
          wsp;
    ff[k] = ff[k] + val;
  } else if (mode == 4) {
    val = (m_spring_const_next_neighbors / 4. * ip[ii][jj] *
               ip[ii - 1][jj - 1] * d[ii - 1][jj - 1] -
           m_spring_const_next_neighbors / 4. * ip[ii][jj] *
               ip[ii + 1][jj + 1] * d[ii + 1][jj + 1] +
           m_spring_const_next_neighbors / 4. * ip[ii][jj] *
               ip[ii + 1][jj - 1] * d[ii + 1][jj - 1] -
           m_spring_const_next_neighbors / 4. * ip[ii][jj] *
               ip[ii - 1][jj + 1] * d[ii - 1][jj + 1]) *
          wsp;
    ff[k] = ff[k] + val;
  }

  return;
}  // compute_u_v_from_wxy

void CGDriver::sort_and_add_matrix_elements(
    const std::size_t& nrow,
    const std::size_t& k,
    std::array<int, column_count + 10>& jcol,
    std::array<double, column_count + 10>& acol,
    double* csr_val,
    int* icsr,
    int* jcsr) {
  // sorting, double bubble sort
  // number of non-zero elements in row
  std::size_t l = static_cast<std::size_t>(jcol[0]);
  for (std::size_t i = 1; i < l; i++) {
    for (std::size_t j = i; j >= 1; j--) {
      if (jcol[j] > jcol[j + 1]) {  // zamieniamy miejscami
        std::swap(acol[j], acol[j + 1]);
        std::swap(jcol[j], jcol[j + 1]);
      }
    }
  }

  if (l < 1) {
    std::cerr << "No elements in matrix row\n";
    std::exit(1);
  }

  // dodajemy elementy do macierzy: csr_val, jcsr, icsr

  // aktualna liczba elementow niezerowych - indeksowane od 0,
  int nnz = icsr[nrow];
  // pozycja nnz jest pusta - od niej zaczynamy wypelnianie
  // wiersza k-tego
  icsr[k] = nnz;
  for (std::size_t i = 1; i <= l; i++) {
    csr_val[nnz] = acol[i];
    jcsr[nnz] = jcol[i];
    nnz++;
  }
  icsr[nrow] = nnz;  // zachowujemy aktualna wartosc nnz

  return;
}  // sort_and_add
