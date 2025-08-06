#include "MathDriver.hpp"

#include <cmath>
#include <iostream>

#include "AtomContainers.hpp"
#include "MatrixTypes.hpp"
#include "Simulator.hpp"

MathDriver::MathDriver(const ConfigMathDriver& conf_md,
                       const ConfigPhysics& conf_ph)
    : m_substrate_lattice_constant{conf_ph.substrate_lattice_constant},
      m_adatom_lattice_constant{conf_ph.adatom_lattice_constant},
      m_vertical_lat_spacing{conf_ph.vertical_lat_spacing},
      m_spring_const_neighbors{conf_ph.spring_const_neighbors},
      m_spring_const_next_neighbors{conf_ph.spring_const_next_neighbors},
      m_misfit_coeff{conf_ph.misfit_coeff},
      m_D{conf_ph.D},
      m_E{conf_ph.E},
      m_bond_energy{conf_ph.bond_energy},
      m_max_iterations{conf_md.max_iterations},
      m_tolerance{conf_md.tolerance},
      m_last_iterations_no{},
      m_last_tolerance{} {}

double MathDriver::compute_displacements(Grid& grid,
                                         const std::size_t imin,
                                         const std::size_t imax_old,
                                         std::size_t jmin,
                                         std::size_t jmax) {
  const std::size_t grid_x = grid.size();
  const std::size_t grid_y = grid[0].size();

  // check boundaries to prevent out-of-memory access
  if (jmin < 1 || jmax > grid_y - 2) {
    jmin = std::max(jmin, 1ul);
    jmax = std::min(jmax, grid_y - 2ul);
  }

  // NOTE: imax can be bigger than (grid_x-1) - index is renormalized
  const std::size_t imax = imin + imax_old;

  // delete old global indexes, insert blockade (-1: Dirichlet boundary
  // condition), numbers (0,1,2,3,...) dictate Neumann boundary condition
  for (auto& row : grid) {
    for (auto& atom : row) {
      atom.boundary1 = -1;
      atom.boundary2 = -1;
    }
  }

  // max no of rows (two directions (x,y))
  const std::size_t nrow_max = (imax_old + 1) * (jmax - jmin + 1) * 2;

  // index table, x_index, y_index, (move in x(3), move in y(4))
  std::vector<std::vector<std::size_t>> indx;
  indx.resize(nrow_max,
              std::vector<std::size_t>(3, std::numeric_limits<size_t>::max()));

  // calculate real no of rows
  std::size_t row_count = 0;
  for (std::size_t i = imin; i <= imax; i++) {
    for (std::size_t j = jmin; j <= jmax; j++) {
      const std::size_t i_wrapped = i % grid_x;
      // tylko komorki zajete przez atomy
      if (grid[i_wrapped][j].type != ATOM_TYPE::NO_ATOM) {
        // blokada zniesiona
        grid[i_wrapped][j].boundary1 = static_cast<int>(row_count);
        indx[row_count][0] = i_wrapped;
        indx[row_count][1] = j;
        // indeks przesuniecie w 'x'
        indx[row_count][2] = 3;
        row_count++;

        // blokada zniesiona
        grid[i_wrapped][j].boundary2 = static_cast<int>(row_count);
        indx[row_count][0] = i_wrapped;
        indx[row_count][1] = j;
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

  // tablice globalne do rozwiazywania ukladu rownan

  // maksymalna liczba niezerowych elementow w wierszu * liczba wierszy
  const std::size_t nmax = row_count * 9 * 2;
  CSRMatrix A(nmax, row_count);
  A.row_index[A.n_rows] = 0;

  FastVector<double> ff(row_count);
  FastVector<double> xx(row_count);
  FastVector<double> bb(row_count);

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

    // atom centralny dla wiersza
    std::size_t i_central = indx[k][0];
    std::size_t j_central = indx[k][1];

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
          // subst-subst
          d1_matrix[i][j] = 0;
          d2_matrix[i][j] = 0;
        } else if (id == 2 || id == 4) {
          // ad-ad (4) or ad-subst (2)
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
    std::size_t number = indx[k][2];

    // zerujemy element wektora wyrazow wolnych - usuwamy smieci z poprzednich
    // iteracji
    ff[k] = 0.;
    std::fill(acol.begin(), acol.end(), 0.);
    std::fill(jcol.begin(), jcol.end(), 0);

    if (number == 3) {
      compute_u_v_from_wxx(number, k, i_central, j_central, grid_x, ip,
                           iboundary, d1_matrix, grid, acol, jcol, ff);
      compute_u_v_from_wxy(number, k, i_central, j_central, grid_x, ip,
                           iboundary, d2_matrix, grid, acol, jcol, ff);
    } else if (number == 4) {
      compute_u_v_from_wxx(number, k, i_central, j_central, grid_x, ip,
                           iboundary, d2_matrix, grid, acol, jcol, ff);
      compute_u_v_from_wxy(number, k, i_central, j_central, grid_x, ip,
                           iboundary, d1_matrix, grid, acol, jcol, ff);
    }
    sort_and_add_matrix_elements(k, jcol, acol, A);
  }  // k=row index

  // rozwiazujemy uklad rownan A*(uv)=ff
  // Conjugate Gradients
  std::size_t itmax0 = m_max_iterations;

  std::fill_n(xx.get(), row_count, 0.0);

  // wektor startowy to poprzednie rozwiazanie
  for (std::size_t k = 0; k < row_count; k++) {  // numer wiersza globalnego
    std::size_t i = indx[k][0];
    std::size_t j = indx[k][1];
    std::size_t number = indx[k][2];  // 3-uij, 4-vij
    // xx[k] = grid[i][j][number - 2];  // number-2: 1-uij, 2-vij
    if (number == 3)
      xx[k] = grid[i][j].u;
    else if (number == 4)
      xx[k] = grid[i][j].v;
  }

  solve_linear_system(A, ff, xx);

  if (m_last_tolerance >= 1.0E-3 || m_last_iterations_no >= itmax0) {
    printf("solution:  iterations,  tolerance  =   %6ld   %15.5E  \n\n",
           m_last_iterations_no, m_last_tolerance);
  }

  // zachowujemy nowe polozenia/przesuniecia atomow
  for (std::size_t k = 0; k < row_count; k++) {  // numer wiersza globalnego
    std::size_t i = indx[k][0];
    std::size_t j = indx[k][1];
    std::size_t number = indx[k][2];  // 3-uij, 4-vij
    if (number == 3)
      grid[i][j].u = xx[k];  // number-2: 1-uij, 2-vij
    else if (number == 4)
      grid[i][j].v = xx[k];  // number-2: 1-uij, 2-vij
  }

  // norma max z wektora reszt - liczymy zawsze: ierr-dowolne
  // return value is important only for local relaxation
  matrix_times_vector(A, xx, bb);
  double biggest_abs_bi{};
  for (std::size_t i = 0; i < row_count; i++) {
    bb[i] -= ff[i];
    if (std::abs(bb[i]) > biggest_abs_bi)
      biggest_abs_bi = std::abs(bb[i]);
  }
  return biggest_abs_bi;
}

void MathDriver::matrix_times_vector(const CSRMatrix& A,
                                     const FastVector<double>& input,
                                     FastVector<double>& output) {
  // iterate over rows
  for (std::size_t i = 0; i < A.n_rows; i++) {
    double sum = 0;
    for (std::size_t j = A.row_index[i]; j <= A.row_index[i + 1] - 1; j++) {
      const std::size_t col = A.col_index[j];
      sum += A.value[j] * input[col];
    }
    output[i] = sum;
  }
  return;
}

double MathDriver::dot(const FastVector<double>& x,
                       const FastVector<double>& y) {
  double res = 0.;
  for (std::size_t i = 0; i < x.size; i++) {
    res += x[i] * y[i];
  }
  return res;
}

void MathDriver::compute_u_v_from_wxx(
    const std::size_t mode,
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
    FastVector<double>& ff) {
  // matrix central i index
  std::size_t ii = 1;
  // matrix central j index
  std::size_t jj = 1;
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

    int lu = jcol[0] + 1;
    jcol[0] = lu;
    jcol[static_cast<std::size_t>(lu)] = grid[i_central][j_central].boundary1;
    acol[static_cast<std::size_t>(lu)] = val;

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
    int lu = jcol[0] + 1;
    jcol[0] = lu;
    jcol[static_cast<std::size_t>(lu)] = grid[i_central][j_central].boundary2;
    acol[static_cast<std::size_t>(lu)] = val;

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
      std::size_t i3 = static_cast<std::size_t>(
          (static_cast<int>(i_central) + im + static_cast<int>(grid_x)) %
          static_cast<int>(grid_x));
      val = m_spring_const_neighbors * ip[ii][jj] * ip[i_offset][j_offset];
      val = val * (-1);  // pochodna wewnetrzna
      if (iboundary[i_offset][j_offset] == 1) {
        std::size_t lu = static_cast<std::size_t>(jcol[0]) + 1;
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
      std::size_t i3 = (i_central + im + grid_x) % grid_x;
      val = m_spring_const_neighbors * ip[ii][jj] * ip[i_offset][j_offset];
      val = val * (-1);  // pochodna wewnetrzna
      if (iboundary[i_offset][j_offset] == 1) {
        std::size_t lu = static_cast<std::size_t>(jcol[0] + 1);
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
      std::size_t i3 = static_cast<std::size_t>(
          (static_cast<int>(i_central) + im + static_cast<int>(grid_x)) %
          static_cast<int>(grid_x));
      val = m_spring_const_next_neighbors / 2. * ip[ii][jj] *
            ip[i_offset][j_offset];
      val = val * (-1);  // pochodna wewnetrzna
      // Neumann
      if (iboundary[i_offset][j_offset] == 1) {
        std::size_t lu = static_cast<std::size_t>(jcol[0] + 1);
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
      } else if (iboundary[i_offset][j_offset] == 0) {  // Dirichlet
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
  }

  return;
}  // compute_u_v_from_wxx

void MathDriver::compute_u_v_from_wxy(
    const std::size_t mode,
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
    FastVector<double>& ff) {
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
          jcol[lu] = crystal[i3][static_cast<std::size_t>(
                                     static_cast<int>(j_central) + jm)]
                         .boundary1;  // oddzialywanie:
                                      // u->v, v->u
        else if (number == 4)
          jcol[lu] = static_cast<int>(
              lround(crystal[i3][static_cast<std::size_t>(
                                     static_cast<int>(j_central) + jm)]
                         .boundary2));  // oddzialywanie:
                                        // u->v, v->u
        acol[lu] = val_local;
      } else if (iboundary[i_offset][j_offset] == 0) {
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

void MathDriver::sort_and_add_matrix_elements(
    const std::size_t k,
    std::array<int, column_count + 10>& jcol,
    std::array<double, column_count + 10>& acol,
    CSRMatrix& A) {
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

  // dodajemy elementy do macierzy A

  // aktualna liczba elementow niezerowych - indeksowane od 0,
  std::size_t nnz = A.row_index[A.n_rows];
  // pozycja nnz jest pusta - od niej zaczynamy wypelnianie
  // wiersza k-tego
  A.row_index[k] = nnz;
  for (std::size_t i = 1; i <= l; i++) {
    A.value[nnz] = acol[i];
    A.col_index[nnz] = static_cast<std::size_t>(jcol[i]);
    nnz++;
  }
  A.row_index[A.n_rows] = nnz;  // zachowujemy aktualna wartosc nnz

  return;
}  // sort_and_add
