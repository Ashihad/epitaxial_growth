#include "linalg.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

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
    double* bmax) {
  auto start_time = std::chrono::high_resolution_clock::now();

  // grid x size
  std::size_t nx = crystal.size();
  // grid y size
  std::size_t ny = crystal[0].size();

  // check boundaries to prevent out-of-memory access
  if (jmin < 1 || jmax > ny - 2) {
    jmin = std::max(jmin, 1ul);
    jmax = std::min(jmax, ny - 2ul);
  }

  // NOTE: imax can be bigger than (nx-1) - index is renormalized
  std::size_t imax = imin + i_nodes;

  // delete old global indexes, insert blockade (-1: Dirichlet boundary
  // condition), numbers (0,1,2,3,...) dictate Neumann boundary condition
  for (std::size_t i = 0; i < nx; i++) {
    for (std::size_t j = 0; j < ny; j++) {
      crystal[i][j][3] = -1;
      crystal[i][j][4] = -1;
    }
  }

  /*
   * numeracja globalna w ukladzie rownan: nrow -liczba zmiennych/wierszy
   * indeksacja wierszy: 0:(nrow-1)
   *
   */

  std::size_t nrow_max = (i_nodes + 1) * (jmax - jmin + 1) *
                         2;  // maksymalna liczba wierszy (dwa kierunki (x,y))
  std::size_t nrow = 0;      // faktyczna liczba wierszy - do ustalenia
  std::vector<std::vector<int>> indx;  // tablica indeksow
  indx.resize(nrow_max, std::vector<int>(3, -1));

  for (std::size_t ii = imin; ii <= imax; ii++) {
    for (std::size_t j = jmin; j <= jmax; j++) {
      std::size_t i = (ii + nx) % (nx);  // fizyczny numer komorki
      if (crystal[i][j][0] > 0.5) {      // tylko komorki zajete przez atomy
        crystal[i][j][3] = static_cast<int>(nrow);  // blokada zniesiona
        indx[nrow][0] = static_cast<int>(i);
        indx[nrow][1] = static_cast<int>(j);
        indx[nrow][2] = 3;  // indeks przesuniecie w 'x'
        nrow++;

        crystal[i][j][4] = static_cast<int>(nrow);  // blokada zniesiona
        indx[nrow][0] = static_cast<int>(i);
        indx[nrow][1] = static_cast<int>(j);
        indx[nrow][2] = 4;  // indeks przesuniecie w 'y'
        nrow++;
      }
    }
  }

  /*******************************************************************************************************
   * tablice w postaci CSR -  wyznaczamy elementy w wierszu i wpisujemy
   *posortowane do tablicy glownej tablica dla wartosci w pojedynczym wierszu -
   *po wypelnieniu sortujemy
   *******************************************************************************************************/

  std::size_t ncol = 9 * 2;  // maksymalna liczba elementow w wierszu - liczba
                             // sasiadow * liczba kierunkow
  std::vector<double> acol;
  std::vector<int> jcol;

  acol.resize(ncol + 10, 0.);
  jcol.resize(ncol + 10,
              0);  // w zerowym indeksie zapisujemy liczbe elementow w wierszu

  // tablice globalne do rozwiazywania ukladu rownan - allokacja jak w C

  std::size_t nmax =
      nrow * 9 *
      2;  // maksymalna liczba niezerowych elementow w wierszu * liczba wierszy
  double* csr_val = static_cast<double*>(malloc(nmax * sizeof(double)));
  int* icsr = static_cast<int*>(malloc((nrow + 1) * sizeof(int)));
  icsr[nrow] = 0;  // aktualna liczba NNZ	w macierzy ukladu
  int* jcsr = static_cast<int*>(malloc(nmax * sizeof(int)));
  double* ff = static_cast<double*>(malloc(nrow * sizeof(double)));
  double* xx = static_cast<double*>(malloc(nrow * sizeof(double)));
  double* bb = static_cast<double*>(malloc(nrow * sizeof(double)));

  /**
   * tworzymy tablice lokalnego otoczenia punktu 3x3
   *   00  01  02    - numeracja wezlow w otoczeniu wezla (i,j) centralnego (11)
   *   10 (11) 12
   *   20  21  22
   *
   */

  // obsadzenie sasiadow - pij
  std::vector<std::vector<int>> ip;
  ip.resize(3, std::vector<int>(3, 0));

  // rodzaj brzegu: 0-Dirichlet, 1-Neumann
  std::vector<std::vector<int>> iboundary;
  iboundary.resize(3, std::vector<int>(3, 0));

  // tablica oddzialywania - d1
  std::vector<std::vector<double>> d1;
  d1.resize(3, std::vector<double>(3, 0));

  // tablica oddzialywania - d2
  std::vector<std::vector<double>> d2;
  d2.resize(3, std::vector<double>(3, 0));

  /*================================================================================================
   * generujemy elementy macierzowe i wektor wyrazow wolnych
   *================================================================================================*/
  for (std::size_t k = 0; k < nrow; k++) {  // numer wiersza globalnego

    // experimental, wcześniej:
    // int i=indx[k][0]; //atom centralny dla wiersza
    // int j=indx[k][1];
    std::size_t i =
        static_cast<std::size_t>(indx[k][0]);  // atom centralny dla wiersza
    std::size_t j = static_cast<std::size_t>(indx[k][1]);

    for (std::size_t ii = 0; ii < 3; ii++) {
      for (std::size_t jj = 0; jj < 3; jj++) {
        ip[ii][jj] = 0;
        d1[ii][jj] = 0;
        d2[ii][jj] = 0;
        iboundary[ii][jj] = 0;
      }
    }

    //----------wypelniamy lokalne macierze
    // pomocnicze------------------------------------------
    for (std::size_t ii = 0; ii < 3; ii++) {
      for (std::size_t jj = 0; jj < 3; jj++) {
        std::size_t i3 = (i + ii - 1 + nx) % (nx);
        std::size_t j3 = j + jj - 1;
        if (lround(crystal[i3][j3][0]) > 0)
          ip[ii][jj] = 1;  // jest atom
        else
          ip[ii][jj] = 0;  // brak atomu

        if (crystal[i3][j3][3] < 0) {
          iboundary[ii][jj] =
              0;  // brzeg: Dirichlet (wyraz przenosimy do wyrazow wolnych)
        } else {
          iboundary[ii][jj] =
              1;  // brzeg: Neumann (wyrazy zostawiamy w macierzy A)
        }

        long id =
            lround(crystal[i][j][0] * crystal[i3][j3][0]);  // typ oddzialywania
        if (id == 1) {                                      // s-s
          d1[ii][jj] = 0;
          d2[ii][jj] = 0;
        } else if (id == 2 || id == 4) {  // g-g (4) lub g-s (2)
          d1[ii][jj] = ag - as;
          d2[ii][jj] = ag - al;
        }
      }
    }

    /*================================================================================
     * *********** liczymy elementy: A, F ******************************
     * A: format CSR - macierz rzadka (csr_val,icsr,jcsr)
     * F=ff[nrow] - wektor wyrazow wolnych
     *
     *================================================================================*/
    jcol[0] =
        0;  // 0-brak elementow: liczbe elementow trzymamy w elemencie  jcol[0]

    std::size_t number =
        static_cast<std::size_t>(indx[k][2]);  // number:  3-uij, 4-vij

    ff[k] = 0.;  // zerujemy element wektora wyrazow wolnych - usuwamy smieci z
                 // poprzednich iteracji
    fill(acol.begin(), acol.end(), 0.0);
    fill(jcol.begin(), jcol.end(), 0.0);

    if (number == 3) {
      compute_u_v_from_wxx(number, k, i, j, nx, skl, skd, ip, iboundary, d1,
                           crystal, acol, jcol, ff);
      compute_u_v_from_wxy(number, k, i, j, nx, skd, ip, iboundary, d2, crystal,
                           acol, jcol, ff);
    } else if (number == 4) {
      compute_u_v_from_wxx(number, k, i, j, nx, skl, skd, ip, iboundary, d2,
                           crystal, acol, jcol, ff);
      compute_u_v_from_wxy(number, k, i, j, nx, skd, ip, iboundary, d1, crystal,
                           acol, jcol, ff);
    }
    sort_and_add_matrix_elements(nrow, k, jcol, acol, csr_val, icsr, jcsr);
  }  // k=row index

  // rozwiazujemy uklad rownan A*(uv)=ff
  // Conjugate Gradients
  std::size_t itmax0 = *iterations;

  for (std::size_t i = 0; i < nrow; i++)
    xx[i] = 0.0;

  // wektor startowy to poprzednie rozwiazanie
  for (std::size_t k = 0; k < nrow; k++) {  // numer wiersza globalnego
    std::size_t i = static_cast<std::size_t>(indx[k][0]);
    std::size_t j = static_cast<std::size_t>(indx[k][1]);
    std::size_t number = static_cast<std::size_t>(indx[k][2]);  // 3-uij, 4-vij
    xx[k] = crystal[i][j][number - 2];  // number-2: 1-uij, 2-vij
  }

  auto s2 = std::chrono::high_resolution_clock::now();

  // ierr=0,1:
  // 0 - rozwiazujemy uklad rownan
  // 1 - liczymy blad lokalny jak w publikacji
  if (ierr == 0) {  // rozwiazujemy uklad rownan

    solve_linear_system_CG_standard(nrow, csr_val, icsr, jcsr, ff, xx,
                                    iterations, tolerance);

    if (*tolerance >= 1.0E-3 || *iterations >= itmax0) {
      printf("solution:  iterations,  tolerance  =   %6ld   %15.5E  \n\n",
             *iterations, *tolerance);
    }

    // zachowujemy nowe polozenia/przesuniecia atomow
    for (std::size_t k = 0; k < nrow; k++) {  // numer wiersza globalnego
      std::size_t i = static_cast<std::size_t>(indx[k][0]);
      std::size_t j = static_cast<std::size_t>(indx[k][1]);
      std::size_t number = static_cast<std::size_t>(indx[k][2]);  // 3-uij,
                                                                  // 4-vij
      crystal[i][j][number - 2] = xx[k];  // number-2: 1-uij, 2-vij
    }
  }

  // norma max z wektora reszt - liczymy zawsze: ierr-dowolne
  compute_sparse_Ax_y(nrow, csr_val, icsr, jcsr, xx, bb);  // bb = csr_val*xx
  *bmax = 0.;
  for (std::size_t i = 0; i < nrow; i++) {
    bb[i] = bb[i] - ff[i];
    if (std::fabs(bb[i]) > *bmax)
      *bmax = std::fabs(bb[i]);
  }

  auto s3 = std::chrono::high_resolution_clock::now();
  auto s21 =
      std::chrono::duration_cast<std::chrono::microseconds>(s2 - start_time);
  auto s32 = std::chrono::duration_cast<std::chrono::microseconds>(s3 - s2);

  free(csr_val);
  free(icsr);
  free(jcsr);
  free(ff);
  free(xx);
  free(bb);

}  // solve Au=F:end

inline void compute_sparse_Ax_y(const std::size_t& n_rows,
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

inline double scalar_product_x_y(const std::size_t& n, double* x, double* y) {
  double res = 0.;
  for (std::size_t i = 0; i < n; i++) {
    res += x[i] * y[i];
  }
  return res;
}

void solve_linear_system_CG_standard(const std::size_t& n,
                                     double* csr_val,
                                     int* csr_row,
                                     int* csr_col,
                                     double* b,
                                     double* x,
                                     std::size_t* itmax,
                                     double* tolerance) {
  // if b is zero vector, return trivial solution x={0,...}
  double b_2 = scalar_product_x_y(n, b, b);
  if (b_2 < 1.0E-10) {
    for (std::size_t i = 0; i < n; i++) {
      x[i] = 0.;
    }
    *itmax = 0;
    *tolerance = 0.;
    return;
  }

  // residual vectors
  double* rj = static_cast<double*>(malloc(n * sizeof(double)));
  double* rj_proposed = static_cast<double*>(malloc(n * sizeof(double)));
  // approximate solution vectors
  double* xj = static_cast<double*>(malloc(n * sizeof(double)));
  double* xj_proposed = static_cast<double*>(malloc(n * sizeof(double)));
  // vector for storing A*x0
  double* A_times_x0 = static_cast<double*>(malloc(n * sizeof(double)));
  // matrix-vector product result (A*pj)
  double* A_times_pj = static_cast<double*>(malloc(n * sizeof(double)));
  // direction vectors
  double* pj = static_cast<double*>(malloc(n * sizeof(double)));
  double* pj_proposed = static_cast<double*>(malloc(n * sizeof(double)));

  // compute initial guess for A*x0, store in tmp1
  compute_sparse_Ax_y(n, csr_val, csr_row, csr_col, x, A_times_x0);

  for (std::size_t i = 0; i < n; i++) {
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
  for (std::size_t j = 0; j < *itmax; j++) {
    // compute A*pj
    compute_sparse_Ax_y(n, csr_val, csr_row, csr_col, pj, A_times_pj);

    // compute step size alpha = r_j * r_j / (p_j * A * p_j)
    rj_2 = scalar_product_x_y(n, rj, rj);
    Apj_2 = scalar_product_x_y(n, A_times_pj, pj);
    alfa = rj_2 / Apj_2;
    if (fabs(alfa) < 1.0E-5)
      // step too small
      printf(
          "Conjugate Gradient method error, step too small:  alfa= %15.5E \n",
          alfa);

    // update approximate solution vector
    for (std::size_t i = 0; i < n; i++) {
      xj_proposed[i] = xj[i] + alfa * pj[i];
    }

    // update residual vector
    for (std::size_t i = 0; i < n; i++) {
      rj_proposed[i] = rj[i] - alfa * A_times_pj[i];
    }

    // compute the update factor beta = r_{j+1} * r_{j+1} / (r_j * r_j)
    rjp1_2 = scalar_product_x_y(n, rj_proposed, rj_proposed);
    double beta = rjp1_2 / rj_2;
    for (std::size_t i = 0; i < n; i++)
      pj_proposed[i] = rj_proposed[i] + beta * pj[i];

    // insert computed tmp solutions to corresponding vectors
    for (std::size_t i = 0; i < n; i++) {
      xj[i] = xj_proposed[i];
      rj[i] = rj_proposed[i];
      pj[i] = pj_proposed[i];
    }

    // compute solution error
    approximation_error = sqrt(rj_2) / sqrt(b_2);
    // if error is satisfying, end procedure
    if (approximation_error < *tolerance && j > 0) {
      *itmax = j;
      break;
    }
  }
  *tolerance = approximation_error;

  // save computed solution
  for (std::size_t i = 0; i < n; i++)
    x[i] = xj[i];

  // free memory
  free(rj);
  free(rj_proposed);
  free(xj);
  free(xj_proposed);
  free(A_times_x0);
  free(A_times_pj);
  free(pj);
  free(pj_proposed);

  return;
}  // CG-standard

inline void compute_u_v_from_wxx(
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
    double* ff) {
  std::size_t ii = 1;
  std::size_t jj = 1;
  std::size_t lu;
  std::size_t i3;
  double val;

  // element diagonalny

  if (number == 3) {  // uij

    val = -skl * ip[ii][jj] * ip[ii + 1][jj] -
          skl * ip[ii][jj] * ip[ii - 1][jj] -
          skd / 2. * ip[ii][jj] * ip[ii + 1][jj + 1] -
          skd / 2. * ip[ii][jj] * ip[ii - 1][jj - 1] -
          skd / 2. * ip[ii][jj] * ip[ii + 1][jj - 1] -
          skd / 2. * ip[ii][jj] * ip[ii - 1][jj + 1];
    val = val * (-1);  // pochodna wewnetrzna

    lu = static_cast<std::size_t>(jcol[0] + 1);
    jcol[0] = static_cast<int>(lu);
    jcol[lu] = static_cast<int>(lround(crystal[i][j][number]));
    acol[lu] = val;

    // element wolny - wxx
    val = skl * ip[ii][jj] * ip[ii + 1][jj] * d[ii + 1][jj] -
          skl * ip[ii][jj] * ip[ii - 1][jj] * d[ii - 1][jj] +
          skd / 2. * ip[ii][jj] * ip[ii + 1][jj + 1] * d[ii + 1][jj + 1] -
          skd / 2. * ip[ii][jj] * ip[ii - 1][jj - 1] * d[ii - 1][jj - 1] +
          skd / 2. * ip[ii][jj] * ip[ii + 1][jj - 1] * d[ii + 1][jj - 1] -
          skd / 2. * ip[ii][jj] * ip[ii - 1][jj + 1] * d[ii - 1][jj + 1];
    val = val * (-1);  // pochodna wewnetrzna
    ff[k] += val;

  } else if (number == 4) {  // vij

    val = -skl * ip[ii][jj] * ip[ii][jj + 1] -
          skl * ip[ii][jj] * ip[ii][jj - 1] -
          skd / 2. * ip[ii][jj] * ip[ii + 1][jj + 1] -
          skd / 2. * ip[ii][jj] * ip[ii - 1][jj - 1] -
          skd / 2. * ip[ii][jj] * ip[ii + 1][jj - 1] -
          skd / 2. * ip[ii][jj] * ip[ii - 1][jj + 1];
    val = val * (-1);  // pochodna wewnetrzna
    lu = static_cast<std::size_t>(jcol[0] + 1);
    jcol[0] = static_cast<int>(lu);
    jcol[lu] = static_cast<int>(lround(crystal[i][j][number]));
    acol[lu] = val;

    // element wolny - wyy
    val = skl * ip[ii][jj] * ip[ii][jj + 1] * d[ii][jj + 1] -
          skl * ip[ii][jj] * ip[ii][jj - 1] * d[ii][jj - 1] +
          skd / 2. * ip[ii][jj] * ip[ii + 1][jj + 1] * d[ii + 1][jj + 1] -
          skd / 2. * ip[ii][jj] * ip[ii - 1][jj - 1] * d[ii - 1][jj - 1] -
          skd / 2. * ip[ii][jj] * ip[ii + 1][jj - 1] * d[ii + 1][jj - 1] +
          skd / 2. * ip[ii][jj] * ip[ii - 1][jj + 1] * d[ii - 1][jj + 1];
    val = val * (-1);  // pochodna wewnetrzna
    ff[k] += val;
  }

  // elementy pozadiagonalne: horyzontalne (number=3) i wertykalne (number=4)
  if (number == 3) {
    for (int im = -1; im <= 1; im += 2) {
      std::size_t jm = 0;
      std::size_t i_offset =
          static_cast<std::size_t>(static_cast<int>(ii) + im);
      std::size_t j_offset = jj + jm;
      i3 = static_cast<std::size_t>(
          (static_cast<int>(i) + im + static_cast<int>(nx)) %
          static_cast<int>(nx));
      val = skl * ip[ii][jj] * ip[i_offset][j_offset];
      val = val * (-1);  // pochodna wewnetrzna
      if (iboundary[i_offset][j_offset] == 1) {
        lu = static_cast<std::size_t>(jcol[0]) + 1;
        jcol[0] = static_cast<int>(lu);
        jcol[lu] = static_cast<int>(lround(crystal[i3][j + jm][number]));
        acol[lu] = val;
      } else if (iboundary[i_offset][j_offset] == 0)
        ff[k] -= val * crystal[i3][j + jm][number - 2];
    }

  } else if (number == 4) {
    for (int jm = -1; jm <= 1; jm += 2) {
      std::size_t im = 0;
      std::size_t i_offset = ii + im;
      std::size_t j_offset =
          static_cast<std::size_t>(static_cast<int>(jj) + jm);
      i3 = (i + im + nx) % nx;
      val = skl * ip[ii][jj] * ip[i_offset][j_offset];
      val = val * (-1);  // pochodna wewnetrzna
      if (iboundary[i_offset][j_offset] == 1) {
        lu = static_cast<std::size_t>(jcol[0] + 1);
        jcol[0] = static_cast<int>(lu);
        jcol[lu] = static_cast<int>(lround(crystal[i3][static_cast<std::size_t>(
            static_cast<int>(j) + jm)][number]));
        acol[lu] = val;
      } else if (iboundary[i_offset][j_offset] == 0)
        ff[k] -= val *
                 crystal[i3][static_cast<std::size_t>(static_cast<int>(j) + jm)]
                        [number - 2];
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
          (static_cast<int>(i) + im + static_cast<int>(nx)) %
          static_cast<int>(nx));
      val = skd / 2. * ip[ii][jj] * ip[i_offset][j_offset];
      val = val * (-1);                          // pochodna wewnetrzna
      if (iboundary[i_offset][j_offset] == 1) {  // Neumann
        lu = static_cast<std::size_t>(jcol[0] + 1);
        jcol[0] = static_cast<int>(lu);
        jcol[lu] = static_cast<int>(lround(crystal[i3][static_cast<std::size_t>(
            static_cast<int>(j) + jm)][number]));
        acol[lu] = val;
      } else if (iboundary[i_offset][j_offset] == 0)  // Dirichlet
        ff[k] -= val *
                 crystal[i3][static_cast<std::size_t>(static_cast<int>(j) + jm)]
                        [number - 2];
    }
  }

  return;
}  // compute_u_v_from_wxx

inline void compute_u_v_from_wxy(
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
    double* ff) {
  std::size_t ii = 1;
  std::size_t jj = 1;
  double wsp = 2.0;  // mnoznik dla wxy w wij
  std::size_t lu;
  std::size_t i3;
  double val;

  std::size_t number_2;

  if (number == 3) {
    number_2 = 4;  // indeks dla elementu v
  } else if (number == 4) {
    number_2 = 3;  // indeks dla elementu u
  }

  for (int im = -1; im <= 1; im += 2) {
    for (int jm = -1; jm <= 1; jm += 2) {
      int sign = im * jm * (-1);
      std::size_t i_offset =
          static_cast<std::size_t>(static_cast<int>(ii) + im);
      std::size_t j_offset =
          static_cast<std::size_t>(static_cast<int>(jj) + jm);
      i3 = static_cast<std::size_t>(
          (static_cast<int>(i) + im + static_cast<int>(nx)) %
          static_cast<int>(nx));
      double val_local =
          sign * skd / 4. * ip[ii][jj] * ip[i_offset][j_offset] * wsp;
      if (iboundary[i_offset][j_offset] == 1) {
        lu = static_cast<std::size_t>(jcol[0] + 1);
        jcol[0] = static_cast<int>(lu);
        jcol[lu] = static_cast<int>(lround(
            crystal[i3][static_cast<std::size_t>(static_cast<int>(j) + jm)]
                   [number_2]));  // oddzialywanie: u->v, v->u
        acol[lu] = val_local;
      } else if (iboundary[i_offset][j_offset] == 0)
        ff[k] -= val_local *
                 crystal[i3][static_cast<std::size_t>(static_cast<int>(j) + jm)]
                        [number_2 - 2];
    }
  }

  // element: vij*uij  - do diagonali w csr_val
  val = (skd / 4. * ip[ii][jj] * ip[ii - 1][jj - 1] +
         skd / 4. * ip[ii][jj] * ip[ii + 1][jj + 1] -
         skd / 4. * ip[ii][jj] * ip[ii + 1][jj - 1] -
         skd / 4. * ip[ii][jj] * ip[ii - 1][jj + 1]) *
        wsp;

  lu = static_cast<std::size_t>(jcol[0] + 1);
  jcol[0] = static_cast<int>(lu);
  i3 = (i + nx) % (nx);
  jcol[lu] = static_cast<int>(lround(crystal[i3][j][number_2]));  // v
  acol[lu] = val;

  // element wolny: f(k)  - wxy

  if (number == 3) {
    val = (skd / 4. * ip[ii][jj] * ip[ii - 1][jj - 1] * d[ii - 1][jj - 1] -
           skd / 4. * ip[ii][jj] * ip[ii + 1][jj + 1] * d[ii + 1][jj + 1] -
           skd / 4. * ip[ii][jj] * ip[ii + 1][jj - 1] * d[ii + 1][jj - 1] +
           skd / 4. * ip[ii][jj] * ip[ii - 1][jj + 1] * d[ii - 1][jj + 1]) *
          wsp;
    ff[k] = ff[k] + val;
  } else if (number == 4) {
    val = (skd / 4. * ip[ii][jj] * ip[ii - 1][jj - 1] * d[ii - 1][jj - 1] -
           skd / 4. * ip[ii][jj] * ip[ii + 1][jj + 1] * d[ii + 1][jj + 1] +
           skd / 4. * ip[ii][jj] * ip[ii + 1][jj - 1] * d[ii + 1][jj - 1] -
           skd / 4. * ip[ii][jj] * ip[ii - 1][jj + 1] * d[ii - 1][jj + 1]) *
          wsp;
    ff[k] = ff[k] + val;
  }

  return;
}  // compute_u_v_from_wxy

inline void sort_and_add_matrix_elements(const std::size_t& nrow,
                                         const std::size_t& k,
                                         std::vector<int>& jcol,
                                         std::vector<double>& acol,
                                         double* csr_val,
                                         int* icsr,
                                         int* jcsr) {
  // sorting
  std::size_t l =
      static_cast<std::size_t>(jcol[0]);  // number of non-zero elements in row
  int l1, l2;
  double a1, a2;

  for (std::size_t i = 1; i < l; i++) {
    for (std::size_t j = i; j >= 1; j--) {
      l1 = jcol[j];  // numery kolumn
      l2 = jcol[j + 1];
      if (l1 > l2) {  // zamieniamy miejscami
        a1 = acol[j];
        a2 = acol[j + 1];
        acol[j] = a2;
        acol[j + 1] = a1;
        jcol[j] = l2;
        jcol[j + 1] = l1;
      }
    }
  }

  if (l < 1) {
    printf("brak elementow w wierszu macierzy\n");
    exit(0);
  }

  // sprawdzamy numery kolumn
  for (std::size_t i = 1; i < l; i++) {
    if (jcol[i] >= jcol[i + 1]) {
      printf("blad w numerach kolumn:  %d   %d\n ",
             jcol[static_cast<std::size_t>(i)],
             jcol[static_cast<std::size_t>(i) + 1]);
      exit(0);
    }
  }

  // dodajemy elementy do macierzy: csr_val, jcsr, icsr
  int nnz =
      icsr[nrow];  // aktualna liczba elementow niezerowych - indeksowane od 0,
  icsr[k] = nnz;   // pozycja nnz jest pusta - od niej zaczynamy wypelnianie
                   // wiersza k-tego
  for (std::size_t i = 1; i <= l; i++) {
    csr_val[nnz] = acol[i];
    jcsr[nnz] = jcol[i];
    nnz++;
  }
  icsr[nrow] = nnz;  // zachowujemy aktualna wartosc nnz

  return;
}  // sort_and_add
