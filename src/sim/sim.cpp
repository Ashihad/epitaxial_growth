#include "sim.hpp"

#include "physics.hpp"
#include "rng.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace std;

/********************************************************************************************************************************
 ********************************************************************************************************************************
 *
 * 				procedura symulacji wzrostu krzystalu: 1+1
 *
 *
 ********************************************************************************************************************************
 ********************************************************************************************************************************/

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
                       int n_global_relaxation) {
  printf(
      "=============================================================== \n\n");
  printf("as [A] =   %15.5f \n", as);
  printf("ag [A] =   %15.5f \n", ag);
  printf("al [A] =   %15.5f \n", al);
  printf("skl [eV/A/A] =  %15.5f \n", skl);
  printf("skd [eV/A/A] =  %15.5f \n", skd);
  printf("mu [A] =   %15.5f \n", mu);
  printf("D [1/s]=    %15.5E \n", D);
  printf("E [eV] =    %15.5f \n", E);
  printf("gamma [eV] = %15.5f \n", gamma);
  printf("nx=    %ld \n", nx);
  printf("ny=    %ld \n", ny);
  printf("nsurf= %ld \n", nsurf);
  printf("irange_min = %ld \n", irange_min);
  printf("irange_max = %ld \n", irange_max);
  printf("time_fluence [s] = %15.3f \n", time_fluence);
  printf("time_simulation [s] = %15.3f \n", time_simulation);
  printf("k_max_step= %ld \n", k_max_step);
  printf("f_deposition [ML/s]= %15.3f \n", f_deposition);
  printf("temperature [K] = %15.0f \n", temperature);
  printf(
      "================================================================ \n\n");

  /***************************************************************************************************
   *  tablice z polozeniami atomow - polozenie (x,y)
   *  odleglosci atomowe zrenormalizowane - oddzialywania/naprezenia skalowane
   *wzgledem typu atomow
   *
   ***************************************************************************************************/
  vector<vector<vector<double>>> crystal;
  crystal.resize(nx, vector<vector<double>>(ny, vector<double>(10, 0.)));

  vector<vector<vector<double>>> crystal_copy;
  crystal_copy.resize(nx, vector<vector<double>>(ny, vector<double>(10, 0.)));

  // ustawiamy atomy podloza
  for (std::size_t i = 0; i < nx; i++) {
    for (std::size_t j = 0; j < nsurf; j++) {
      crystal[i][j][0] = 1.0;  // 0-empty, 1-Si, 2-Ge
      crystal[i][j][1] = 0.0;  // u
      crystal[i][j][2] = 0.0;  // v
    }
  }

  std::size_t istart = 0;        //(nx/2)-10;
  std::size_t ikoniec = nx - 1;  //(nx/2+10);

  // kladziemy 1 monowarstwe atomow typu-2
  for (std::size_t i = istart; i <= ikoniec; i++) {
    std::size_t j = nsurf;

    crystal[i][j][0] = 2.0;  // 0-empty, 1-Si, 2-Ge
    crystal[i][j][1] = 0.0;  // u
    crystal[i][j][2] = 0.0;  // v

    j = nsurf + 1;
    crystal[i][j][0] = 2.0;  // 0-empty, 1-Si, 2-Ge
    crystal[i][j][1] = 0.0;  // u
    crystal[i][j][2] = 0.0;  // v

    j = nsurf + 2;
    crystal[i][j][0] = 2.0;  // 0-empty, 1-Si, 2-Ge
    crystal[i][j][1] = 0.0;  // u
    crystal[i][j][2] = 0.0;  // v

    j = nsurf + 3;
    crystal[i][j][0] = 2.0;  // 0-empty, 1-Si, 2-Ge
    crystal[i][j][1] = 0.0;  // u
    crystal[i][j][2] = 0.0;  // v

    j = nsurf + 4;
    crystal[i][j][0] = 2.0;  // 0-empty, 1-Si, 2-Ge
    crystal[i][j][1] = 0.0;  // u
    crystal[i][j][2] = 0.0;  // v
    /*
     */
  }

  /*

                  FILE *fpp;
                  fpp=fopen("atoms_zle_old.dat","r");
                  for(int j=0;j<ny;j++){
                          for(int i=0;i<nx;i++){
                                  double a0,a1,a2,a3,a4,a5,a6;
                                  fscanf(fpp,"%lf %lf  %lf  %lf  %lf  %lf  %lf
     ",&a0,&a1,&a2,&a3,&a4,&a5,&a6); crystal[i][j][0]=a0; crystal[i][j][1]=a1;
                                  crystal[i][j][2]=a2;
                                  crystal[i][j][3]=a3;
                                  crystal[i][j][4]=a4;
                                  crystal[i][j][5]=a5;
                                  crystal[i][j][6]=a6;
                          }
                  }
                  fclose(fpp);
                  printf("WCZYTANE\n");

  */

  // ewolucja
  double time{};
  double dt;

  double ev{1.602E-19};  // 1-elektronowolt
  double k_boltzmann{1.38E-23};
  double kbt{k_boltzmann * temperature / ev};  // energia termiczna w [eV]

  double fluence;
  double rdep;
  std::size_t n_at_diff;
  double ep, ri, r0;
  int neigh;
  int number_atoms = 0;
  double n_proposed = 1.;  // liczba wylosowanych procesow dyfuzji
  double n_accepted = 1.;  // liczba zaakceptowanych procesow dyfuzji

  vector<vector<double>> atoms_diff;  // tu trzymac bedziemy informacje o
                                      // atomach podlegajacych dyfuzji
  atoms_diff.resize(nx + 1, vector<double>(10, 0));

  int iter{};
  auto czas_start = std::chrono::high_resolution_clock::now();

  while (time < time_simulation) {
    iter++;
    // zapis do pliku co n_write_files iteracji
    if (iter % n_write_files == 0) {
      // liczymy energie sperezystosci dla kazdego atomu
      for (std::size_t i = 0; i < nx; i++) {
        for (std::size_t j = 1; j < ny - 1; j++) {
          crystal[i][j][8] =
              compute_elastic_energy_wij(i, j, crystal, as, ag, al, skl, skd);
        }
      }

      FILE* fp;
      fp = fopen("atoms.dat", "w");
      for (std::size_t i = 0; i < nx; i++) {
        for (std::size_t j = nsurf - 1; j < ny; j++) {
          if (crystal[i][j][0] > 0.5)
            fprintf(fp, "%8ld %8ld  %8ld    %15.5E\n", i, j,
                    lround(crystal[i][j][0]), crystal[i][j][8]);
        }
      }
      fclose(fp);

      fp = fopen("elastic_en.dat", "w");
      for (std::size_t i = 0; i < nx; i++) {
        for (std::size_t j = 1; j < ny - 1; j++) {
          double grad_loc_2 =
              pow(crystal[i][j][6], 2) + pow(crystal[i][j][7], 2);
          fprintf(fp, "%8ld %8ld  %8ld    %15.5E   %15.5E\n", i, j,
                  lround(crystal[i][j][0]), crystal[i][j][8], grad_loc_2);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);

      // dane tymczasowe
      FILE* fpp;
      fpp = fopen("atoms_data.dat", "w");
      for (std::size_t j = 0; j < ny; j++) {
        for (std::size_t i = 0; i < nx; i++) {
          fprintf(fpp, "%15.5E  ", crystal[i][j][0]);
          fprintf(fpp, "%15.5E  ", crystal[i][j][1]);
          fprintf(fpp, "%15.5E  ", crystal[i][j][2]);
          fprintf(fpp, "%15.5E  ", crystal[i][j][3]);
          fprintf(fpp, "%15.5E  ", crystal[i][j][4]);
          fprintf(fpp, "%15.5E  ", crystal[i][j][5]);
          fprintf(fpp, "%15.5E  ", crystal[i][j][6]);
          fprintf(fpp, "\n");
        }
      }
      fclose(fpp);

      //********* INFO ************************************************
      auto czas_teraz = std::chrono::high_resolution_clock::now();
      auto czas = std::chrono::duration_cast<std::chrono::milliseconds>(
          czas_teraz - czas_start);

      double duv_max = 0.;
      for (std::size_t j = 0; j < ny; j++) {
        for (std::size_t i = 0; i < nx; i++) {
          double duv = pow(crystal[i][j][1], 2) + pow(crystal[i][j][2], 2);
          duv = sqrt(duv);
          if (duv > duv_max && crystal[i][j][0] > 0.5)
            duv_max = duv;
        }
      }

      printf(
          "iter, time, new_atoms, pr_acc, (real time [s]), duv_max =    %6d    "
          "%12.6f  %6d  %12.4f   (%12.1f)    %10.3f\n",
          iter, time, number_atoms, n_accepted / n_proposed,
          static_cast<double>(czas.count()) * 1.0E-3, duv_max);
    }

    /******************************************************************************
     * globalna relaksacja naprezen w sieci
     ******************************************************************************/
    if (iter % n_global_relaxation == 0 || iter == 1) {
      int iloc = 0;  // relaksacja globalna
      std::size_t iter0 = iterations;
      double tol0 = tolerance;
      conduct_relaxation(0, 0, irange_min, irange_max, &iter0, &tol0, crystal,
                         tol_local_max, iloc, as, ag, al, skl, skd);
    }

    /*********************************************************************************************************************
     * okreslamy tempo depozycji atomow: deposition rate
     * rate= fluence*nx -> tempo (prawdopodobienstwo/sekunde) depozycji jednego
     *atomu w calym ukladzie (gdziekolwiek)
     *
     *********************************************************************************************************************/
    if (time < time_fluence)
      fluence = f_deposition * static_cast<double>(nx);
    else
      fluence = 0.;
    rdep = (static_cast<double>(k_max_step) + 1.0) *
           (2 * static_cast<double>(k_max_step) + 1.0) * fluence / 6.;
    atoms_diff[0][5] = rdep;  // tempo depozycji atomow z wiazki
    atoms_diff[0][6] = rdep;

    /**********************************************************************************************************************
     * szukamy atomow powierzchniowych podlegajacych dyfuzji -> tworzymy liste,
     *z ktorej wybierzemy jeden lub depozycje idiffusion=1,2:   1-atomy podloza
     *i zdeponowane, 2-tylko zdeponowane
     *
     **********************************************************************************************************************/
    n_at_diff = 0;  // liczba atomow powierzchniowych podlegajacych dyfuzji
    for (std::size_t i = 0; i < nx; i++) {
      for (std::size_t j = ny - 1; j >= 1; j--) {
        long int ll =
            lround(crystal[i][j][0]);  // ll=0(brak), 1(podloze), 2(atom wiazki)
        if (ll > 0) {
          if (ll >= idiffusion) {
            n_at_diff++;
            atoms_diff[n_at_diff][0] =
                static_cast<double>(lround(crystal[i][j][0]));  // typ atomu
            atoms_diff[n_at_diff][1] = static_cast<double>(i);  // polozenie x
            atoms_diff[n_at_diff][2] = static_cast<double>(j);  // polozenie y

            neigh = 0;  // sasiedzi: NN i NNN (NN=poziomo+pionowo,
                        // NNN=diagonala+antydiagonala)
            for (int ii = -1; ii <= 1; ii++) {
              for (int jj = -1; jj <= 1; jj++) {
                if ((abs(ii) + abs(jj)) > 0) {
                  std::size_t i2 = static_cast<std::size_t>(
                      (static_cast<int>(i) + ii + static_cast<int>(nx)) %
                      static_cast<int>(nx));
                  std::size_t j2 = static_cast<std::size_t>(
                      j + static_cast<std::size_t>(jj));
                  if (crystal[i2][j2][0] > 0.5)
                    neigh++;  // jest sasiad - dodajemy
                }
              }
            }
            atoms_diff[n_at_diff][3] = neigh;
            ep = compute_elastic_energy_wij(i, j, crystal, as, ag, al, skl,
                                            skd);  // local elastic energy
            atoms_diff[n_at_diff][4] = ep;

            r0 = 12. * D / (static_cast<double>(k_max_step) + 1.0) /
                 (2 * static_cast<double>(k_max_step) + 2.);
            double delta_w = ep;  // neigh<=2
            if (neigh == 3) {
              delta_w = ep * 1.5;
            } else if (neigh == 4) {
              delta_w = ep * 2.0;
            } else if (neigh >= 5) {
              delta_w = ep * 3.5;
            }
            ri = r0 * exp((-gamma * neigh + delta_w + E) /
                          kbt);  // szacowane tempo dyfuzji atomu
            atoms_diff[n_at_diff][5] = ri;
            atoms_diff[n_at_diff][6] =
                ri;  // kopia dla porownania dokladnego prawdopodobienstwa
            atoms_diff[n_at_diff][7] = neigh;
            atoms_diff[n_at_diff][8] = delta_w;
          }
          break;  // przerywamy sprawdzanie - natrafilismy na atom idac od gory
        }
      }
    }

    // krok czasowy -> zmiana czasu [odwrotnosc aktualnej sumy czestosci
    // procesow dt=1/sum_{i}(Gamma_i)]
    double sum_ri = 0;
    for (std::size_t i = 0; i <= n_at_diff; i++) {
      sum_ri += atoms_diff[i][5];
    }
    dt = 1.0 / sum_ri;
    time += dt;

    // losujemy proces: depozycja lub dyfuzja atomu
    for (std::size_t i = 1; i <= n_at_diff; i++) {
      atoms_diff[i][5] +=
          atoms_diff[i - 1][5];  // suma Gamma_i do generatora dyskretnego
    }

    double u1 = gen_uniform() * sum_ri;
    std::size_t ktory =
        n_at_diff;  // zabezpieczenie na wypadek gdyby petla nie zadziala
    for (std::size_t i = 0; i <= n_at_diff; i++) {
      if (u1 <= atoms_diff[i][5]) {
        ktory = i;
        break;
      }
    }

    if (ktory == 0) {
      /***************************************************************************************************
       * 			   ktory=0:  	losowa depozycja atomu
       *
       ***************************************************************************************************/
      number_atoms++;  // zwiekszamy liczbe atomow w ukladzie

      std::size_t ipos =
          gen_discrete_1_K(nx) - 1;  // polozenie atomu w kierunku x: 0-(nx-1)
      std::size_t jpos = 0;
      for (std::size_t j = ny - 1; j >= 1; j--) {
        long int ll = lround(crystal[ipos][j][0]);
        if (ll > 0) {
          jpos = j + 1;
          if (jpos < ny) {
            crystal[ipos][jpos][0] = 2.;  // dodajemy atom typu 2
            crystal[ipos][jpos][1] = 0.;  // u
            crystal[ipos][jpos][2] = 0.;  // v
            break;
          } else {
            printf("za duzo atomow w tablicy: jpos=ny\n");
            printf("omijamy punkt i=%ld\n", ipos);
            break;
          }
        }
      }

      std::size_t iter0 = iterations;
      double tol0 = tolerance;
      int iloc = 1;  // relaksacja lokalna
      conduct_relaxation(ipos, jpos, irange_min, irange_max, &iter0, &tol0,
                         crystal, tol_local_max, iloc, as, ag, al, skl, skd);

    }  // ktory==0: depozycja atomu
    else {
      /****************************************************************************************************
       * 				ktory>0: dyfuzja losowego atomu
       *
       ****************************************************************************************************/

      std::size_t ii_shift = gen_discrete_1_K_multiply_sign(
          k_max_step);  // losowe przesuniecie lewo-prawo
      std::size_t i_old = static_cast<std::size_t>(lround(
          atoms_diff[ktory][1]));  // aktualna pozycja atomu dyfundujacego
      std::size_t j_old = static_cast<std::size_t>(
          lround(atoms_diff[ktory][2]));                   // aktualna pozycja
      std::size_t i_new = (i_old + ii_shift + nx) % (nx);  // nowa pozycja atomu
      double typ;
      std::size_t j_new;

      std::size_t iter0;
      double tol0;

      int irange = 25;

      auto s1 = std::chrono::high_resolution_clock::now();
      //**** liczymy stara energie [-irange,irange]:  atom-on
      //*****************************
      double en_old = 0.;
      for (int i = -irange; i <= irange; i++) {
        for (int j = -irange; j <= irange; j++) {
          std::size_t ii = static_cast<std::size_t>(
              (static_cast<int>(i_old) + i + static_cast<int>(nx)) %
              static_cast<int>(nx));
          std::size_t jj =
              static_cast<std::size_t>(static_cast<int>(j_old) + j);
          if (jj > 1 && jj < (ny - 1) && crystal[ii][jj][0] > 0.5) {
            en_old += compute_elastic_energy_wij(ii, jj, crystal, as, ag, al,
                                                 skl, skd);  // elastic energy
          }
        }
      }

      auto s2 = std::chrono::high_resolution_clock::now();
      //**** liczymy energie po usunieciu atomu (i_old,j_old):   atom-off
      //*************
      crystal_copy = crystal;              // kopia: atom-on
      crystal_copy[i_old][j_old][0] = 0.;  // usuwamy atom

      //***** relaksacja naprezen w sieci atom-off
      //*************************************
      iter0 = iterations;
      tol0 = tolerance;
      int iloc = 1;  // relaksacja lokalna
      conduct_relaxation(i_old, j_old, irange_min, irange_max, &iter0, &tol0,
                         crystal_copy, tol_local_max, iloc, as, ag, al, skl,
                         skd);

      auto s3 = std::chrono::high_resolution_clock::now();
      //**** liczymy nowa energie - atom-off
      //[-range,irange]*****************************
      double en_new = 0.;
      for (int i = -irange; i <= irange; i++) {
        for (int j = -irange; j <= irange; j++) {
          int ii =
              static_cast<int>((i_old + static_cast<std::size_t>(i) + nx) % nx);
          int jj = static_cast<int>(j_old + static_cast<std::size_t>(j));
          if (jj > 1 && static_cast<std::size_t>(jj) < (ny - 1) &&
              crystal_copy[static_cast<std::size_t>(ii)]
                          [static_cast<std::size_t>(jj)][0] > 0.5) {
            en_new += compute_elastic_energy_wij(
                static_cast<std::size_t>(ii), static_cast<std::size_t>(jj),
                crystal_copy, as, ag, al, skl, skd);  // elastic energy
          }
        }
      }

      auto s4 = std::chrono::high_resolution_clock::now();
      auto s21 = std::chrono::duration_cast<std::chrono::microseconds>(s2 - s1);
      auto s32 = std::chrono::duration_cast<std::chrono::microseconds>(s3 - s2);
      auto s43 = std::chrono::duration_cast<std::chrono::microseconds>(s4 - s3);
      // printf("%15.5E   %15.5E   %15.5E
      // %d\n",s21*1.0E-6,s32*1.0E-6,s43*1.0E-6,iter0);

      //******* sprawdzamy czy prawdopodobienstwo r_new<r_old=r_approx - jesli
      // tak to atom dyfunduje   ***********************
      double neigh_d = std::round(
          atoms_diff[ktory][7]);  // aktualna liczba najblizszych sasiadow
      r0 = 12. * D / (static_cast<double>(k_max_step) + 1.0) /
           (2 * static_cast<double>(k_max_step) + 2.);
      double ri_new = r0 * exp((-gamma * neigh_d + (en_old - en_new) / 2 + E) /
                               kbt);  // dokladne tempo dyfuzji atomu
      double ri_old =
          atoms_diff[ktory][6];  // szacowane prawodpodobienstwo dyfuzji

      n_proposed++;

      if (gen_uniform() < (ri_new / ri_old)) {
        n_accepted++;
        typ = crystal[i_old][j_old][0];  // zachowujemy typ atomu
        crystal[i_old][j_old][0] = 0;    // kasujemy atom w starej pozycji

        crystal = crystal_copy;
        // osadzamy atom w j_new
        for (std::size_t j = ny - 1; j >= 1; j--) {
          long ll = lround(crystal[i_new][j][0]);
          if (ll > 0) {
            if (j < (ny - 1)) {
              j_new = j + 1;  // tu przesuwamy atom na wolne miejsce
              crystal[i_new][j_new][0] = typ;
              break;
            } else {
              printf(
                  "nie mozemy przesunac atomu do nowej pozycji - brak miejsca "
                  "w tabilicy \n\n");
              break;
            }
          }
        }
        // lokalna relaksacja polozen atomow w starym i w nowym polozeniu
        iter0 = iterations;
        tol0 = tolerance;
        iloc = 1;  // relaksacja lokalna
        conduct_relaxation(i_new, j_new, irange_min, irange_max, &iter0, &tol0,
                           crystal, tol_local_max, iloc, as, ag, al, skl, skd);
      }

    }  // ktory>0: dyfuzja

  }  // time<time_simulation

}  // growth_simulation
