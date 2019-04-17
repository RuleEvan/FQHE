#include "hilbert.h"

void setup_hamiltonian(double* ham) {
  int n_states = 2*M_S + 1;
  int L_landau = M_S; 

  ham = (double*) malloc(sizeof(double)*pow(n_states, 2));
  

}

void build_two_body_jumps(int l) {
  FILE* out_file;
  out_file = fopen("two_body_jumps.dat", "w");
  int dim = 2*l + 1;
  for (int m1 = -l; m1 <= l; m1++) {
    for (int m2 = -l; m2 <= l; m2++) {
      for (int m3 = -l; m3 <= l; m3++) {
        for (int m4 = -l; m4 <= l; m4++) {
          if (m1 + m2 != m3 + m4) {continue;}
          gsl_complex jump = double_lebedev(wigner_coulomb, l, l, l, l, m1, m2, m3, m4);
          fprintf(out_file, "%d %d %d %d %d %g %g\n", l, m1, m2, m3, m4, GSL_REAL(jump), GSL_IMAG(jump));
          printf("%d %d %d %d %d %g %g\n", l, m1, m2, m3, m4, GSL_REAL(jump), GSL_IMAG(jump));
        }
      }
    }
  }
  fclose(out_file);
  return;
}

void generate_interaction_file(int l) {

  for (int j = 0; j <= 2*l; j++) {
    double h = 0;
    if ((2*l + j) % 2 == 0) {continue;}
    for (int n = 0; n <= 2*l; n++) {
      h += pow(-1.0, j)/sqrt(l)*pow(2*l + 1, 2)*gsl_sf_coupling_6j ((int) 2*l, (int) 2*l, (int) 2*n, (int) 2*l, (int) 2*l, (int) 2*j)*pow(gsl_sf_coupling_3j((int) 2*l, (int) 2*n, (int) 2*l, (int) -2*l, 0, (int) 2*l), 2.0);
    }
    
//    printf("%d %d %d %d %d %d %g\n", 1, 1, 1, 1, j, 0, );
    printf("%d %d %d %d %d %d %g\n", 1, 1, 1, 1, j, 1, h);
  }

  return;
}
