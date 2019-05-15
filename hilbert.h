#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H
#include "lanczos.h"

typedef struct slater_det
{
  int np;
  int* m_val;
} slater_det;

typedef struct wfn
{
  int np;
  slater_det* basis;
  double* coeff;
} wfn;

void build_two_body_jumps(int l);
void generate_wave_function();
void generate_interaction_file(double l); 
void compute_1body_energy(int np, double m, double e_shift);
double xme(int m1, int m1p, int m2, int m2p);
gsl_complex u_spin(double m, double theta, double phi);
gsl_complex u_spin_l(double l, double m, double theta, double phi);
void basis_coeff(int m, int n);
void excited_state(int m, int n);
void laughlin(int m, int n);
void hierarchy_mult(double l, double s, int n);
void hierarchy(double l, double s, int n);
int perm_sign(int* perm, int n);
void generate_rs_matrix(double mtot1, double mtot2, double mtot3);
void order_perm(int** perm, int n);
void perm_compress(int* perm, int** result, int n);
unsigned int next_perm(unsigned int t);
void increment(unsigned int** perm, unsigned int* perm_min, unsigned int* perm_max, double** s_ar, int n, int i, int* done);
#endif
