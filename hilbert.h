#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H
#include "lanczos.h"

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
void hierarcy(double l, double s, int n);
int perm_sign(int* perm, int n);
#endif
