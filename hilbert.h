#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H
#include "lanczos.h"


void generate_interaction_file(double q, int n, int iv); 
void compute_1body_energy(int np, double m, double e_shift);
double xme(int m1, int m1p, int m2, int m2p);
void laughlin(int m, int n);
wfnData* hierarchy(double l, double s, int n);
wfnData* composite(double l, double s, int n);
wfnData* deform_hierarchy(double l, double s, int n);
wfnData* shift_op(wfnData* wd, double m, int i_type);
void normalize_wfn(wfnData* wd);
void print_wfn(wfnData* wd);
void wfn_sum(wfnData* wd1, wfnData* wf2, double fact);
double wfn_dot(wfnData* wd1, wfnData* wf2);
void wfn_mult(wfnData* wd, double fact);
wfnData* u_dot_d(wfnData* wd);
wfnData* u_tensor_d(wfnData* wd, double j_tot, double m_tot);
wfnData* u_tensor_u(wfnData* wd, double j_tot, double m_tot);
wfnData* d_tensor_d(wfnData* wd, double j_tot, double m_tot);
wfnData* uu_tensor_dd(wfnData* wd, double j_u, double j_d, double j_tot, double m_tot);
wfnData* ud_tensor_ud(wfnData* wd, double j_u, double j_d, double j_tot, double m_tot);
wfnData* order_wfn(wfnData* wf);
wfnData* coherent_shift(wfnData* wd, double theta, double phi, int iv);
void generate_multilevel_interaction_file(double m, int iv);
#endif
