#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H
#include "lanczos.h"


void generate_interaction_file(double l); 
void compute_1body_energy(int np, double m, double e_shift);
double xme(int m1, int m1p, int m2, int m2p);
void laughlin(int m, int n);
void hierarchy(double l, double s, int n);
wfnData* shift_op(wfnData* wd, double m);
void normalize_wfn(wfnData* wd);
void print_wfn(wfnData* wd);
#endif
