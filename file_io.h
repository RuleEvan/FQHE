#ifndef FILE_IO_H
#define FILE_IO_H
#include "angular.h"

typedef struct slater_det
{
  float* m_val;
} slater_det;

typedef struct wfnData
{
  int n_spin_up, n_spin_down;
  slater_det** basis;
  long long int n_states;
  int n_shells, n_orbits, n_data;
  int n_eig;
  float b_mag;
  unsigned int n_sds_up, n_sds_down;
  double *bc;
  int *n_shell, *l_shell, *j_shell, *jz_shell, *tz_shell;
  int *n_orb, *l_orb;
  float *j_orb;
  float jz;
  float *e_nuc, *j_nuc, *t_nuc;
} wfnData;

wfnData* read_binary_wfn_data(char *wfn_file, char* basis_file, int i_eig);
int perm_sign(int* perm, int n);
void order_perm(int** perm, int n);
void perm_compress(int* perm, int** result, int n);
unsigned int next_perm(unsigned int t);
void increment(unsigned int** perm, unsigned int* perm_min, unsigned int* perm_max, double** s_ar, int n, int i, int* done);

#endif
