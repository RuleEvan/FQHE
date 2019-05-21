#include "hilbert.h"

wfnData* u_dot_d(wfnData* wd) {
  // Acts the dot product of the raising shift operator
  // with the lowering shift operator on the given wavefunction
  int n = wd->n_spin_up;
  double m_max = n/2.0;
  wfnData* wd_f = shift_op(shift_op(wd, m_max, 1), -m_max, 0);
  for (int k = 0; k < n; k++) {
    double m = k - m_max;
    wfn_sum(wd_f, shift_op(shift_op(wd, m, 1), -m, 0), pow(-1, (n/2.0 - m)));
  }

  return wd_f;
}

wfnData* ud_tensor_ud(wfnData* wd, double j_u, double j_d, double j_tot, double m_tot) {
  wfnData* wd_f;
  int done_one = 0;
  for (int im1 = 0; im1 < 2*j_u + 1; im1++) {
    double m_u = im1 - j_u;
    double m_d = m_tot - m_u;
    if (fabs(m_d) > j_d) {continue;}
    if (clebsch_gordan(j_u, j_d, j_tot, m_u, m_d, m_tot) == 0) {continue;}
    if (done_one == 0) {
      wd_f = u_tensor_d(u_tensor_d(wd, j_d, m_d), j_u, m_u);
      wfn_mult(wd_f, clebsch_gordan(j_u, j_d, j_tot, m_u, m_d, m_tot));
      done_one = 1;
    } else {
      wfn_sum(wd_f, u_tensor_d(u_tensor_d(wd, j_d, m_d), j_u, m_u), clebsch_gordan(j_u, j_d, j_tot, m_u, m_d, m_tot));
    }
  }
  return wd_f;
}

wfnData* uu_tensor_dd(wfnData* wd, double j_u, double j_d, double j_tot, double m_tot) {
  wfnData* wd_f;
  int done_one = 0;
  for (int im1 = 0; im1 < 2*j_u + 1; im1++) {
    double m_u = im1 - j_u;
    double m_d = m_tot - m_u;
    if (fabs(m_d) > j_d) {continue;}
    if (clebsch_gordan(j_u, j_d, j_tot, m_u, m_d, m_tot) == 0) {continue;}
    if (done_one == 0) {
      wd_f = u_tensor_u(d_tensor_d(wd, j_d, m_d), j_u, m_u);
      wfn_mult(wd_f, clebsch_gordan(j_u, j_d, j_tot, m_u, m_d, m_tot));
      done_one = 1;
    } else {
      wfn_sum(wd_f, u_tensor_u(d_tensor_d(wd, j_d, m_d), j_u, m_u), clebsch_gordan(j_u, j_d, j_tot, m_u, m_d, m_tot));
    }
  }
  return wd_f;
}

wfnData* u_tensor_d(wfnData* wd, double j_tot, double m_tot) {
  int n = wd->n_spin_up;
  if (j_tot > n) {printf("Impossibly large total angular momentum J = %d, J_max = %d\n", j_tot, n); exit(0);}
  if (fabs(m_tot) > j_tot) {printf("Magnetic quantum number exceeds total angular momentum\n"); exit(0);}
  wfnData* wd_f;
  int done_one = 0;
  for (int im1 = 0; im1 < n + 1; im1++) {
    double m1 = im1 - n/2.0;
    double m2 = m_tot - m1;
    if (fabs(m2) > n/2.0) {continue;}
    if (clebsch_gordan(n/2.0, n/2.0, j_tot, m1, m2, m_tot) == 0) {continue;}
    if (done_one == 0) {
      wd_f = shift_op(shift_op(wd, m2, 1), m1, 0);
      wfn_mult(wd_f, clebsch_gordan(n/2.0, n/2.0, j_tot, m1, m2, m_tot));
      done_one = 1;
    } else {
      wfn_sum(wd_f, shift_op(shift_op(wd, m2, 1), m1, 0), clebsch_gordan(n/2.0, n/2.0, j_tot, m1, m2, m_tot));
    }
  }
  return wd_f;
}

wfnData* u_tensor_u(wfnData* wd, double j_tot, double m_tot) {
  int n = wd->n_spin_up;
  if (j_tot > n) {printf("Impossibly large total angular momentum J = %d, J_max = %d\n", j_tot, n); exit(0);}
  if (fabs(m_tot) > j_tot) {printf("Magnetic quantum number exceeds total angular momentum\n"); exit(0);}
  wfnData* wd_f;
  int done_one = 0;
  for (int im1 = 0; im1 < n + 1; im1++) {
    double m1 = im1 - n/2.0;
    double m2 = m_tot - m1;
    if (fabs(m2) > n/2.0) {continue;}
    if (clebsch_gordan(n/2.0, n/2.0, j_tot, m1, m2, m_tot) == 0) {continue;}
    if (done_one == 0) {
      wd_f = shift_op(shift_op(wd, m2, 0), m1, 0);
      wfn_mult(wd_f, clebsch_gordan(n/2.0, n/2.0, j_tot, m1, m2, m_tot));
      done_one = 1;
    } else {
      wfn_sum(wd_f, shift_op(shift_op(wd, m2, 0), m1, 0), clebsch_gordan(n/2.0, n/2.0, j_tot, m1, m2, m_tot));
    }
  }
  return wd_f;
}


wfnData* d_tensor_d(wfnData* wd, double j_tot, double m_tot) {
  int n = wd->n_spin_up;
  if (j_tot > n) {printf("Impossibly large total angular momentum J = %d, J_max = %d\n", j_tot, n); exit(0);}
  if (fabs(m_tot) > j_tot) {printf("Magnetic quantum number exceeds total angular momentum\n"); exit(0);}
  wfnData* wd_f;
  int done_one = 0;
  for (int im1 = 0; im1 < n + 1; im1++) {
    double m1 = im1 - n/2.0;
    double m2 = m_tot - m1;
    if (fabs(m2) > n/2.0) {continue;}
    if (clebsch_gordan(n/2.0, n/2.0, j_tot, m1, m2, m_tot) == 0) {continue;}
    if (done_one == 0) {
      wd_f = shift_op(shift_op(wd, m2, 1), m1, 1);
      wfn_mult(wd_f, clebsch_gordan(n/2.0, n/2.0, j_tot, m1, m2, m_tot));
      done_one = 1;
    } else {
      wfn_sum(wd_f, shift_op(shift_op(wd, m2, 1), m1, 1), clebsch_gordan(n/2.0, n/2.0, j_tot, m1, m2, m_tot));
    }
  }
  return wd_f;
}


void wfn_sum(wfnData* wd1, wfnData* wd2, double fact) {

  if (wd1->n_states != wd2->n_states) {printf("Incompatible wfns %d %d\n", wd1->n_states, wd2->n_states); exit(0);}
  for (int i_state = 0; i_state < wd1->n_states; i_state++) {
    wd1->bc[i_state] += wd2->bc[i_state]*fact;
  }
  return;
}

double wfn_dot(wfnData* wd1, wfnData* wd2) {
  double dot = 0.0;
  if (wd1->n_states != wd2->n_states) {printf("Incompatible wfns %d %d\n", wd1->n_states, wd2->n_states); exit(0);}
  for (int i_state = 0; i_state < wd1->n_states; i_state++) {
    dot += wd1->bc[i_state]*wd2->bc[i_state];
  }
  
  return dot;
}

void wfn_mult(wfnData* wd, double fact) {
  for (int i_state = 0; i_state < wd->n_states; i_state++) {
    wd->bc[i_state] *= fact;
  }
  return;
}


wfnData* order_wfn(wfnData* wd) {
  // Orders the wfn basis according to my prescription
  clock_t start, end;
  int n = wd->n_spin_up;
  int b_dim_f = 2*wd->b_mag + 1;
  int* compress = (int*) calloc(n, sizeof(int));
  int *m_array = (int*) malloc(sizeof(int)*n);
  wfnData* wd_f = malloc(sizeof(*wd_f)); 
  wd_f->b_mag = wd->b_mag;
  wd_f->n_spin_up = n;
  wd_f->n_states = gsl_sf_choose(b_dim_f, n);
  wd_f->basis = (slater_det**) calloc(wd_f->n_states, sizeof(slater_det*));
  wd_f->bc = (double*) calloc(wd_f->n_states, sizeof(double));
  double* coeff = (double*) calloc(pow(b_dim_f, n), sizeof(double));
  
  for (unsigned int i_state = 0; i_state < wd->n_states; i_state++) {
    unsigned int index = 0;
    for (int i = 0; i < n; i++) {
      m_array[i] = wd->basis[i_state]->m_val[i] + wd->b_mag;
    }
    perm_compress(m_array, &compress, n);
    int phase = perm_sign(compress, n);
    order_perm(&m_array, n);
    for (int i = 0; i < n; i++) {
      index += m_array[i];
      if (i != n-1) {index *= b_dim_f;}
    }
    coeff[index] = wd->bc[i_state]*phase;
  }
  unsigned int m_num_max = pow(b_dim_f, n); 
  unsigned int i_state = 0;
  for (unsigned int im = 0; im < m_num_max; im++) {
    unsigned int m_store = im;
    int skip = 0;
    unsigned int index = 0;
    for (int i = 0; i < n; i++) {
      m_array[i] = (m_store % b_dim_f);
      if (i > 0) {if (m_array[i] <= m_array[i - 1]) {skip = 1; break;}}
      m_store -= m_array[i];
      m_store /= b_dim_f;
      index += m_array[i];
      if (i != n-1) {index *= b_dim_f;}
    }
    // Enforce anti-symmetry
    if (skip) {continue;}
    wd_f->basis[i_state] = malloc(sizeof(slater_det*));
    wd_f->basis[i_state]->m_val = malloc(sizeof(float)*n);
    for (int i = 0; i < n; i++) {
      wd_f->basis[i_state]->m_val[i] = m_array[i] - wd_f->b_mag;
    }
    wd_f->bc[i_state] = coeff[index];
    i_state++;
  }

  return wd_f;
}   

wfnData* shift_op(wfnData* wd, double m, int i_type) {
  // Acts the shift operator on the given FQHE wfn
  // The standard shift operator raises the magnetic field by one unit b -> b + 1/2
  clock_t start, end;
  double cpu_time;
  start = clock();
  int n = wd->n_spin_up;
  int r_up = m + n/2.0;
  int ph_sym = 0;
  if (r_up == 0) {r_up = n; ph_sym = 1;}
  unsigned int perm_min = pow(2, r_up) - 1;
  unsigned int perm_max = pow(2, n) - pow(2, n - r_up);
  // B field is raised by one unit
  double b_initial = wd->b_mag;
  double b_final = wd->b_mag + 0.5;
  if (i_type == 1) {b_final -= 1;}
  int b_dim_f = (int) (2*b_final + 1);
  double* coeff = (double*) calloc(pow(b_dim_f, n), sizeof(double));
  double pre_fact = sqrt(gsl_sf_fact(n/2.0 - m)*gsl_sf_fact(n/2.0 + m)/gsl_sf_fact(n));
  for (unsigned int perm = perm_min; perm <= perm_max;) {
    for (int i_state = 0; i_state < wd->n_states; i_state++) {
      unsigned int index = 0;
      unsigned int p_store = perm;
      int bad_rep = 0;
      double cg_fact = pre_fact;
      if (wd->bc[i_state] == 0.0) {continue;}
      for (int i_part = 0; i_part < n; i_part++) {
        double m1 = (p_store % 2) - 0.5;
        if (ph_sym) {m1 *= -1;}
        double m2 = wd->basis[i_state]->m_val[i_part];
        float mf = m1 + m2;
        if (fabs(mf) > b_final) {bad_rep = 1; break;}
        p_store -= (p_store % 2);
        p_store /= 2;
        index += mf + b_final;
        if (i_part != n - 1) {index *= b_dim_f;}
        if (i_type == 0) {
          cg_fact *= clebsch_gordan(0.5, b_initial, b_final, m1, m2, mf);
        } else if (i_type == 1 && m1 == 0.5) {
          cg_fact *= -sqrt(2*b_initial*(b_initial - m2));
        } else {
          cg_fact *= sqrt(2*b_initial*(b_initial + m2));
        }
//        cg_fact *= clebsch_gordan(0.5, b_initial, b_final, m1, m2, mf);
      }
      if (bad_rep) {continue;}
      coeff[index] += cg_fact*wd->bc[i_state];
    }
    perm = next_perm(perm);
  }
  double norm = 0.0;
  unsigned int m_num_max = pow(b_dim_f, n);
  int* compress = (int*) calloc(n, sizeof(int));
  int* bseen = (int*) calloc(b_dim_f, sizeof(int));
  int *m_array = (int*) malloc(sizeof(int)*n);
  wfnData* wd_f = malloc(sizeof(*wd_f)); 
  wd_f->b_mag = b_final;
  wd_f->n_spin_up = n;
  for (unsigned int im = 0; im < m_num_max; im++) {
    unsigned int m_store = im;
    int i_perm = 0;
    int repeat = 0;
    unsigned int index = 0;
    for (int i = 0; i < b_dim_f; i++) {bseen[i] = 0;}
    for (int i = 0; i < n; i++) {
      m_array[i] = (m_store % b_dim_f);
      if (bseen[m_array[i]]) {repeat = 1; break;}
      if (i > 0) {if (m_array[i] <= m_array[i - 1]) {i_perm = 1;}}
      m_store -= m_array[i];
      m_store /= b_dim_f;
      bseen[m_array[i]] = 1;
      index += m_array[i];
      if (i != n-1) {index *= b_dim_f;}
    }
    if (repeat) {continue;}
    if (coeff[index] == 0.0) {continue;}
    if (i_perm) {
      perm_compress(m_array, &compress, n);
      int phase = perm_sign(compress, n);
      order_perm(&m_array, n);
      int indexp = 0;
      for (int i = 0; i < n; i++) {
        indexp += m_array[i];
        if (i != n-1) {indexp *= b_dim_f;}
      }
      coeff[indexp] += coeff[index]*phase;
    }
  }
  wd_f->n_states = gsl_sf_choose(b_dim_f, n);
  for (unsigned int im = 0; im < m_num_max; im++) {
    unsigned int m_store = im;
    int skip = 0;
    unsigned int index = 0;
    for (int i = 0; i < n; i++) {
      m_array[i] = (m_store % b_dim_f);
      if (i > 0) {if (m_array[i] <= m_array[i - 1]) {skip = 1;}}
      m_store -= m_array[i];
      m_store /= b_dim_f;
      bseen[m_array[i]] = 1;
      index += m_array[i];
      if (i != n-1) {index *= b_dim_f;}
    }
    if (skip) {continue;}
  //  if (fabs(coeff[index]) > pow(10, -8)) {(wd_f->n_states)++;}
  }

  wd_f->basis = (slater_det**) calloc(wd_f->n_states, sizeof(slater_det*));
  wd_f->bc = (double*) calloc(wd_f->n_states, sizeof(double));
  unsigned int i_state = 0;
  for (unsigned int im = 0; im < m_num_max; im++) {
    unsigned int m_store = im;
    int skip = 0;
    unsigned int index = 0;
    for (int i = 0; i < n; i++) {
      m_array[i] = (m_store % b_dim_f);
      if (i > 0) {if (m_array[i] <= m_array[i - 1]) {skip = 1; break;}}
      m_store -= m_array[i];
      m_store /= b_dim_f;
      index += m_array[i];
      if (i != n-1) {index *= b_dim_f;}
    }
    // Enforce anti-symmetry
    if (skip) {continue;}
//    if (fabs(coeff[index]) < pow(10, -8)) {continue;}
    wd_f->basis[i_state] = malloc(sizeof(slater_det*));
    wd_f->basis[i_state]->m_val = malloc(sizeof(float)*n);
    for (int i = 0; i < n; i++) {
      wd_f->basis[i_state]->m_val[i] = m_array[i] - wd_f->b_mag;
    }
    wd_f->bc[i_state] = coeff[index];
    i_state++;
  }
  return wd_f;
}      

void normalize_wfn(wfnData* wd) {
  double norm = 0.0;
  for (unsigned int i_state = 0; i_state < wd->n_states; i_state++) {
    norm += pow(wd->bc[i_state], 2);
  }
  norm = sqrt(norm);
  printf("Norm: %g\n", norm);
  if (norm == 0) {printf("Error: state has zero norm.\n"); exit(0);}
  for (unsigned int i_state = 0; i_state < wd->n_states; i_state++) {
    wd->bc[i_state] /= norm;
  }
  return;
}

void print_wfn(wfnData* wd) {
  char indices[1000];
  char form[100];
  int n = wd->n_spin_up;
  for (unsigned int i_state = 0; i_state < wd->n_states; i_state++) {
    if (fabs(wd->bc[i_state]) < pow(10, -8)) {continue;}
    strcpy(indices, "");
    for (int i = 0; i < n; i++) {
      strcpy(form, "");
      sprintf(form, "%g ", wd->basis[i_state]->m_val[i]);
      strcat(indices, form);
    }
    printf("%s %g\n", indices, wd->bc[i_state]);
  //   printf("%g\n", wd->bc[i_state]);

  }
  return;
}

wfnData* hierarchy(double l, double s, int n) {
  clock_t start, end;
  double cpu_time;
  start = clock();
  double b = n - 1 + l - s;
  int bdim = 2*b + 1;
  printf("Magnetic field b = %g\n", b);
  int *k_array = (int*) malloc(sizeof(int)*n);
  double *m_array = (double*) malloc(sizeof(double)*n);
  double *q_array = (double*) malloc(sizeof(double)*n);
  double *m_final = (double*) calloc(n, sizeof(double));
  double *me_array = (double*) calloc(n, sizeof(double));
  double *coeff = (double*) calloc(pow(bdim, n), sizeof(double));
  unsigned int* perm_array = (unsigned int*) malloc(sizeof(unsigned int)*n);
  unsigned int* perm_min_array = (unsigned int*) malloc(sizeof(unsigned int)*n);
  unsigned int* perm_max_array = (unsigned int*) malloc(sizeof(unsigned int)*n);
  // Populate m and q arrays
  int i_hold = 0;
  int* bseen = (int*) calloc(bdim, sizeof(int));
  // Setup m and q values
  for (int i = 0; i < (int) 2*l + 1; i++) {
    double l_i = -l + i;
    for (int j = 0; j < (int) 2*s + 1; j++) {
      double s_j = -s + j;
      m_array[i_hold] = l_i;
      q_array[i_hold] = s_j;
      i_hold++;
    }
  }
  // Choose one permutation of m and q values
  for (int i = 0; i < n; i++) {
    k_array[i] = i;
  }
  // Don't really need this prefactor because we normalize anyway
  double pre_fact = sqrt((2*s + 1)*(n - 2*s + 2* l)*(2*n + 2*l - 2*s - 1));
  pre_fact *= gsl_sf_gamma(n - 2*s)*sqrt(gsl_sf_gamma(2*s + 1.0)*gsl_sf_gamma(2*l + 1.0)*gsl_sf_gamma(n + 2*l - 2*s));
  pre_fact /= sqrt(gsl_sf_gamma(n + 1.0)*gsl_sf_gamma(n - 2*s + 2*l + 1)*gsl_sf_gamma(2*n + 2*l - 2*s)*gsl_sf_gamma(n)); 
//  double pre_fact = 1.0;
  // Loop over all vortex coupling values

  int me_dim = n - 2*s;
  unsigned int ime_max = pow(me_dim, n);

  printf("IME MAX: %u\n", ime_max);
  // ime is an integer that encodes the m_e(i) values for all i,
  // m_e(i) is the magnetic quantum number of the electron spinor [u_i]^{(N-1)/2 - s}
  for (unsigned int ime = 0; ime < ime_max; ime++) {
    if (ime % 10000 == 0) {printf("%g %\n", 100.0*((double) ime)/((double) ime_max));}
    unsigned int ime_store = ime;
    // Some realizations of m_e(i) are not physical, bad_rep keeps track of this
    int bad_rep = 0;
    // Loop over each electron
    double me_tot = 0.0;
    for (int i = 0; i < n; i++) {
      // Set the value of m_e(i) by modding ime by the appropriate factor
      me_array[i] = ime_store % me_dim;
      ime_store -= me_array[i];
      ime_store /= me_dim;
      // Translate between number of spin-up electrons and physical m_e values
      me_array[i] -= ((n - 1)/2.0 - s);
      // r_up gives the number of spin-up electrons in the vortex spinor of particle i
      int r_up = q_array[k_array[i]] - me_array[i] + (n-1)/2.0;
      if (r_up < 0) {bad_rep = 1; break;}
      perm_min_array[i] = pow(2, r_up) - 1;
      perm_array[i] = perm_min_array[i];
      perm_max_array[i] = pow(2, n - 1) - pow(2, n - 1 - r_up);
      m_final[i] = m_array[k_array[i]] + me_array[i] -(n-1)/2.0;
      me_tot += me_array[i];
    }
    if (bad_rep) {continue;}
    int phase = pow(-1.0, 1.5*n - 3.0*pow(n, 2)/2.0 + 3.0*me_tot + n*s);
    for (int i = 0; i < n; i++) {
      int mv = perm_array[i];
      for (int j = 0; j < n - 1; j++) {
        if (mv % 2) {
          if (j >= i) {
            m_final[j + 1]++;
          } else {
            m_final[j]++;
          }
  
        }
        mv -= (mv % 2);
        mv /= 2;
      }
    }

    int done = 0;
    while (!done) {
      unsigned int index = 0;
      double cg_fact = pre_fact*phase;
      int anti_sym = 1;
      for (int i = 0; i < bdim; i++) {bseen[i] = 0;};
      for (int i = 0; i < n; i++) {
        double m_f = m_final[i];
        int imf = m_f + b;
        if (bseen[imf]) {anti_sym = 0; break;}  
        bseen[imf] = 1;
      }
      if (!anti_sym) {
      //  printf("Not anti-sym\n");
        increment(&perm_array, perm_min_array, perm_max_array, &m_final, n, 0, &done);
        continue;
      }
      for (int i = 0; i < n; i++) {
        double q = q_array[k_array[i]];
        double me = me_array[i];
        double ri = q - me;
        double m = m_array[k_array[i]];
        double m_f = m_final[i];
    //    printf("q: %g m: %g ri: %g si: %g me: %g m_final: %g\n", q, m, ri, si, me, m_final);
        index += m_f + b;
        index *= bdim;
        //cg_fact *= pow(-1.0, 1.5 - 3.0*n/2.0 - ri + q + s + 2*(m_f - me -m));
        cg_fact *= gsl_sf_fact(n/2.0 - 0.5 - ri)*gsl_sf_fact(n/2.0 - 0.5 + ri);
        cg_fact *= sqrt(gsl_sf_fact(b + m_f)*gsl_sf_fact(b - m_f)); 
        cg_fact *= 1.0/(gsl_sf_fact(n/2.0 - 0.5 + q - ri - s)*gsl_sf_fact(n/2.0 - 0.5 - q + ri - s));
        //cg_fact *= 1.0/sqrt(gsl_sf_gamma(s - q + 1)*gsl_sf_gamma(s + q + 1)*gsl_sf_gamma(l - m + 1)*gsl_sf_gamma(l + m + 1));
      }
      index /= bdim;
     // if (done) {exit(0);}
      increment(&perm_array, perm_min_array, perm_max_array, &m_final, n, 0, &done);
      coeff[index] += cg_fact;
    }
  }
//  for (int i = 0; i < pow(bdim, n); i++) {
//    if (fabs(coeff[i]) > pow(10, -10)) {printf("%g\n", coeff[i]);}
//  }
  double norm = 0.0;
  unsigned int m_num_max = pow(bdim, n);
  int* compress = (int*) malloc(sizeof(int)*n);
  wfnData* wd_f = malloc(sizeof(*wd_f)); 
  wd_f->b_mag = b;
  wd_f->n_spin_up = n;
  printf("m_num: %u\n", m_num_max);
  for (unsigned __int128 ik = 0; ik < m_num_max; ik++) {
    unsigned __int128 k_store = ik;
    int perm = 0;
    int repeat = 0;
    unsigned __int128 index = 0;
    for (int i = 0; i < bdim; i++) {bseen[i] = 0;}

    for (int i = 0; i < n; i++) {
      k_array[i] = (k_store % bdim);
      if (bseen[k_array[i]]) {repeat = 1; break;}
      if (i > 0) {if (k_array[i] <= k_array[i - 1]) {perm = 1;}}
      k_store -= k_array[i];
      k_store /= bdim;
      bseen[k_array[i]] = 1;
      index += k_array[i];
      if (i != n-1) {index *= bdim;}
    }
    if (repeat) {continue;}
    if (coeff[index] == 0.0) {continue;}
    if (perm) {
      perm_compress(k_array, &compress, n);
      int phase = perm_sign(compress, n);
      order_perm(&k_array, n);
      unsigned __int128 indexp = 0;
      for (int i = 0; i < n; i++) {
        indexp += k_array[i];
        if (i != n-1) {indexp *= bdim;}
      }
      coeff[indexp] += coeff[index]*phase;
    }
  }
  wd_f->n_states = gsl_sf_choose(bdim, n);
  wd_f->basis = (slater_det**) calloc(wd_f->n_states, sizeof(slater_det*));
  wd_f->bc = (double*) calloc(wd_f->n_states, sizeof(double));
  unsigned int i_state = 0;
  for (unsigned int im = 0; im < m_num_max; im++) {
    unsigned int m_store = im;
    int skip = 0;
    unsigned int index = 0;
    for (int i = 0; i < n; i++) {
      k_array[i] = (m_store % bdim);
      if (i > 0) {if (k_array[i] <= k_array[i - 1]) {skip = 1; break;}}
      m_store -= k_array[i];
      m_store /= bdim;
      index += k_array[i];
      if (i != n-1) {index *= bdim;}
    }
    // Enforce anti-symmetry
    if (skip) {continue;}
//    if (fabs(coeff[index]) < pow(10, -8)) {continue;}
    wd_f->basis[i_state] = malloc(sizeof(slater_det*));
    wd_f->basis[i_state]->m_val = malloc(sizeof(float)*n);
    for (int i = 0; i < n; i++) {
      wd_f->basis[i_state]->m_val[i] = k_array[i] - wd_f->b_mag;
    }
    wd_f->bc[i_state] = coeff[index];
    i_state++;
  }
  return wd_f;
}   
 
 
void compute_1body_energy(int np, double m, double e_shift) {
  
  double coeff = (-1)*np/m;
  double val1 = (-1)*e_shift/np;
  
  for (int i = 0; i < (int) 2*m + 1; i++) {
    //double mi = - m + i;
    double val2 = 0;
    for (int j = 0; j < (int) 2*m + 1; j++) {
      //double mj = -m + j;
      val2 += xme(i, j, i, j);
    }
    val2 = val2*coeff + val1;
    printf("%g\n", val2);
  }

  return;
}

double xme(int m1p, int m2p, int m1, int m2) {

  double* factlog = (double*) malloc(201*sizeof(double));
  double* xint = (double*) malloc(201*sizeof(double));

  factlog[0] = 0.0;
  factlog[1] = 0.0;
  for (int i = 2; i < 201; i++) {
    double x = i;
    factlog[i] = factlog[i - 1] + log(x);
  }
  xint[0] = log(0.886226925);
  for (int i = 1; i < 201; i++) {
    double x = i;
    xint[i] = xint[i-1] + log(x - 0.5);
  }
  double factor = exp(-(m1 + m2)*log(2.0)+0.5*factlog[m1p]+0.5*factlog[m2p] + 0.5*factlog[m1] + 0.5*factlog[m2]);

  int i1p = m1p + 1;
  int i2p = m2p + 1;
  int i1 = m1 + 1;
  int i2 = m2 + 1;
  double val = 0.0;
  for (int j1p = 0; j1p < i1p; j1p++) {
    for (int j2p = 0; j2p < i2p; j2p++) {
      for (int j1 = 0; j1 < i1; j1++) {
        for (int j2 = 0; j2 < i2; j2++) {
          if (i1 + i2 - i1p - i2p != 0) {continue;}
//          printf("%d, %d, %d, %d\n", j1p, j2p, j1, j2);
          val += pow(-1.0, j2 + j2p)*exp(factlog[m1 + m2 - j1 - j2] - factlog[m1p - j1p] - factlog[j1p] -factlog[m2p - j2p] - factlog[j2p] - factlog[m1 - j1] - factlog[j1] -factlog[m2 - j2] - factlog[j2] + xint[j1 + j2]);
        }
      }
    }
  }
  double xme = val*factor;

  return xme;
}

void generate_interaction_file(double q, int ll) {
  double l = q + ll; 
  for (int j = 0; j <= (int) 2*l; j++) {
    double h = 0;
    if ((int) ( 2*l + j) % 2 == 0) {continue;}
    for (int n = 0; n <= (int) 2*l; n++) {
      h += pow(-1.0, 2*q + j)/sqrt(q)*pow(2*l + 1, 2)*six_j(j, l, l, n, l, l)*pow(three_j(l, n, l, -q, 0, q), 2.0);
    }
    
//    printf("%d %d %d %d %d %d %g\n", 1, 1, 1, 1, j, 0, );
    printf("%d %d %d %d %d %d %g\n", 1, 1, 1, 1, j, 1, h);
  }

  return;
}

void generate_multilevel_interaction_file(double m) {
  for (int ij1 = 0; ij1 <= 1; ij1++) {
    double j1 = m + ij1;
    for (int ij2 = 0; ij2 <= 1; ij2++) {
      double j2 = m + ij2;
      for (int ij3 = 0; ij3 <= 1; ij3++) {
        double j3 = m + ij3;
        for (int ij4 = 0; ij4 <= 1; ij4++) {
          double j4 = m + ij4;
           
          for (int j = 0; j <= (int) 2*(m+1); j++) {
            if (j > j1 + j2 || j > j3 + j4 || j < fabs(j1 - j2) || j < fabs(j3 - j4)) {continue;}
            double h = 0.0;
            for (int l = 0; l <= (int) 2*(m+1); l++) {
              h += pow(-1.0, 2*m + j + 2*j3 + j4 + j2)*sqrt((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))*three_j(j1, l, j3, -m, 0, m)*three_j(j2, l, j4, -m, 0, m)*six_j(j, j2, j1, l, j3, j4);
              h += pow(-1.0, 1.0 + j3 + j4 - j + 2*m + j + 2*j4 + j3 + j2)*sqrt((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))*three_j(j1, l, j4, -m, 0, m)*three_j(j2, l, j3, -m, 0, m)*six_j(j, j2, j1, l, j4, j3);
            }
            h *= 1/sqrt(m);
            h *= 0.5; 
            if (j1 != j2) {h *= sqrt(2.0);}
            if (j3 != j4) {h *= sqrt(2.0);}
            if (fabs(h) < pow(10, -8)) {continue;}
            printf("%d %d %d %d %d %d %g\n", 1 + ij1, 1 + ij2, 1 + ij3, 1 + ij4, j, 1, h);
          }
        }
      }
    }
  }

  return;
}

