#include "hilbert.h"


void hierarchy_mult(double l, double s, int n) {
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
      double cg_fact = phase;
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
     //   cg_fact *= pow(-1.0, 1.5 - 3.0*n/2.0 - ri + q + s + 2*(m + me -m_f));
        cg_fact *= gsl_sf_fact(n/2.0 - 0.5 - ri)*gsl_sf_fact(n/2.0 - 0.5 + ri);
        cg_fact *= sqrt(gsl_sf_fact(b + m_f)*gsl_sf_fact(b - m_f)); 
        cg_fact *= 1.0/(gsl_sf_fact(n/2.0 - 0.5 + q - ri - s)*gsl_sf_fact(n/2.0 - 0.5 - q + ri - s));
      //  cg_fact *= 1.0/sqrt(gsl_sf_gamma(s - q + 1)*gsl_sf_gamma(s + q + 1)*gsl_sf_gamma(l - m + 1)*gsl_sf_gamma(l + m + 1));
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
  unsigned int k_num_max = pow(bdim, n);
  int* compress = (int*) malloc(sizeof(int)*n);
  for (unsigned __int128 ik = 0; ik < k_num_max; ik++) {
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
/*      printf("Before\n");
      for (int h = 0; h < n; h++) {
        printf("%d\n", k_array[h]);
      }
*/
      perm_compress(k_array, &compress, n);
      int phase = perm_sign(compress, n);
      order_perm(&k_array, n);
/*
      printf("After phase = %d\n", phase);

      for (int h = 0; h < n; h++) {
        printf("%d\n", k_array[h]);
      }
*/
      int indexp = 0;
      for (int i = 0; i < n; i++) {
        indexp += k_array[i];
        if (i != n-1) {indexp *= bdim;}
      }
  //    printf("%d %d %g %g\n", (int) index, (int) indexp, coeff[index], coeff[indexp]);
      coeff[indexp] += coeff[index]*phase;
    }
  }
  for (unsigned __int128 ik = 0; ik < k_num_max; ik++) {
    unsigned __int128 k_store = ik;
    int skip = 0;
    unsigned __int128 index = 0;
    for (int i = 0; i < n; i++) {
      k_array[i] = (k_store % bdim);
      if (i > 0) {if (k_array[i] <= k_array[i - 1]) {skip = 1; break;}}
      k_store -= k_array[i];
      k_store /= bdim;
      index += k_array[i];
      if (i != n-1) {index *= bdim;}
    }
    // Enforce anti-symmetry
    if (skip) {continue;}
   //   printf("index: %d\n", index);

    if (fabs(coeff[index]) > pow(10, -10)) {norm += pow(coeff[index], 2);}
  }
  end = clock();
  norm = sqrt(norm);
  printf("Norm: %g\n", norm);
  char indices[1000];
  char form[100];
  for (unsigned __int128 ik = 0; ik < k_num_max; ik++) {
    unsigned __int128 k_store = ik;
    int skip = 0;
    unsigned __int128 index = 0;
    strcpy(indices, "");
    for (int i = 0; i < n; i++) {
      strcpy(form, "");
      k_array[i] = (k_store % bdim);
      if (i > 0) {if (k_array[i] <= k_array[i - 1]) {skip = 1; break;}}
      k_store -= k_array[i];
      k_store /= bdim;
      index += k_array[i];
      sprintf(form, "%g ", k_array[i] - b);
      strcat(indices, form);
      if (i != n-1) {index *= bdim;}
    }
    // Enforce anti-symmetry
    if (skip) {continue;}
   //   printf("index: %d\n", index);

    if (fabs(coeff[index]) > pow(10, -10)) {printf("%s %d %g\n", indices, (int) index, coeff[index]/norm);}
  }
  printf("Time: %g\n", (double) (end - start)/CLOCKS_PER_SEC);
  return;
}
 
unsigned int next_perm(unsigned int v) {
  if (v == 0) {return 0;}
  unsigned int t = v | (v - 1);
  unsigned int w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));
//  printf("Next %d %d\n", (int) v, (int) w);
  return w;
}

void increment(unsigned int **perm, unsigned int *perm_min, unsigned int *perm_max, double** m_final, int n, int i, int* done) {
  if ((*perm)[i] == perm_max[i]) {
    int pi = (*perm)[i];
    int pf = perm_min[i];
    (*perm)[i] = pf;
    if (i < n - 1) {
      for (int j = 0; j < n - 1; j++) {
        if (pi % 2) {
          if (j >= i) {
            (*m_final)[j + 1]--;
          } else {
            (*m_final)[j]--;
          }
        }
        if (pf % 2) {
          if (j >= i) {
            (*m_final)[j + 1]++;
          } else {
            (*m_final)[j]++;
          }
        }
        pi -= (pi % 2);
        pi /= 2;
        pf -= (pf % 2);
        pf /= 2;
      } 

      increment(perm, perm_min, perm_max, m_final, n, i + 1, done);
    } else {
      *done = 1;
      return;
    }
  } else {
    int pi = (*perm)[i];
    int pf = next_perm(pi);
    (*perm)[i] = pf;
    for (int j = 0; j < n - 1; j++) {
      if (pi % 2) {
        if (j >= i) {
          (*m_final)[j + 1]--;
        } else {
          (*m_final)[j]--;
        }
      }
      if (pf % 2) {
        if (j >= i) {
          (*m_final)[j + 1]++;
        } else {
          (*m_final)[j]++;
        }
      }
      pi -= (pi % 2);
      pi /= 2;
      pf -= (pf % 2);
      pf /= 2;
    } 
  }
  //exit(0);
  return;
}


void hierarchy(double l, double s, int n) {
  clock_t start, end;
  double cpu_time;
  start = clock();
  double b = n - 1 + l - s;
  int bdim = 2*b + 1;
  printf("Magnetic field b = %g\n", b);
  int *k_array = (int*) malloc(sizeof(int)*n);
  unsigned __int128 k_num_max = pow(n, n);
  double *m_array = (double*) malloc(sizeof(double)*n);
  double *q_array = (double*) malloc(sizeof(double)*n);
  double *s_ar = (double*) calloc(n, sizeof(double));
  double *r_ar = (double*) calloc(n, sizeof(double));
  double *me_array = (double*) calloc(n, sizeof(double));
  int *i_couple = (int*) calloc(n, sizeof(int));
  int i_hold = 0;
  double *m_vortex = (double*) malloc(sizeof(double)*n*(n-1));
  unsigned __int128 mv_num_max = pow(2, n*(n-1));
  double *coeff = (double*) calloc(pow(bdim, n), sizeof(double));
  int* seen = (int*) calloc(n, sizeof(int));
  int* bseen = (int*) calloc(bdim, sizeof(int));
  printf("%u, %u\n", (int) k_num_max, (int) mv_num_max);
  rs_list** rs_mult = card_u_fqhe(n);
  // Populate m and q arrays
  for (int i = 0; i < (int) 2*l + 1; i++) {
    double l_i = -l + i;
    for (int j = 0; j < (int) 2*s + 1; j++) {
      double s_j = -s + j;
      m_array[i_hold] = l_i;
      q_array[i_hold] = s_j;
      i_hold++;
    }
  }
  for (int i = 0; i < n; i++) {
    k_array[i] = i;
  }
  double pre_fact = sqrt((2*s + 1)*(n - 2*s + 2* l)*(2*n + 2*l - 2*s - 1));
  pre_fact *= gsl_sf_gamma(n - 2*s)*sqrt(gsl_sf_gamma(2*s + 1.0)*gsl_sf_gamma(2*l + 1.0)*gsl_sf_gamma(n + 2*l - 2*s));
  pre_fact /= sqrt(gsl_sf_gamma(n + 1.0)*gsl_sf_gamma(n - 2*s + 2*l + 1)*gsl_sf_gamma(2*n + 2*l - 2*s)*gsl_sf_gamma(n)); 
//  double pre_fact = 1.0;
  // Loop over all vortex coupling values
  for (unsigned __int128 imv = 0; imv < mv_num_max; imv++) {
   // if (imv % 10000 == 0) {printf("%u\n", (int) imv);}
    unsigned __int128 m_store = imv;
    int ivor = 0;
    for (int i = 0; i < n; i++) {
      r_ar[i] = 0;
      s_ar[i] = 0;
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n - 1; j++) {
        int mv = (m_store % 2);
        m_store -= mv;
        m_store /= 2;
        r_ar[i] += (mv - 0.5);
        if (j >= i) {
          s_ar[j + 1] += (mv - 0.5);
        } else {
          s_ar[j] += (mv - 0.5);
        }
      }
    }
    int bad_rep = 0;
    for (int i = 0; i < n; i++) {
      if (fabs(q_array[k_array[i]] - r_ar[i]) > (n - 1.0)/2.0 - s) {bad_rep = 1; break;}
    }
          // Enforce anti-symmetry
    if (bad_rep) {continue;}
    unsigned __int128 index = 0;
 /*     printf("\n");
      for (int i = 0; i < n; i++) {
        printf("%d %g %g\n", k_array[i], m_array[k_array[i]], q_array[k_array[i]]);
      }*/
    double cg_fact = pre_fact;
    for (int i = 0; i < n; i++) {
      double q = q_array[k_array[i]];
      double me = q - r_ar[i];
      double m = m_array[k_array[i]];
      double m_final = m + me + s_ar[i]; 
      index += m_final + b;
      if (i != n - 1) {index *= bdim;}
      cg_fact *= pow(-1.0, 1.5 - 3.0*n/2.0 + q - r_ar[i] + +s + 2*s_ar[i]);
      cg_fact *= gsl_sf_gamma(n/2.0 + 0.5 - r_ar[i])*gsl_sf_gamma(n/2.0 + 0.5 + r_ar[i]);
      cg_fact *= sqrt(gsl_sf_gamma(b + m_final + 1.0)*gsl_sf_gamma(b - m_final + 1.0)); 
      cg_fact *= 1.0/(gsl_sf_gamma(n/2.0 + 0.5 + q - r_ar[i] - s)*gsl_sf_gamma(n/2.0 + 0.5 - q + r_ar[i] - s));
      cg_fact *= 1.0/sqrt(gsl_sf_gamma(s - q + 1)*gsl_sf_gamma(s + 1 + 1)*gsl_sf_gamma(l - m + 1)*gsl_sf_gamma(l + m + 1));
      //if (cg_fact == 0.0) {break;}
/*
      cg_fact *= clebsch_gordan((n-1)/2.0 - s, (n-1)/2.0, s, me, r_ar[i], q_array[k_array[i]]);
      cg_fact *= sqrt(gsl_sf_gamma((n-1)/2.0 - r_ar[i] + 1.0)*gsl_sf_gamma((n-1)/2.0 + r_ar[i] + 1.0)/gsl_sf_gamma(n));
      cg_fact *= clebsch_gordan(l, (n-1)/2.0 - s, l + (n-1)/2.0 - s, m_array[k_array[i]], me, m_array[k_array[i]] + me);
      cg_fact *= clebsch_gordan(l + (n-1)/2.0 - s, (n-1)/2.0, l - s + n - 1, m_array[k_array[i]] + me, s_ar[i], m_final);
      cg_fact *= sqrt(gsl_sf_gamma((n-1)/2.0 - s_ar[i] + 1.0)*gsl_sf_gamma((n-1)/2.0 + s_ar[i] + 1.0)/gsl_sf_gamma(n));
*/
    }
    //if (cg_fact == 0.0) {continue;}
    //int phase = perm_sign(k_array, n);
    coeff[index] += cg_fact;//*phase;
  }
//  for (int i = 0; i < pow(bdim, n); i++) {
//    if (fabs(coeff[i]) > pow(10, -10)) {printf("%g\n", coeff[i]);}
//  }
  double norm = 0.0;
  k_num_max = pow(bdim, n);
  int* compress = (int*) malloc(sizeof(int)*n);
  for (unsigned __int128 ik = 0; ik < k_num_max; ik++) {
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
/*      printf("Before\n");
      for (int h = 0; h < n; h++) {
        printf("%d\n", k_array[h]);
      }
*/
      perm_compress(k_array, &compress, n);
      int phase = perm_sign(compress, n);
      order_perm(&k_array, n);
/*
      printf("After phase = %d\n", phase);

      for (int h = 0; h < n; h++) {
        printf("%d\n", k_array[h]);
      }
*/
      int indexp = 0;
      for (int i = 0; i < n; i++) {
        indexp += k_array[i];
        if (i != n-1) {indexp *= bdim;}
      }
  //    printf("%d %d %g %g\n", (int) index, (int) indexp, coeff[index], coeff[indexp]);
      coeff[indexp] += coeff[index]*phase;
    }
  }
  for (unsigned __int128 ik = 0; ik < k_num_max; ik++) {
    unsigned __int128 k_store = ik;
    int skip = 0;
    unsigned __int128 index = 0;
    for (int i = 0; i < n; i++) {
      k_array[i] = (k_store % bdim);
      if (i > 0) {if (k_array[i] <= k_array[i - 1]) {skip = 1; break;}}
      k_store -= k_array[i];
      k_store /= bdim;
      index += k_array[i];
      if (i != n-1) {index *= bdim;}
    }
    // Enforce anti-symmetry
    if (skip) {continue;}
   //   printf("index: %d\n", index);

    if (fabs(coeff[index]) > pow(10, -10)) {norm += pow(coeff[index], 2);}
  }
  end = clock();
  norm = sqrt(norm);
  printf("Norm: %g\n", norm);
  char indices[1000];
  char form[100];
  for (unsigned __int128 ik = 0; ik < k_num_max; ik++) {
    unsigned __int128 k_store = ik;
    int skip = 0;
    unsigned __int128 index = 0;
    strcpy(indices, "");
    for (int i = 0; i < n; i++) {
      strcpy(form, "");
      k_array[i] = (k_store % bdim);
      if (i > 0) {if (k_array[i] <= k_array[i - 1]) {skip = 1; break;}}
      k_store -= k_array[i];
      k_store /= bdim;
      index += k_array[i];
      sprintf(form, "%g ", k_array[i] - b);
      strcat(indices, form);
      if (i != n-1) {index *= bdim;}
    }
    // Enforce anti-symmetry
    if (skip) {continue;}
   //   printf("index: %d\n", index);

    if (fabs(coeff[index]) > pow(10, -10)) {printf("%s %d %g\n", indices, (int) index, coeff[index]/norm);}
  }
  printf("Time: %g\n", (double) (end - start)/CLOCKS_PER_SEC);
  return;
}

void generate_rs_matrix(double mtot1, double mtot2, double mtot3) {
  int rs_index = 0;
  for (int r1 = -1; r1 <= 1; r1++) {
    for (int r2 = -1; r2 <= 1; r2++) {
      for (int r3 = -1; r3 <= 1; r3++) {
        for (int s1 = -1; s1 <= 1; s1++) {
          for (int s2 = -1; s2 <= 1; s2++) {
            for (int s3 = -1; s3 <= 1; s3++) {
              if (r1 + r2 + r3 != s1 + s2 + s3) {continue;}
              double mat = 0.0; 
              //int rs_index = r1 + 1 + 3*(r2 + 1 + 3*(r3 + 1 + 3*(s1 + 1 + 3*(s2 + 1 + 3*(s3 + 1)))));
              rs_index++;
              for (int m1 = -1; m1 <= 1; m1++) {
                if (mtot1 != m1 - r1 + s1) {continue;}
                for (int m2 = -1; m2 <= 1; m2++) {
                  if (m2 == m1) {continue;}
                  if (mtot2 != m2 - r2 + s2) {continue;}
                  for (int m3 = -1; m3 <= 1; m3++) {
                    if ((m3 == m1) || (m3 == m2)) {continue;}
                    if (mtot3 != m3 - r3 + s3) {continue;}
                    double cg_fact = clebsch_gordan(1.0, 1.0, 0.0, -r1, r1, 0.0);
                    cg_fact *= clebsch_gordan(1.0, 1.0, 0.0, -r2, r2, 0.0);
                    cg_fact *= clebsch_gordan(1.0, 1.0, 0.0, -r3, r3, 0.0);
                    cg_fact *= sqrt(gsl_sf_gamma(2.0 - r1)*gsl_sf_gamma(2.0 + r1)/gsl_sf_gamma(3.0));
                    cg_fact *= sqrt(gsl_sf_gamma(2.0 - r2)*gsl_sf_gamma(2.0 + r2)/gsl_sf_gamma(3.0));
                    cg_fact *= sqrt(gsl_sf_gamma(2.0 - r3)*gsl_sf_gamma(2.0 + r3)/gsl_sf_gamma(3.0));
                    cg_fact *= clebsch_gordan(1.0, 1.0, 2.0, m1, -r1, m1 - r1);
                    cg_fact *= clebsch_gordan(1.0, 1.0, 2.0, m2, -r2, m2 - r2);
                    cg_fact *= clebsch_gordan(1.0, 1.0, 2.0, m3, -r3, m3 - r3);
                    cg_fact *= clebsch_gordan(2.0, 1.0, 3.0, m1 - r1, s1, mtot1);
                    cg_fact *= clebsch_gordan(2.0, 1.0, 3.0, m2 - r2, s2, mtot2);
                    cg_fact *= clebsch_gordan(2.0, 1.0, 3.0, m3 - r3, s3, mtot3);
                    cg_fact *= sqrt(gsl_sf_gamma(2.0 - s1)*gsl_sf_gamma(2.0 + s1)/gsl_sf_gamma(3.0));
                    cg_fact *= sqrt(gsl_sf_gamma(2.0 - s2)*gsl_sf_gamma(2.0 + s2)/gsl_sf_gamma(3.0));
                    cg_fact *= sqrt(gsl_sf_gamma(2.0 - s3)*gsl_sf_gamma(2.0 + s3)/gsl_sf_gamma(3.0));
                    if ((m1 == -1) && (m2 == 1) && (m3 == 0)) {cg_fact *= -1;}
                    if ((m1 == 1) && (m2 == 0) && (m3 == -1)) {cg_fact *= -1;}
                    if ((m1 == 0) && (m2 == -1) && (m3 == 1)) {cg_fact *= -1;}
                    mat += cg_fact;
                  } 
                }
              }
              printf("%g, ", mat);
            }
          }
        }
      }
    }
  }
  return;
}
int perm_sign(int* perm, int n) {
  int *seen = (int*) calloc(n, sizeof(int));
  int sign = 1;

  for (int k = 0; k < n; k++) {
    if (!seen[k]) {
      int next = k;
      int l = 0;
      while (!seen[next]) {
        l++;
        seen[next] = 1;
        next = perm[next];
      }
      if (l % 2 == 0) {
        sign *= -1;
      }
    }
  }
  return sign;
}

void order_perm(int** perm, int n) {

  int* seen = (int*) calloc(n, sizeof(int));
  int* ord = (int*) calloc(n, sizeof(int));

  for (int i = 0; i < n; i++) {
    int min = 1000;
    int j_min;
    for (int j = 0; j < n; j++) {
      if (seen[j]) {continue;}
      if ((*perm)[j] < min) {
        min = (*perm)[j];
        j_min = j;
      }
    }
    ord[i] = (*perm)[j_min];
    seen[j_min] = 1;
  }
  *perm = ord;
  return;
}

void perm_compress(int* perm, int** result, int n) { 
  int* seen = (int*) calloc(n, sizeof(int));

  for (int i = 0; i < n; i++) {
    int min = 1000;
    int j_min;
    for (int j = 0; j < n; j++) {
      if (seen[j]) {continue;}
      if (perm[j] < min) {
        min = perm[j];
        j_min = j;
      }
    }
    (*result)[j_min] = i;
    seen[j_min] = 1;
  }

  return;
}
 
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

void generate_interaction_file(double l) {

  for (int j = 0; j <= (int) 2*l; j++) {
    double h = 0;
    if ((int) ( 2*l + j) % 2 == 0) {continue;}
    for (int n = 0; n <= (int) 2*l; n++) {
      h += pow(-1.0, 2*l + j)/sqrt(l)*pow(2*l + 1, 2)*gsl_sf_coupling_6j ((int) 2*l, (int) 2*l, (int) 2*n, (int) 2*l, (int) 2*l, (int) 2*j)*pow(gsl_sf_coupling_3j((int) 2*l, (int) 2*n, (int) 2*l, (int) -2*l, 0, (int) 2*l), 2.0);
    }
    
//    printf("%d %d %d %d %d %d %g\n", 1, 1, 1, 1, j, 0, );
    printf("%d %d %d %d %d %d %g\n", 1, 1, 1, 1, j, 1, h);
  }

  return;
}

gsl_complex u_spin(double m, double theta, double phi) {
  gsl_complex u = gsl_complex_rect(0.0, 0.0);
  if (m == 0.5) {
    gsl_complex i_phi = gsl_complex_rect(0.0, phi/2.0);
    u = gsl_complex_mul_real(gsl_complex_exp(i_phi), cos(theta/2.0));
  } else if (m == -0.5) {
    gsl_complex i_phi = gsl_complex_rect(0.0, -phi/2.0);
    u = gsl_complex_mul_real(gsl_complex_exp(i_phi), sin(theta/2.0));
  } else {
    printf("Not a valid m-value %g\n", m); exit(0);
  }

  return u;
}

gsl_complex u_spin_l(double l, double m, double theta, double phi) {
  gsl_complex u = gsl_complex_rect(0.0, 0.0);
  if (l > 1.0) {
    for (int im = -1; im <= 1; im += 2) {
      double m1 = im/2.0;
      double m2 = m - m1;
      if ((m2 > l - 0.5) || (m2 < -1.0*(l - 0.5))) {continue;}
      double cg = pow(-1.0, l - 1.0 - m)*sqrt(2*l + 1)*gsl_sf_coupling_3j(1.0, (int) 2*(l - 0.5), (int) 2*l, (int) 2*m1, (int) 2*m2, (int) -2.0*m);
      if (cg == 0.0) {continue;}
      u = gsl_complex_add(u, gsl_complex_mul_real(gsl_complex_mul(u_spin(m1, theta, phi), u_spin_l(l - 0.5, m2, theta, phi)), cg));
    }
  } else if (l == 1.0) {
    for (int im = -1; im <= 1; im += 2) {
      double m1 = im/2.0;
      double m2 = m - m1;
      if ((m2 > 0.5) || (m2 < -0.5)) {continue;}
      double cg = pow(-1.0, l - 1.0 - m)*sqrt(2*l + 1)*gsl_sf_coupling_3j(1.0, 1.0, 2.0, (int) 2*m1, (int) 2*m2, (int) -2.0*m);
      if (cg == 0.0) {continue;}
      u = gsl_complex_add(u, gsl_complex_mul_real(gsl_complex_mul(u_spin(m1, theta, phi), u_spin(m2, theta, phi)), cg));
    }
  }
  return u;
}


void basis_coeff(int m, int n) {
  int dim = 2*m+1;
  double *basis = (double*) calloc(pow(dim, 3), sizeof(double));

  for (int k12 = 0; k12 <= m; k12++) {
    for (int k13 = 0; k13 <= m; k13++) {
      for (int k23 = 0; k23 <= m; k23++) {
        int im1 = k12 + k13 - m + m;
        int im2 = k23 - k12  + m;
        int im3 = 3 - k13 - k23 + m;
        basis[im1 + dim*(im2 + dim*im3)] += pow(-1.0, k12 + k13 + k23)*gsl_sf_choose(m, k12)*gsl_sf_choose(m, k13)*gsl_sf_choose(m, k23)*clebsch_gordan(1.5, 1.5, 3.0, k12 - 1.5, k13 - 1.5, k12 + k13 - 3)*clebsch_gordan(1.5, 1.5, 3.0, 1.5 - k12, k23 - 1.5, k23 - k12)*clebsch_gordan(1.5, 1.5, 3.0, 1.5 - k13, 1.5 - k23, 3 - k13 - k23)*gsl_sf_gamma(k12 + 1)*gsl_sf_gamma(k13 + 1)*gsl_sf_gamma(k23 + 1)*gsl_sf_gamma(3 - k12 + 1)*gsl_sf_gamma(3 - k13 + 1)*gsl_sf_gamma(3 - k23 + 1)/pow(gsl_sf_gamma(4.0), 3);
      }
    }
  }
  double norm = 0.0;
  for (int m1 = -3; m1 <= 3; m1++) {
    int im1 = m1 + 3;
    for (int m2 = -3; m2 <= m1; m2++) {
      int im2 = m2 + 3;
      for (int m3 = -3; m3 <= m2; m3++) {
        int im3 = m3 + 3;
        if (fabs(basis[im1 + dim*(im2 + dim*im3)]) > pow(10, -10)) {norm += pow(basis[im1 + dim*(im2 + dim*im3)], 2);}
      }
    }
  }
  norm = sqrt(norm);
  for (int m1 = -3; m1 <= 3; m1++) {
    int im1 = m1 + 3;
    for (int m2 = -3; m2 <= m1; m2++) {
      int im2 = m2 + 3;
      for (int m3 = -3; m3 <= m2; m3++) {
        int im3 = m3 + 3;
        if (fabs(basis[im1 + dim*(im2 + dim*im3)]) > pow(10, -10)) {printf("%d, %d, %d, %g\n", m1, m2, m3, basis[im1 + dim*(im2 + dim*im3)]/norm);}
      }
    }
  }

  return;
}

void excited_state(int m, int n) {
  int dim = 2*m+1;
  int dim2 = 2*m + 2;
  double *basis = (double*) calloc(pow(dim, 3), sizeof(double));
  double *excite = (double*) calloc(pow(dim2, 3), sizeof(double));
  for (int k12 = 0; k12 <= m; k12++) {
    for (int k13 = 0; k13 <= m; k13++) {
      for (int k23 = 0; k23 <= m; k23++) {
        int im1 = k12 + k13 - m + m;
        int im2 = k23 - k12  + m;
        int im3 = 3 - k13 - k23 + m;
        basis[im1 + dim*(im2 + dim*im3)] += pow(-1.0, k12 + k13 + k23)*clebsch_gordan(1.5, 1.5, 3.0, k12 - 1.5, k13 - 1.5, k12 + k13 - 3)*clebsch_gordan(1.5, 1.5, 3.0, 1.5 - k12, k23 - 1.5, k23 - k12)*clebsch_gordan(1.5, 1.5, 3.0, 1.5 - k13, 1.5 - k23, 3 - k13 - k23);
      }
    }
  }
  double norm = 0.0;
  for (int im1 = 0; im1 <= 6; im1++) {
    for (int im2 = 0; im2 <= im1; im2++) {
      for (int im3 = 0; im3 <= im2; im3++) {
        if (fabs(basis[im1 + dim*(im2 + dim*im3)]) > pow(10, -10)) {norm += pow(basis[im1 + dim*(im2 + dim*im3)], 2);}
      }
    }
  }
  norm = sqrt(norm); 
  for (int im1 = 0; im1 <= 6; im1++) {
    for (int im2 = 0; im2 <= 6; im2++) {
      for (int im3 = 0; im3 <= 6; im3++) {
        if (fabs(basis[im1 + dim*(im2 + dim*im3)]) > pow(10, -10)) {basis[im1 + dim*(im2 + dim*im3)] *= 1/norm;}
      }
    }
  }
  for (int im1 = 0; im1 <= 6; im1++) {
    double m1 = im1 - 3.0;
    for (int im2 = 0; im2 <= 6; im2++) {
      double m2 = im2 - 3.0;
      for (int im3 = 0; im3 <= 6; im3++) {
        double m3 = im3 - 3.0;
        for (int im1p = -1; im1p <= 1; im1p += 2) {
          double m1p = im1p/2.0;
          for (int im2p = -1; im2p <= 1; im2p += 2) {
            double m2p = im2p/2.0;
            for (int im3p = -1; im3p <= 1; im3p += 2) {
              double m3p = im3p/2.0;
              double fact = clebsch_gordan(0.5, 0.5, 1, m2p, m3p, m2p + m3p)*clebsch_gordan(0.5, 1.0, 1.5, m1p, m2p + m3p, 0.5);
              if (fact == 0.0) {continue;}
              int im1pp = m1 + m1p + 3.5;
              int im2pp = m2 + m2p + 3.5;
              int im3pp = m3 + m3p + 3.5;
              excite[im1pp + dim2*(im2pp + dim2*im3pp)] += basis[im1 + dim*(im2 + dim*im3)]*fact*clebsch_gordan(0.5, 3.0, 3.5, m1p, m1, m1p + m1)*clebsch_gordan(0.5, 3.0, 3.5, m2p, m2, m2p + m2)*clebsch_gordan(0.5, 3.0, 3.5, m3p, m3, m3p + m3);
            }
          }
        }
      }
    }
  }
  norm = 0.0;
  for (int im1 = 0; im1 <= 7; im1++) {
    double m1 = (im1 - 3.5);
    for (int im2 = 0; im2 <= im1; im2++) {
      double m2 = (im2 - 3.5);
      for (int im3 = 0; im3 <= im2; im3++) {
        double m3 = (im3 - 3.5);
        if (fabs(excite[im1 + dim2*(im2 + dim2*im3)]) > pow(10, -10)) {norm += pow(excite[im1 + dim2*(im2 + dim2*im3)], 2.0);}
      }
    }
  }
  norm = sqrt(norm);
  printf("%g\n",1.0/ norm);
  for (int im1 = 0; im1 <= 7; im1++) {
    double m1 = (im1 - 3.5);
    for (int im2 = 0; im2 <= im1; im2++) {
      double m2 = (im2 - 3.5);
      for (int im3 = 0; im3 <= im2; im3++) {
        double m3 = (im3 - 3.5);
        if (fabs(excite[im1 + dim2*(im2 + dim2*im3)]) > pow(10, -10)) {printf("%g, %g, %g, %g\n", m1, m2, m3, excite[im1 + dim2*(im2 + dim2*im3)]/norm);}
      }
    }
  }

  return;
}

void laughlin(int m, int n) {

  int num_k = gsl_sf_choose(n, 2);
  int m_dim = 2*m + 1;
  double m_tot = m/2.0*(n-1);
  int m_tot_dim = 2*m_tot + 1;
  int *k_array = (int*) malloc(sizeof(int)*num_k);
  double *m_array = (double*) calloc(n, sizeof(double));
  unsigned int k_num_max = pow(m + 1, num_k);
  double *basis = (double*) calloc(pow(m_tot_dim, n), sizeof(double));
  for (int ik = 0; ik < k_num_max; ik++) {
    int k_store = ik;
    for (int i = 0; i < num_k; i++) {
      k_array[i] = (k_store % (m + 1));
      k_store -= k_array[i];
      k_store /= (m + 1);
    }
    double cg_fact = 1.0;
    for (int i = 0; i < n; i++) {
      m_array[i] = 0.0;
    }
    for (int i = 1; i <= n; i++) {
      int i_couple = 0;
      for (int j = 1; j <= i - 1; j++) {
        double m_prev = m_array[i - 1];
        m_array[i - 1] += m/2.0 - k_array[(j - 1) + (i - 1)*(i - 2)/2];
        if (i_couple == 0) {i_couple++; continue;}
        cg_fact *= clebsch_gordan(i_couple*m/2.0, m/2.0, (i_couple + 1)*m/2.0, m_prev, m/2.0 - k_array[(j - 1) + (i - 1)*(i - 2)/2], m_array[i - 1]);
        i_couple++;
      }
      for (int j = i + 1; j <= n; j++ ) {
        double m_prev = m_array[i - 1];
        m_array[i - 1] += k_array[(i - 1) + (j - 1)*(j - 2)/2] - m/2.0;
        if (i_couple == 0) {i_couple++; continue;}
        cg_fact *= clebsch_gordan(i_couple*m/2.0, m/2.0, (i_couple + 1)*m/2.0, m_prev, k_array[(i - 1) + (j - 1)*(j - 2)/2] - m/2.0, m_array[i - 1]);
        i_couple++;
      }
    }
    if (cg_fact == 0.0) {continue;}
    int index = 0;
    for (int i = 1; i <= n - 1; i++) {
      index += m_array[i - 1] + m_tot;
      index *= m_tot_dim;
    }
    index += m_array[n - 1] + m_tot;
    int k_sum = 0;
    for (int i = 0; i < num_k; i++) {
      k_sum += k_array[i];
    }
    basis[index] += pow(-1.0, k_sum)*cg_fact;
  }
  double norm = 0.0;
  for (int im1 = 0; im1 < m_tot_dim; im1++) {
    for (int im2 = 0; im2 <= im1; im2++) {
      for (int im3 = 0; im3 <= im2; im3++) {
        for (int im4 = 0; im4 <= im3; im4++) {
          for (int im5 = 0; im5 <= im4; im5++) {
            if (fabs(basis[im1 + m_tot_dim*(im2 + m_tot_dim*(im3 + m_tot_dim*(im4 + m_tot_dim*im5)))]) > pow(10, -10)) {norm += pow(basis[im1 + m_tot_dim*(im2 + m_tot_dim*(im3 + m_tot_dim*(im4 + m_tot_dim*im5)))], 2.0);}
          }
        }
      }
    }
  }
  norm = sqrt(norm);
  for (int im1 = 0; im1 < m_tot_dim; im1++) {
    double m1 = (im1 - m_tot);
    for (int im2 = 0; im2 <= im1; im2++) {
      double m2 = (im2 - m_tot);
      for (int im3 = 0; im3 <= im2; im3++) {
        double m3 = (im3 - m_tot);
        for (int im4 = 0; im4 <= im3; im4++) {
          double m4 = (im4 - m_tot);
          for (int im5 = 0; im5 <= im4; im5++) {
            double m5 = (im5 - m_tot);
            if (fabs(basis[im1 + m_tot_dim*(im2 + m_tot_dim*(im3 + m_tot_dim*(im4 + m_tot_dim*im5)))]) > pow(10, -10)) {printf("%g, %g, %g, %g, %g, %g\n", m1, m2, m3, m4, m5, basis[im1 + m_tot_dim*(im2 + m_tot_dim*(im3 + m_tot_dim*(im4 + m_tot_dim*im5)))]/norm);}
          }
        }
      }
    }
  }
 

  return;
}

