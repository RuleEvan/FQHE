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

void hierarchy(double l, double s, int n) {
  double b = n - 1 + l - s;
  int bdim = 2*b + 1;
  printf("Magnetic field b = %g\n", b);
  int *k_array = (int*) malloc(sizeof(int)*n);
  unsigned int k_num_max = pow(n, n);
  double *m_array = (double*) malloc(sizeof(double)*n);
  double *q_array = (double*) malloc(sizeof(double)*n);
  double *m_tot_array = (double*) calloc(n, sizeof(double));
  double *me_array = (double*) calloc(n, sizeof(double));
  int *i_couple = (int*) calloc(n, sizeof(int));
  int i_hold = 0;
  double *m_vortex = (double*) malloc(sizeof(double)*n*(n-1));
  unsigned int mv_num_max = pow(2, n*(n-1));
  double *coeff = (double*) calloc(pow(bdim, n), sizeof(double));
  int* seen = (int*) calloc(n, sizeof(int));
  for (int i = 0; i < (int) 2*l + 1; i++) {
    double l_i = -l + i;
    for (int j = 0; j < (int) 2*s + 1; j++) {
      double s_j = -s + j;
      m_array[i_hold] = l_i;
      q_array[i_hold] = s_j;
      i_hold++;
    }
  }
  for (int ik = 0; ik < k_num_max; ik++) {
    int k_store = ik;
    int repeat = 0;
    for (int i = 0; i < n; i++) {seen[i] = 0;}
    printf("\n");
    for (int i = 0; i < n; i++) {
      k_array[i] = (k_store % n);
      if (seen[k_array[i]]) {repeat = 1; break;}
      seen[k_array[i]] = 1;
      k_store -= k_array[i];
      k_store /= n;
    }
    if (repeat) {continue;}
    int phase = perm_sign(k_array, n);
    for (int i = 0; i < n; i++) {
  //    printf("%d %g %g\n", k_array[i], m_array[k_array[i]], q_array[k_array[i]]);
    }
    for (unsigned int imv = 0; imv < mv_num_max; imv++) {
      unsigned int m_store = imv;
      for (int i = 0; i < n*(n-1); i++) {
        m_vortex[i] = (m_store % 2);
        m_store -= m_vortex[i];
        m_store /= 2;
        m_vortex[i] -= 0.5;
      }
      double cg_fact = 1.0;
      for (int i = 0; i < n; i++) {
        i_couple[i] = 0;
        m_tot_array[i] = 0.0;
      }
      int bad_rep = 0;
      for (int i = 0; i < n; i++) {
        double mv_tot = 0;
        for (int j = 0; j < n - 1; j++) {
          if (j > 0) {
            cg_fact *= clebsch_gordan(j/2.0, 0.5, (j + 1)/2.0, mv_tot, m_vortex[j + i*(n - 1)], mv_tot + m_vortex[j + i*(n - 1)]); 
          }
          if (j >= i) {
            double m_prev = m_tot_array[j + 1];
            m_tot_array[j + 1] += m_vortex[j + i*(n - 1)];
            if (i_couple[j + 1] != 0) {
              cg_fact *= clebsch_gordan(i_couple[j + 1]/2.0, 0.5, (i_couple[j + 1] + 1)/2.0, m_prev, m_vortex[j + i*(n - 1)], m_tot_array[j + 1]);
            }
            i_couple[j + 1]++; 
          } else {
            double m_prev = m_tot_array[j];
            m_tot_array[j] += m_vortex[j + i*(n-1)];
            if (i_couple[j] != 0) {
              cg_fact *= clebsch_gordan(i_couple[j]/2.0, 0.5, (i_couple[j] + 1)/2.0, m_prev, m_vortex[j + i*(n - 1)], m_tot_array[j]);
            }
            i_couple[j]++; 

          }
          mv_tot += m_vortex[j + i*(n - 1)];
        }
        me_array[i] = q_array[i] - mv_tot;
        if (fabs(me_array[i]) > (n-1)/2.0 -s) {bad_rep = 1; break;}
        cg_fact *= clebsch_gordan((n-1)/2.0 - s, (n-1)/2.0, s, me_array[i], mv_tot, q_array[i]);
      }
      if (bad_rep) {continue;}
      int index = 0;
      for (int i = 0; i < n; i++) {
        double m_final = m_array[i] + me_array[i] + m_tot_array[i]; 
        index += m_final + b;
        if (i != n - 1) {index *= bdim;}
        cg_fact *= clebsch_gordan(l, (n-1)/2.0 - s, l + (n-1)/2.0 - s, m_array[i], me_array[i], m_array[i] + me_array[i]);
        cg_fact *= clebsch_gordan(l + (n-1)/2.0 - s, (n-1)/2.0, l - s + n - 1, m_array[i] + me_array[i], m_tot_array[i], m_array[i] + me_array[i] + m_tot_array[i]);
        coeff[index] += cg_fact*phase;
      }
    }
  }
//  for (int i = 0; i < pow(bdim, n); i++) {
//    if (fabs(coeff[i]) > pow(10, -10)) {printf("%g\n", coeff[i]);}
//  }
  double norm = 0.0;
  for (int im1 = 0; im1 < bdim; im1++) {
    for (int im2 = 0; im2 <= im1; im2++) {
      for (int im3 = 0; im3 <= im2; im3++) {
        for (int im4 = 0; im4 <= im3; im4++) {
          if (fabs(coeff[im1 + bdim*(im2 + bdim*(im3 + bdim*im4))]) > pow(10, -10)) {norm += pow(coeff[im1 + bdim*(im2 + bdim*(im3 + bdim*im4))], 2.0);}
        }
      }
    }
  }
  norm = sqrt(norm);
  printf("Norm: %g\n", norm);
  for (int im1 = 0; im1 < bdim; im1++) {
    double m1 = (im1 - b);
    for (int im2 = 0; im2 <= im1; im2++) {
      double m2 = (im2 - b);
      for (int im3 = 0; im3 <= im2; im3++) {
        double m3 = (im3 - b);
        for (int im4 = 0; im4 <= im3; im4++) {
          double m4 = (im4 - b);
          if (fabs(coeff[im1 + bdim*(im2 + bdim*(im3 + bdim*im4))]) > pow(10, -10)) {printf("%g, %g, %g, %g, %g\n", m1, m2, m3, m4, coeff[im1 + bdim*(im2 + bdim*(im3 + bdim*im4))]);}
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
