#include "lanczos.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
double characteristic_polynomial(double lambda, double *alpha, double *beta, int dim, int* sign_diff) {
  double q1 = alpha[0] - lambda;
  double delta = pow(10, -16);
  *sign_diff = 0;
  if (q1 < 0) {*sign_diff = 1;}
  for (int k = 1; k < dim; k++) {
    double q2 = q1;
    if (q2 == 0.0) {q2 = delta;}
    q1 = (alpha[k] - lambda) - pow(beta[k-1], 2.0)/q2;
    if (q1 < 0) {*sign_diff += 1;}
  }
  return q1;
}
void lanczos_eigenvalues(double *h, int dim, int n_iter) {
  double tol = 0.01;
  double *alpha = (double*) malloc(sizeof(double)*n_iter);
  double *beta = (double*) malloc(sizeof(double)*(n_iter - 1));
  printf("Performing Lanczos tri-diagonalization...\n"); 
  lanczos(h, alpha, beta, dim, n_iter);
  printf("Done.\n");
  printf("Computing eigenvalues...\n");
  for (int k = 1; k <= n_iter; k++) {
    double b = alpha[0] + abs(beta[0]);
    for (int j = 1; j < n_iter; j++) {
      double val = 0.0;
      if (j == n_iter - 1) {
        val = alpha[j] + abs(beta[j-1]);
      } else {
        val = alpha[j] + abs(beta[j-1]) + abs(beta[j]);
      }
      if (val > b) {
        b = val;
      }
    } 
    double a = alpha[0] - abs(beta[0]);
    for (int j = 1; j < n_iter; j++) {
      double val = 0.0;
      if (j == n_iter - 1) {
        val = alpha[j] - abs(beta[j-1]);
      } else {
        val = alpha[j] - abs(beta[j-1]) - abs(beta[j]);
      }
      if (val < a) {
        a = val;
      }
    } 
        
    int k_min;
    int k_max;
    double lambda_mid = 0.0;
    double lambda_a = characteristic_polynomial(a, alpha, beta, n_iter, &k_min);
    double lambda_b = characteristic_polynomial(b, alpha, beta, n_iter, &k_max);
    while ((k != k_max) || (k != k_min + 1)) {
      double mid = (a + b)/2.0;
      int k_mid = 0.0;
      lambda_mid = characteristic_polynomial(mid, alpha, beta, n_iter, &k_mid);
      if (k_mid >= k) {
        k_max = k_mid;
        lambda_b = lambda_mid;
        b = mid;
      } else {
        k_min = k_mid;
        lambda_a = lambda_mid;
        a = mid;
      }
      if (k_min == k_max) {break;}
    }
    
    if (lambda_mid == 0) {
      printf("Eigenvalue exact: %g\n", b);
      continue;
    }
    while ((b - a) > tol) {
      double mid = (a+b)/2;
      int k_mid = 0.0;
      lambda_mid = characteristic_polynomial(mid, alpha, beta, n_iter, &k_mid);
      if (lambda_mid == 0.0) {
        a = mid; 
        b = mid;
      } else {
        if (k_mid == k-1) {
          a = mid; 
        } else {
          b = mid; 
        }
      }
    }
    printf("k: %d a:%g b:%g\n", k, a, b);
  }
  return;
}

void get_orthogonal_vector(double *h_tri, double *v, int dim, int n_done) {
  for (int i = 0; i < dim; i++) {
    v[i] = 0.0;
  }
  v[5] = 1.0;
  for (int j = 0; j < n_done; j++) {
    double proj = 0.0;
    for (int i = 0; i < dim; i++) {
      proj += v[i]*h_tri[j*dim + i];
    }
    for (int i = 0; i < dim; i++) {
      v[i] -= proj*h_tri[j*dim + i];
    }
  }
  return;
}

void square_matrix_times_vector(double* mat, double *vec, double *result, int dim) {
  for (int i = 0; i < dim; i++) {
    result[i] = 0.0;
    for (int j = 0; j < dim; j++) {
      result[i] += mat[i*dim + j]*vec[j];
    }
  }
  return;
}

double vector_dot_product(double* v1, double* v2, int dim) {
  double dot = 0.0;
  for (int i = 0; i < dim; i++) {
    dot += v1[i]*v2[i];
  }
  return dot;
}

void vector_difference(double *v1, double *v2, double *v3, int dim) {
  for (int i = 0; i < dim; i++) {
    v3[i] = v1[i] - v2[i];
  }
  return;
}

void vector_times_scalar(double *v1, double c1, int dim) {
  for (int i = 0; i < dim; i++) {
    v1[i] *= c1;
  }
  return;
}

void lanczos(double* h, double* alpha, double* beta, int dim, int n_iter) {
  double *vi = (double*) malloc(sizeof(double)*dim);
  double *h_tri = (double*) malloc(sizeof(double)*dim*n_iter);
  const gsl_rng_type *T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  for (int i = 0; i < dim; i++) {
    vi[i] = gsl_ran_gaussian(r, 1.0);
  }
  double norm = sqrt(vector_dot_product(vi, vi, dim));
  for (int i = 0; i < dim; i++) {
    vi[i] /= norm;
    h_tri[i] = vi[i];
  }
  double *wi = (double*) malloc(sizeof(double)*dim);
  square_matrix_times_vector(h, vi, wi, dim);
  alpha[0] = vector_dot_product(vi, wi, dim);
  vector_times_scalar(vi, alpha[0], dim);
  vector_difference(wi, vi, wi, dim);
  for (int i = 1; i < n_iter; i++) {
    beta[i-1] = sqrt(vector_dot_product(wi, wi, dim));
    if (beta[i-1] != 0.0) {
      vector_times_scalar(wi, 1.0/beta[i-1], dim);
      } else {
      printf("zero length\n");
      get_orthogonal_vector(h_tri, wi, dim, i);
      beta[i-1] = sqrt(vector_dot_product(wi, wi, dim));
      vector_times_scalar(wi, 1.0/beta[i-1], dim); 
    }
    for (int j = 0; j < dim; j++ ) {
      vi[j] = wi[j];
      h_tri[i*dim + j] = vi[j];
    }
    square_matrix_times_vector(h, vi, wi, dim);
    alpha[i] = vector_dot_product(vi, wi, dim);
    vector_times_scalar(vi, alpha[i], dim);
    vector_difference(wi, vi, wi, dim);
    for (int k = 0; k < i; k++) {
      double proj = 0.0;
      for (int j = 0; j < dim; j++) {
        proj += wi[j]*h_tri[k*dim + j];
      }
      for (int j = 0; j < dim; j++) {
        wi[j] -= proj*h_tri[k*dim + j];
      }
    }
    for (int k = 0; k < i; k++) {
      double proj = 0.0;
      for (int j = 0; j < dim; j++) {
        proj += wi[j]*h_tri[k*dim + j];
      }
      for (int j = 0; j < dim; j++) {
        wi[j] -= proj*h_tri[k*dim + j];
      }
    }
  }
  return;
}
