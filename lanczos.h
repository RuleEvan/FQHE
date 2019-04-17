#ifndef LANCZOS_H
#define LANCZOS_H
#include "lebedev.h"
void square_matrix_times_vector(double *mat, double *vec, double *result, int dim);
void vector_times_scalar(double *v1, double c1, int dim);
void vector_difference(double *v1, double *v2, double *v3, int dim);
void lanczos(double *h, double *alpha, double *beta, int dim, int n_iter);
double vector_dot_product(double *v1, double *v2, int dim);
void matrix_elements(double *h, double *h_tri, int dim, int n_iter, double* alpha, double* beta);
void get_orthogonal_vector(double *h_tri, double *v, int dim, int n_done);
double characteristic_polynomial(double lambda, double *alpha, double *beta, int dim, int *num_sign);
void vector_max(double* v, double* max, int dim, int* n_max);
void vector_min(double* v, double* min, int dim, int* n_min);
void lanczos_eigenvalues(double*h, int dim, int n_iter);
#endif
