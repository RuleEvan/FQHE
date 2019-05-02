#ifndef LEBEDEV_H
#include "angular.h"

#define LEBEDEV_H

gsl_complex lebedev (gsl_complex (*f) (int, int, int, int, double, double), int l1, int l2, int m1, int m2);

gsl_complex y1sph(double theta, double phi);

gsl_complex wigner_d_max(int l, int m, double theta, double phi);

gsl_complex wigner_ortho(int l1, int l2, int m1, int m2, double theta, double phi);

gsl_complex wigner_coulomb(int l1, int l2, int l3, int l4, int m1, int m2, int m3, int m4, double theta_1, double phi_1, double theta_2, double phi_2);

gsl_complex double_lebedev (gsl_complex (*f) (int, int, int, int, int, int, int, int, double, double, double, double), int l1, int l2, int l3, int l4, int m1, int m2, int m3, int m4);

double double_y1sph(double theta_1, double phi_1, double theta_2, double phi_2);


#endif
