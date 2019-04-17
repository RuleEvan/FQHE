#include "lebedev.h"

gsl_complex y1sph(double theta, double phi) {

  gsl_complex i_phi = gsl_complex_rect(0, phi);
  gsl_complex y1 = gsl_complex_mul_real(gsl_complex_exp(i_phi), sin(theta));
  i_phi = gsl_complex_rect(0, -phi);
  y1 = gsl_complex_mul(y1, gsl_complex_exp(i_phi));
  y1 = gsl_complex_mul_real(y1, sin(theta));


  return y1;

}

double double_y1sph(double theta_1, double phi_1, double theta_2, double phi_2) {

  double y1 = pow(cos(theta_1), 2.0)*pow(cos(theta_2), 2.0);

  return y1;

}

gsl_complex wigner_d_max(int l, int m, double theta, double phi) {
 
  gsl_complex im_phi = gsl_complex_rect(0, m*phi);
  gsl_complex d = gsl_complex_rect(sqrt(gsl_sf_gamma(2*l + 1)/(gsl_sf_gamma(l + m + 1)*gsl_sf_gamma(l - m + 1))), 0.0);

  d = gsl_complex_mul_real(d, pow(cos(theta/2.0), l + m)*pow(sin(theta/2.0), l - m));
  d = gsl_complex_mul(d, gsl_complex_exp(im_phi));

  return d;
}

gsl_complex wigner_ortho(int l1, int l2, int m1, int m2, double theta, double phi) {

  gsl_complex d = gsl_complex_mul(wigner_d_max(l1, m1, theta, phi), gsl_complex_conjugate(wigner_d_max(l2, m2, theta, phi))); 

  return d;
} 

gsl_complex wigner_coulomb(int l1, int l2, int l3, int l4, int m1, int m2, int m3, int m4, double theta_1, double phi_1, double theta_2, double phi_2) {
  double softening = 0.00;
  gsl_complex d = gsl_complex_mul(gsl_complex_conjugate(wigner_d_max(l1, m1, theta_1, phi_1)), gsl_complex_conjugate(wigner_d_max(l2, m2, theta_2, phi_2)));
  d = gsl_complex_mul(d, wigner_d_max(l3, m3, theta_1, phi_1));
  d = gsl_complex_mul(d, wigner_d_max(l4, m4, theta_2, phi_2));
  d = gsl_complex_mul_real(d, 1.0/(sqrt(2)*sqrt(1.0 - sin(theta_1)*sin(theta_2)*cos(phi_1 - phi_2) - cos(theta_1)*cos(theta_2)) + softening));  
  return d;
}

gsl_complex lebedev( gsl_complex (*f) (int, int, int, int, double, double) , int l1, int l2, int m1, int m2) {

  gsl_complex result = gsl_complex_rect(0.0, 0.0);
  FILE *in_file;
  double theta, phi, coeff;
  in_file = fopen(LEB_FILE, "r");
  while (fscanf(in_file, "%lf %lf %lf\n", &phi, &theta, &coeff) == 3) {
    theta = theta/180.0*M_PI;
    phi = phi/180.0*M_PI;
    result = gsl_complex_add(result, gsl_complex_mul_real(f(l1, l2,m1, m2, theta, phi), coeff));
  }
  result = gsl_complex_mul_real(result, 4.0*M_PI);

  return result;
}


gsl_complex double_lebedev(gsl_complex (*f) (int, int, int, int, int, int, int, int, double, double, double, double), int l1, int l2, int l3, int l4, int m1, int m2, int m3, int m4 ) {
  // Integrates f(theta_1, phi_1, theta_2, phi_2) over the unit sphere

  gsl_complex result = gsl_complex_rect(0.0, 0.0);
  FILE *in_file;
  double *theta, *phi, *coeff;
  in_file = fopen(LEB_FILE, "r");
  theta = (double*) malloc(sizeof(double)*LEBEDEV_ORDER);
  phi = (double*) malloc(sizeof(double)*LEBEDEV_ORDER);
  coeff = (double*) malloc(sizeof(double)*LEBEDEV_ORDER);

  for (int i = 0; i < LEBEDEV_ORDER; i++) {
    fscanf(in_file, "%lf %lf %lf\n", &phi[i], &theta[i], &coeff[i]);
    theta[i] = theta[i]/180.0*M_PI;
    phi[i] = phi[i]/180.0*M_PI;
    
  }
  for (int i = 0; i < LEBEDEV_ORDER; i++) {
    double phi_1 = phi[i];
    double theta_1 = theta[i];
    double coeff_1 = coeff[i];

    for (int j = 0; j < LEBEDEV_ORDER; j++) {
      if (i == j) {continue;}
      double phi_2 = phi[j];
      double theta_2 = theta[j];
      double coeff_2 = coeff[j];

      result = gsl_complex_add(result, gsl_complex_mul_real( f(l1, l2, l3, l4, m1, m2, m3, m4, theta_1, phi_1, theta_2, phi_2), coeff_1*coeff_2) );
//      printf("%g %g\n", GSL_REAL(result), GSL_IMAG(result));
    }
  }
  result = gsl_complex_mul_real(result, 4.0*M_PI*4.0*M_PI);
  return result;
}
