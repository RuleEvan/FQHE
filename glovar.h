#include <stdio.h>
#include "math.h"
#include "stdlib.h"
#include "stdint.h"
#include "string.h"
#include "time.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>

// Isotope setup
#define M_S 6 // Number of magnetic flux quanta
#define N_PART 3 // Number of particles
#define N_LANCZOS 20

#define HASH_SIZE 1001467

#define LEBEDEV_ORDER 5810
#define LEB_FILE "lebedev_131.txt"

#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))



