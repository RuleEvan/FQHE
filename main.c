#include "charge.h"

int main(int argc, char *argv[]) {
//  wfnData* wd_shift =   tensor_composite(1.5, 0.0, 4, 1.5, 4.5, 4, 0);
//  wfnData* wd_shift =   ph2_tensor_composite(1.5, 0.0, 4, 1.5, 1.5, 2.5, 2.5, 0, 4, 4, 0);

//  wfnData* wd_shift =   composite(2.5, 0.0, 6);
//    wfnData* wd = read_binary_wfn_data("fqhe_b4.5_n6.wfn", "fqhe_b4.5_n6.bas", 0);
//    wfnData *wd_shift = shift_op(wd, 0.0, 0);
//  wfnData *wd_shift = coherent_shift(wd, M_PI/2.0, 0.0, 1);
//  wd_shift = coherent_shift(wd_shift, 0, 0.0, 0);
//generate_lookup_table(7, 8.5);
  wfnData* wd_shift =  composite_fast(2, 0.0, 5);
  normalize_wfn(wd_shift);
//    wfnData* wd2 = read_binary_wfn_data("fqhe_b11.5_n9_m1.5.wfn", "fqhe_b11.5_n9_m1.5.bas", 0);
//    print_charge_density(wd, 100);
//    wd = order_wfn(wd);
//    print_wfn(wd_shift);
//   printf("dot: %g\n", pow(wfn_dot(wd, wd2), 2));

 // wfnData *wd; //  wfnData *wd_shift2 = u_tensor_d(wd, 4,0);
//printf("%g\n", composite_coulomb(wd, wd, 100));
//  print_coupled_charge_density("fqhe_b11.5_n9", 4.5, 4.5, 2, 2, 100, 10);
//  wfnData *wd_shift = ud_tensor_ud(wd, 3, 3, 6, 0);
//  wfnData *wd_shift3 = ud_tensor_ud(wd, 3, 4, 6, 0);
//  normalize_wfn(wd);
//  normalize_wfn(wd_shift);
//  normalize_wfn(wd_shift2);
//  normalize_wfn(wd_shift3);
//  wfn_mult(wd_shift1, -3.78448);
//  wfn_sum(wd_shift1, wd_shift2, 4.15488);
//  wfn_mult(wd_shift1, 1.02183);
//  wfn_sum(wd_shift1, wd_shift3, -1.82618);
//  wd = order_wfn(wd);
/*  int n_steps = 100;
  double total = 0.0;

  for (int i = 0; i < 100; i++) {
    double theta = i*M_PI/n_steps;
    double rho_couple = 2*M_PI*charge_density(wd, theta);
;
    total += rho_couple*sin(theta)*M_PI/n_steps;
    printf("%g, %g, %g, %g\n",  theta, rho_couple, sin(theta)*rho_couple, total);

  }
*/
/*
  double j_tot = 8;
  int i_pow = 2;
  for (int i = 0; i < 100; i++) {
    double theta = i*M_PI/n_steps;
    wd = read_binary_wfn_data("fqhe_b11.5_n9_m4.5.wfn", "fqhe_b11.5_n9_m4.5.bas", 0);
    double rho_couple = 0.0;

    double rho1 = charge_density(wd, theta);
    wd = read_binary_wfn_data("fqhe_b11.5_n9_m3.5.wfn", "fqhe_b11.5_n9_m3.5.bas", 0);
    double rho2 = charge_density(wd, theta);
    rho_couple += (rho1 + rho2)*pow(clebsch_gordan(4.5, 4.5, j_tot, 4.5, 3.5, 8), i_pow);

    wd = read_binary_wfn_data("fqhe_b12.5_n9_m3.5.wfn", "fqhe_b12.5_n9_m3.5.bas", 0);
    rho1 = charge_density(wd, theta);
    wd = read_binary_wfn_data("fqhe_b12.5_n9_m2.5.wfn", "fqhe_b12.5_n9_m2.5.bas", 0);
    rho2 = charge_density(wd, theta);
    rho_couple += (rho1 + rho2)*pow(clebsch_gordan(4.5, 4.5, j_tot, 3.5, 2.5, 6), i_pow);

    wd = read_binary_wfn_data("fqhe_b11.5_n9_m2.5.wfn", "fqhe_b11.5_n9_m2.5.bas", 0);
    rho1 = charge_density(wd, theta);
    wd = read_binary_wfn_data("fqhe_b11.5_n9_m1.5.wfn", "fqhe_b11.5_n9_m1.5.bas", 0);
    rho2 = charge_density(wd, theta);
    rho_couple += (rho1 + rho2)*pow(clebsch_gordan(4.5, 4.5, j_tot, 2.5, 1.5, 4), i_pow);

    wd = read_binary_wfn_data("fqhe_b12.5_n9_m1.5.wfn", "fqhe_b12.5_n9_m1.5.bas", 0);
    rho1 = charge_density(wd, theta);
    wd = read_binary_wfn_data("fqhe_b12.5_n9_mn1.5.wfn", "fqhe_b12.5_n9_mn1.5.bas", 0);
    rho2 = charge_density(wd, theta);
    rho_couple += (rho1 + rho2)*pow(clebsch_gordan(4.5, 4.5, j_tot, 1.5, -1.5, 0), i_pow);

    wd = read_binary_wfn_data("fqhe_b12.5_n9_m0.5.wfn", "fqhe_b12.5_n9_m0.5.bas", 0);
    rho1 = charge_density(wd, theta);
    wd = read_binary_wfn_data("fqhe_b12.5_n9_mn0.5.wfn", "fqhe_b12.5_n9_mn0.5.bas", 0);
    rho2 = charge_density(wd, theta);
    rho_couple += (rho1 + rho2)*pow(clebsch_gordan(4.5, 4.5, j_tot, 0.5, -0.5, 0), i_pow);


    rho_couple = 2*2*M_PI*10*(1.0/(4.0*M_PI) - rho_couple);

    total += rho_couple*sin(theta)*M_PI/n_steps;
    printf("%g, %g, %g, %g\n",  theta, rho_couple, sin(theta)*rho_couple, total);
  }
*/
//  print_wfn(wd);
//  printf("\n");
  wfnData *wde; 
  for (int i = 0; i < 10; i++) {
    wde = read_binary_wfn_data("fqhe_b6_n5.wfn", "fqhe_b6_n5.bas", i);
    double j = wde->j_nuc[i];
   wde = order_wfn(wde);
   printf("dot: %d %g %g\n",i, j, pow(wfn_dot(wd_shift, wde), 2));
  }
 //  printf("\n");
//  print_wfn(wde);
//  generate_multilevel_interaction_file(3.5, 2);
//  printf("\n");
//  generate_interaction_file(3.5, 1, 0);

  return 0;
}
