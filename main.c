#include "hilbert.h"

int main(int argc, char *argv[]) {
//  wfnData* wd_shift =   hierarchy(0.5, 1.0, 6);
//  wfnData *wd = read_binary_wfn_data("fqhe_b5.5_n6.wfn", "fqhe_b5.5_n6.bas", 0);
//  wfnData *wd_shift = ud_tensor_ud(wd, 3, 3, 3, 0);
//  wfnData *wd_shift2 = ud_tensor_ud(wd, 3, 3, 6, 0);
//  wfnData *wd_shift3 = ud_tensor_ud(wd, 3, 4, 6, 0);
  //wfnData *wd_shift = shift_op(wd, 0.5, 0);
//  normalize_wfn(wd_shift);
//  normalize_wfn(wd_shift2);
//  normalize_wfn(wd_shift3);
//  wfn_mult(wd_shift1, -3.78448);
//  wfn_sum(wd_shift1, wd_shift2, 4.15488);
//  wfn_mult(wd_shift1, 1.02183);
//  wfn_sum(wd_shift1, wd_shift3, -1.82618);

//  print_wfn(wd_shift);
//  printf("\n");

//  wfnData *wde = read_binary_wfn_data("fqhe_b5.5_n6.wfn", "fqhe_b5.5_n6.bas", 3);

//  wde = order_wfn(wde);
//  print_wfn(wde);
//  printf("dot: %g\n", wfn_dot(wd_shift, wde));
//  generate_multilevel_interaction_file(5.5);
//  printf("\n");
  generate_interaction_file(7, 1);

  return 0;
}
