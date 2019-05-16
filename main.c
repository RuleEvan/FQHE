#include "hilbert.h"

int main(int argc, char *argv[]) {
//  hierarchy(1.0, 0.5, 6);
  wfnData *wd = read_binary_wfn_data("fqhe_b6_n5_laugh.wfn", "fqhe_b6_n5.bas", 0);
  wfnData *wd_shift = shift_op(wd, 0.5);
  printf("state: %g\n", wd_shift->bc[0]);
  normalize_wfn(wd_shift);
  print_wfn(wd_shift);
  return 0;
}
