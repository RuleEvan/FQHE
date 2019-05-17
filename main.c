#include "hilbert.h"

int main(int argc, char *argv[]) {
//  hierarchy(2.0, 0.0, 5);
  wfnData *wd = read_binary_wfn_data("fqhe_b3_n4_laugh.wfn", "fqhe_b3_n4.bas", 0);
  wfnData *wd_shift = ud_tensor_ud(wd, 0, 4, 4, 0);
  //wfnData *wd_shift = shift_op(wd, 0.5, 0);
  normalize_wfn(wd_shift);
  print_wfn(wd_shift);

  return 0;
}
