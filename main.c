#include "hilbert.h"

int main(int argc, char *argv[]) {
//  hierarchy(1.0, 0.5, 6);
  wfnData *wd = read_binary_wfn_data("fqhe_b6_n5_laugh.wfn", "fqhe_b6_n5.bas", 0);
  shift_op(wd, 0.5);
  return 0;
}
