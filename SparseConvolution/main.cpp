#include "conv3d.h"
#include <iostream>
#include <stdlib.h>
using namespace std;

int main() {
  DataImporter di;
  SparseTensor st = di.getTensor();
  // genIdxMatrix(&st);
  Kernel *kl = new Kernel(3, 3, 3);
  ComputeUnit *ct_l[CT_NUM] = {NULL};
  ComputeUnit *ct = 0;
  initCTList(ct_l, &st, kl);

  // conv
  system("rm ../SparseConvolution/data/conv7.csv");
  // iterate over ct's
  for (int m = 0; m < CT_NUM; m++) {
    if (ct_l[m]) {
      if (m == 122) {
        printf("\n");
      }
      ct = ct_l[m];
      ct->write("../SparseConvolution/data/conv7.csv");
    }
  }
}
