#include "conv3d.h"
#include <iostream>
using namespace std;

int main() {
  DataImporter di;
  SparseTensor st = di.getTensor();
  // genIdxMatrix(&st);
  Kernel *kl = new Kernel;
  initKernel(kl, 3, 3, 3);
  ComputeTensor *ct_l[CT_NUM] = {NULL};
  initCTList(ct_l, &st, kl);
  // float cm[SH_X][SH_Y][SH_Z] = {0};

  conv4(ct_l, kl, "../SparseConvolution/data/conv4.csv");
  //  RuleBook *rb = genRuleBook(&st, kl);
  //  conv3(&st, kl, rb, "../SparseConvolution/data/conv3.csv");
  //   float ***c_m = conv2(&st, kl);
  //    matToCsv(c_m, st.sh[0], st.sh[1], st.sh[2],
  //            "../SparseConvolution/data/conv.csv");
  //    delete kl;
  //   genIdxPairs1(&st, &kl);
  //    Voxel ***c_m = submconv1(&st, &kl);
  //    for testing
  //    matToCsv(c_m, st.sh[0], st.sh[1], st.sh[2],
  //            "../SparseConvolution/data/submconv.csv");
}
