#include "conv3d.h"
#include <iostream>
using namespace std;

int main() {
  DataImporter di;
  SparseTensor st = di.getTensor();
  genIdxMatrix(&st);
  Kernel kl = {3, 3, 3, 0, 0, 0, 0};
  kl.m = allocate3dMatrix(kl.x_l, kl.y_l, kl.z_l, 1);
  Voxel ***c_m = submconv1(&st, &kl);
  // for testing
  matToCsv(c_m, st.sh[0], st.sh[1], st.sh[2],
           "../SparseConvolution/data/submconv.csv");
}
