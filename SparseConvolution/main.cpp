#include "conv3d.h"
#include <iostream>
#include <stdlib.h>
using namespace std;

int main() {
  DataImporter *di = new DataImporter;
  SparseTensor st = di->getTensor();
  // genIdxMatrix(&st);
  Kernel *kl = new Kernel(3, 3, 3);
  ComputeUnitList cu_l(kl);
  // iterate voxels
  int m;
  for (m = 0; m < st.num_vox; m++) {
    cu_l.addVox(&st.vox[m]);
  }
  cu_l.display();

  // conv
  system("rm ../SparseConvolution/data_test/conv8.csv");
  cu_l.conv("../SparseConvolution/data_test/conv8.csv");
  delete kl;
}
