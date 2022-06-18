#include "conv3d.h"
#include <iostream>
#include <stdlib.h>
using namespace std;

int main(int argc, char *argv[]) {
  // input args
  if (argc < 2) {
    cout << "not enough arguments" << endl;
    return 1;
  }
  string out = argv[1];
  // loading data
  DataImporter *di = new DataImporter();
  SparseTensor st = di->getTensor();
  Kernel *kl = new Kernel();
  ComputeUnitList cu_l(kl);
  // iterate voxels to dispatch in units
  int m;
  for (m = 0; m < st.num_vox; m++) {
    cu_l.addVox(&st.vox[m]);
  }
  cu_l.display();
  // conv
  system(("rm -f " + out).c_str());
  cu_l.conv((out).c_str());
  delete kl;
}
