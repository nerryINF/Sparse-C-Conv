#include "conv3d.h"
#include <iostream>
#include <stdlib.h>
using namespace std;

int main(int argc, char *argv[]) {
  // input args
  if (argc < 3) {
    cout << "not enough arguments" << endl;
    return 1;
  }
  string folder = argv[1];
  folder += "/";
  string out = argv[2];
  //
  DataImporter *di = new DataImporter(folder);
  SparseTensor st = di->getTensor();
  // genIdxMatrix(&st);
  Kernel *kl = new Kernel(folder);
  ComputeUnitList cu_l(kl);
  // iterate voxels
  int m;
  for (m = 0; m < st.num_vox; m++) {
    cu_l.addVox(&st.vox[m]);
  }
  cu_l.display();

  // conv
  system(("rm -f " + folder + out).c_str());
  cu_l.conv((folder + out).c_str());
  delete kl;
}
