/*
 * conv3d.h
 *
 *  Created on: Nov 26, 2021
 *      Author: user
 */

#ifndef CONV3D_H_
#define CONV3D_H_
#include "computeunit.h"
#include <omp.h>
#include <string.h>

bool checkOOB(Position in) {
  bool ret = 1;
  if (in.x >= 0 && in.x < SH_X)
    if (in.y >= 0 && in.y < SH_Y)
      if (in.z >= 0 && in.z < SH_Z)
        ret = 0;
  return ret;
}

int voxToCTListIdx(Position pos) {
  Position ret;
  int ct_l_i;
  ret.x = pos.x / CT_SH_X;
  ret.y = pos.y / CT_SH_Y;
  ret.z = pos.z / CT_SH_Z;
  // x + (y*x_max) + (z*x_max*y_max)
  ct_l_i = ret.x + (ret.y * CT_DIV_X) + (ret.z * CT_DIV_X * CT_DIV_Y);
  return ct_l_i;
}

Position voxToCTListIn(Position pos) {
  Position ret;
  ret.x = pos.x / CT_SH_X;
  ret.y = pos.y / CT_SH_Y;
  ret.z = pos.z / CT_SH_Z;
  return ret;
}

void addInToCTList(Position ct_in, Voxel *vox, ComputeUnit *ct[CT_NUM],
                   int ctl_i) {
  // ADD VOXEL TO CT
  // ct_n = ct_l[ct_l_i]->n;
  // check if ct is initialized
  if (!ct[ctl_i]) {
    // init ct
    Position pos_ct(ct_in.x * CT_SH_X, ct_in.y * CT_SH_Y, ct_in.z * CT_SH_Z);
    ct[ctl_i] = new ComputeUnit(pos_ct);
  }
  ct[ctl_i]->addVox(*vox);
}

// 27* 8 checks
void addVoxToCTList(Voxel *vox, ComputeUnit *ct_l[CT_NUM], Kernel *kl) {
  short ct_in[8] = {0}; // buffer for non-reapeating ct's for each voxel
  Position in_buf = {0, 0, 0};
  int ctl_i = 0;
  bool inList = 0;
  short i, m, n = 0;
  // iterate over kernel
  for (i = 0; i < KL_VOL; i++) {
    // in_buf contains the shifted indices
    in_buf = vox->pos + kl->off[i];
    // ctl_i contains the CTList index for shifted indices
    ctl_i = voxToCTListIdx(in_buf);
    // checking if already in ct_in list
    inList = 0;
    for (m = 0; m < n; m++) {
      if (ct_in[m] == ctl_i) {
        inList = 1;
        break;
      }
    }
    // if not in list, check if out of bounds
    if (!inList) {
      if (!checkOOB(in_buf)) { // if not oob
        addInToCTList(voxToCTListIn(in_buf), vox, ct_l, ctl_i);
        ct_in[n] = ctl_i;
        n++;
      }
    }
  }
}

// 6+27 checks
void addVoxToCTList2(Voxel *vox, ComputeUnit *ct_l[CT_NUM], Kernel *kl) {}

void initCTList(ComputeUnit *ct_l[CT_NUM], SparseTensor *st, Kernel *kl) {
  unsigned short m, p;
  // iterate over list of voxels to place in ct
  for (m = 0; m < st->num_vox; m++) {
    addVoxToCTList(&st->vox[m], ct_l, kl);
  }
  for (p = 0; p < CT_NUM; p++) {
    if (ct_l[p]) {
      printf("%d,%d: %d,%d,%d\n", p, ct_l[p]->getVoxN(), ct_l[p]->getCULoc().x,
             ct_l[p]->getCULoc().y, ct_l[p]->getCULoc().z);
    }
  }
}

#endif /* CONV3D_H_ */
