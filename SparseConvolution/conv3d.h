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

Position voxToCTListPos(Position pos) {
  Position ret;
  ret.x = pos.x / CT_SH_X;
  ret.y = pos.y / CT_SH_Y;
  ret.z = pos.z / CT_SH_Z;
  return ret;
}

void addInToCTList(Position pbuf, Voxel *vox, ComputeUnit *ct[CT_NUM]) {
  // ADD VOXEL TO CT
  // ct_n = ct_l[ct_l_i]->n;
  // check if ct is initialized
  Position ct_in = voxToCTListPos(pbuf);
  int ctl_i = voxToCTListIdx(pbuf);
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
        addInToCTList(in_buf, vox, ct_l);
        ct_in[n] = ctl_i;
        n++;
      }
    }
  }
}

int boundCheck(Voxel *vox) {
  int scen = 0;
  int px = 1, py = 1, pz = 1;
  // find cu coords
  Position cuPos = voxToCTListPos(vox->pos);
  // find relative vox coords
  Position relPos = vox->pos - cuPos;
  // if <=0 or >=shape-1  we are at an edge
  switch (relPos.x) {
  case 0:
    px = 0;
    break;
  case CT_SH_X - 1:
    px = 2;
    break;
  }
  switch (relPos.y) {
  case 0:
    py = 0;
    break;
  case CT_SH_Y - 1:
    py = 2;
    break;
  }
  switch (relPos.z) {
  case 0:
    pz = 0;
    break;
  case CT_SH_Z - 1:
    pz = 2;
    break;
  }
  scen = px * 9 + py * 3 + pz;
  return scen;
}

Position comb(Position pos) {
  // if (pos.x != 1 && pos.y != 1 && pos.z != 1)
  std::cout << pos.x << "," << pos.y << "," << pos.z << std::endl;
  Position pos1 = pos;
  Position pos2 = pos;
  if (pos.x == 0 || pos.x == 2) {
    pos.x = 1;
    comb(pos);
  }
  if (pos.y == 0 || pos.y == 2) {
    pos1.y = 1;
    comb(pos1);
  }
  if (pos.z == 0 || pos.z == 2) {
    pos2.z = 1;
    comb(pos2);
  }
  return pos;
}

// 6+27 checks
void addVoxToCTList2(Voxel *vox, ComputeUnit *ct_l[CT_NUM], Kernel *kl) {
  // declare pos buffer
  Position pbase = vox->pos;
  Position pbuf = pbase;
  // add to main CU
  addInToCTList(pbuf, vox, ct_l);
  // XYZ, 0 is in bounds, 1 is at edge below, 2 above
  int scen = boundCheck(vox);
  switch (scen) {
  case 0:
    // 0 1 3 4 9 10 12
    pbuf = pbase + kl->off[0];
    addInToCTList(pbuf, vox, ct_l);
    pbuf = pbase + kl->off[1];
    addInToCTList(pbuf, vox, ct_l);
    pbuf = pbase + kl->off[0];
    addInToCTList(pbuf, vox, ct_l);
    pbuf = pbase + kl->off[1];
    addInToCTList(pbuf, vox, ct_l);
    pbuf = pbase + kl->off[0];
    addInToCTList(pbuf, vox, ct_l);
    pbuf = pbase + kl->off[1];
    addInToCTList(pbuf, vox, ct_l);

    break;
  }
}

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
