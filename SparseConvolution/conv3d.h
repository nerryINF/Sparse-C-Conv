/*
 * conv3d.h
 *
 *  Created on: Nov 26, 2021
 *      Author: user
 */

#ifndef CONV3D_H_
#define CONV3D_H_
#include "computeunit.h"
#include <string.h>

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

void addInToCTList(Position pbuf, Voxel *vox, ComputeUnit *ct[CT_NUM],
                   Kernel *kl) {
  // ADD VOXEL TO CT
  // ct_n = ct_l[ct_l_i]->n;
  // check if ct is initialized
  Position ct_in = voxToCTListPos(pbuf);
  int ctl_i = voxToCTListIdx(pbuf);
  if (!ct[ctl_i]) {
    // init ct
    Position pos_ct(ct_in.x * CT_SH_X, ct_in.y * CT_SH_Y, ct_in.z * CT_SH_Z);
    ct[ctl_i] = new ComputeUnit(pos_ct, kl);
  }
  ct[ctl_i]->procVox(*vox);
}

void multiAddVox(Position pbase, Voxel *vox, ComputeUnit *ct_l[CT_NUM],
                 Kernel *kl, short p0, short p1, short p2, short p3, short p4,
                 short p5, short p6) {
  int i;
  short posb[7] = {p0, p1, p2, p3, p4, p5, p6};
  for (i = 0; i < 7; i++)
    addInToCTList(pbase + kl->off[posb[i]], vox, ct_l, kl);
}

void multiAddVox(Position pbase, Voxel *vox, ComputeUnit *ct_l[CT_NUM],
                 Kernel *kl, short p0, short p1, short p2) {
  int i;
  short posb[3] = {p0, p1, p2};
  for (i = 0; i < 3; i++)
    addInToCTList(pbase + kl->off[posb[i]], vox, ct_l, kl);
}

void multiAddVox(Position pbase, Voxel *vox, ComputeUnit *ct_l[CT_NUM],
                 Kernel *kl, short p0) {
  addInToCTList(pbase + kl->off[p0], vox, ct_l, kl);
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

// 6+27 checks
void addVoxToCTList(Voxel *vox, ComputeUnit *ct_l[CT_NUM], Kernel *kl) {
  // declare pos buffer
  Position pbase = vox->pos;
  Position pbuf = pbase;
  // add to main CU
  addInToCTList(pbuf, vox, ct_l, kl);
  // XYZ, 0 is in bounds, 1 is at edge below, 2 above
  int scen = boundCheck(vox);
  switch (scen) {
  case 0: {
    multiAddVox(pbase, vox, ct_l, kl, 0, 1, 3, 4, 9, 10, 12);
    break;
  }
  case 1: {
    multiAddVox(pbase, vox, ct_l, kl, 1, 4, 10);
    break;
  }
  case 2: {
    multiAddVox(pbase, vox, ct_l, kl, 1, 2, 4, 5, 10, 11, 14);
    break;
  }
  case 3: {
    multiAddVox(pbase, vox, ct_l, kl, 3, 4, 12);
    break;
  }
  case 4: {
    multiAddVox(pbase, vox, ct_l, kl, 4);
    break;
  }
  case 5: {
    multiAddVox(pbase, vox, ct_l, kl, 4, 5, 14);
    break;
  }
  case 6: {
    multiAddVox(pbase, vox, ct_l, kl, 4, 3, 6, 7, 12, 15, 16);
    break;
  }
  case 7: {
    multiAddVox(pbase, vox, ct_l, kl, 4, 7, 16);
    break;
  }
  case 8: {
    multiAddVox(pbase, vox, ct_l, kl, 8, 17, 14, 16, 5, 4, 7);
    break;
  }
  case 9: {
    multiAddVox(pbase, vox, ct_l, kl, 9, 12, 10);
    break;
  }
  case 10: {
    multiAddVox(pbase, vox, ct_l, kl, 10);
    break;
  }
  case 11: {
    multiAddVox(pbase, vox, ct_l, kl, 11, 14, 10);
    break;
  }
  case 12: {
    multiAddVox(pbase, vox, ct_l, kl, 12);
    break;
  }
  case 14: {
    multiAddVox(pbase, vox, ct_l, kl, 14);
    break;
  }
  case 15: {
    multiAddVox(pbase, vox, ct_l, kl, 15, 12, 16);
    break;
  }
  case 16: {
    multiAddVox(pbase, vox, ct_l, kl, 16);
    break;
  }
  case 17: {
    multiAddVox(pbase, vox, ct_l, kl, 17, 14, 16);
    break;
  }
  case 18: {
    multiAddVox(pbase, vox, ct_l, kl, 18, 9, 12, 10, 21, 19, 22);
    break;
  }
  case 19: {
    multiAddVox(pbase, vox, ct_l, kl, 19, 10, 22);
    break;
  }
  case 20: {
    multiAddVox(pbase, vox, ct_l, kl, 20, 11, 14, 10, 23, 22, 19);
    break;
  }
  case 21: {
    multiAddVox(pbase, vox, ct_l, kl, 21, 12, 22);
    break;
  }
  case 22: {
    multiAddVox(pbase, vox, ct_l, kl, 22);
    break;
  }
  case 23: {
    multiAddVox(pbase, vox, ct_l, kl, 23, 14, 22);
    break;
  }
  case 24: {
    multiAddVox(pbase, vox, ct_l, kl, 24, 15, 12, 16, 21, 25, 22);
    break;
  }
  case 25: {
    multiAddVox(pbase, vox, ct_l, kl, 25, 16, 22);
    break;
  }
  case 26: {
    multiAddVox(pbase, vox, ct_l, kl, 26, 17, 14, 16, 23, 22, 25);
    break;
  }
  }
}

void initCTList(ComputeUnit *ct_l[CT_NUM], SparseTensor *st, Kernel *kl) {
  unsigned short m, p;
  // iterate over list of voxels to place in ct
  for (m = 0; m < st->num_vox; m++) {
    if (m == 5765) {
      st->vox[m].pos.print();
      printf("\n");
    }
    addVoxToCTList(&st->vox[m], ct_l, kl);
  }
  //  for (p = 0; p < CT_NUM; p++) {
  //    if (ct_l[p]) {
  //      printf("%d,%d: %d,%d,%d\n", p, ct_l[p]->getVoxN(),
  //      ct_l[p]->getCULoc().x,
  //             ct_l[p]->getCULoc().y, ct_l[p]->getCULoc().z);
  //    }
  //  }
}

#endif /* CONV3D_H_ */
