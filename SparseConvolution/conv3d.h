/*
 * conv3d.h
 *
 *  Created on: Nov 26, 2021
 *      Author: user
 */

#ifndef CONV3D_H_
#define CONV3D_H_
#include "preproc.h"
#include <omp.h>
#include <string.h>

struct Kernel {
  int x_l, y_l, z_l;             // xyz sizes of kernel
  int x_h, y_h, z_h;             // xyz half-lenghts
  int vol;                       // vol of kernel
  int x = 0, y = 0, z = 0;       // xyz location of kernel
  int m[MAX_KL][MAX_KL][MAX_KL]; // kernel matrix
  Position off[KL_VOL];
};

bool checkOOB(Position in) {
  bool ret = 1;
  if (in.x >= 0 && in.x < SH_X)
    if (in.y >= 0 && in.y < SH_Y)
      if (in.z >= 0 && in.z < SH_Z)
        ret = 0;
  return ret;
}

void initKernel(Kernel *kl, int x_l, int y_l, int z_l) {
  int i, j, k, n = 0;
  kl->x_l = x_l;
  kl->y_l = y_l;
  kl->z_l = z_l;
  kl->x_h = (x_l - 1) / 2;
  kl->y_h = (y_l - 1) / 2;
  kl->z_h = (z_l - 1) / 2;
  kl->vol = x_l * y_l * z_l;
  // default kernel
  for (i = 0; i < x_l; i++)
    for (j = 0; j < y_l; j++)
      for (k = 0; k < z_l; k++)
        kl->m[i][j][k] = 1;

  for (i = -kl->x_h; i <= kl->x_h; i++)
    for (j = -kl->y_h; j <= kl->y_h; j++)
      for (k = -kl->z_h; k <= kl->z_h; k++) {
        kl->off[n].x = i;
        kl->off[n].y = j;
        kl->off[n].z = k;
        n++;
      }
}

void addRule(int x, int y, int z, int i, int j, int k, int voxIdx, Kernel *kl,
             RuleBook *rb) {
  int p, l, m, n;
  // check that rule doesn't exist (according to vout coords)
  int ruleIdx = -1;
  for (p = 0; p <= rb->n; p++)
    if (x == rb->r[p].x && y == rb->r[p].y && z == rb->r[p].z)
      ruleIdx = p;
  // if rule doesn't exist, reallocate and add new rule
  if (ruleIdx < 0) {
    rb->n++;
    // assign vout coords for new rule
    rb->r[rb->n].x = x;
    rb->r[rb->n].y = y;
    rb->r[rb->n].z = z;
    // init rule matrix
    for (l = 0; l < MAX_KL; l++)
      for (m = 0; m < MAX_KL; m++)
        for (n = 0; n < MAX_KL; n++)
          rb->r[rb->n].m[l][m][n] = -1;
    // assign in rule matrix
    rb->r[rb->n].m[i][j][k] = voxIdx;
  } else {
    // assign in rule matrix
    rb->r[ruleIdx].m[i][j][k] = voxIdx;
  }
}

void voxToRuleBook(int voxIdx, ComputeTensor *ct, Kernel *kl, RuleBook *rb) {
  int x, y, z, i, j, k;
  for (i = -kl->x_h; i <= kl->x_h; i++)
    for (j = -kl->y_h; j <= kl->y_h; j++)
      for (k = -kl->z_h; k <= kl->z_h; k++) {
        // ALL ABSOLUTE COORDINATES COME FROM INVERTED KERNEL COORDINATES
        x = kl->x - i;
        y = kl->y - j;
        z = kl->z - k;
        // check out of bounds
        if (x >= ct->loc.x && x < ct->loc.x + CT_SH_X)
          if (y >= ct->loc.y && y < ct->loc.y + CT_SH_Y)
            if (z >= ct->loc.z &&
                z < ct->loc.z + CT_SH_Z) // if in bounds add rule
              addRule(x, y, z, i + kl->x_h, j + kl->y_h, k + kl->z_h, voxIdx,
                      kl, rb);
      }
}

void genRuleBook(ComputeTensor *ct, Kernel *kl) {
  int l;
  // counter for Vouts (Vout key/index)
  ct->rb.n = 0;
  // init first rule
  ct->rb.r[0].x = ct->in[0].x;
  ct->rb.r[0].y = ct->in[0].y;
  ct->rb.r[0].z = ct->in[0].z;
  // fill rulebook
  for (l = 0; l < ct->n + 1; l++) {
    // center kernel at voxel
    kl->x = ct->in[l].x;
    kl->y = ct->in[l].y;
    kl->z = ct->in[l].z;
    voxToRuleBook(l, ct, kl, &ct->rb);
  }
  printf("rb->n : %d\n", ct->rb.n);
}

void conv4(ComputeTensor *ct_l[CT_NUM], Kernel *kl, const char *f) {
  int i, j, k, p, m;
  float res;
  std::ofstream of(f, std::ios_base::app);
  ComputeTensor *ct = 0;
  // iterate over cts
  for (m = 0; m < CT_NUM; m++) {
    // check that ct is initialized
    if (ct_l[m]) {
      ct = ct_l[m];
      // printf("%d\n", ct->n);
      //  generate rulebook
      genRuleBook(ct, kl);
      // do convolution
      for (p = 0; p <= ct->rb.n; p++) { // iterate rules
        res = 0;
        for (i = 0; i < kl->x_l; i++)
          for (j = 0; j < kl->y_l; j++)
            for (k = 0; k < kl->z_l; k++) {
              // check that idx isn't -1
              if (ct->rb.r[p].m[i][j][k] >= 0) {
                // calculate result
                res += ((ct->vox[ct->rb.r[p].m[i][j][k]].r_m +
                         ct->vox[ct->rb.r[p].m[i][j][k]].x_m +
                         ct->vox[ct->rb.r[p].m[i][j][k]].y_m +
                         ct->vox[ct->rb.r[p].m[i][j][k]].z_m) *
                        kl->m[i][j][k]);
              }
            }
        of << ct->rb.r[p].x << ',' << ct->rb.r[p].y << ',' << ct->rb.r[p].z
           << ',' << res << '\n';
      }
    }
  }
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

void addInToCTList(Position ct_in, Position pos, Voxel *vox,
                   ComputeTensor *ct[CT_NUM], int ctl_i) {
  // ADD VOXEL TO CT
  // ct_n = ct_l[ct_l_i]->n;
  // check if ct is initialized
  if (!ct[ctl_i]) {
    // init ct
    ct[ctl_i] = (ComputeTensor *)malloc(sizeof(ComputeTensor));
    // set ct loc
    ct[ctl_i]->loc.x = ct_in.x * CT_SH_X;
    ct[ctl_i]->loc.y = ct_in.y * CT_SH_Y;
    ct[ctl_i]->loc.z = ct_in.z * CT_SH_Z;
  }
  ct[ctl_i]->vox[ct[ctl_i]->n] = *vox;
  ct[ctl_i]->in[ct[ctl_i]->n] = pos;
  ct[ctl_i]->n++;
}

void addVoxToCTList(Position pos, Voxel *vox, ComputeTensor *ct_l[CT_NUM],
                    Kernel *kl) {
  short ct_in[8] = {0}; // buffer for non-reapeating ct's for each voxel
  Position in_buf = {0, 0, 0};
  int ctl_i = 0;
  bool inList = 0;
  short i, m, n = 0;
  // iterate over kernel
  for (i = 0; i < KL_VOL; i++) {
    // in_buf contains the shifted indices
    in_buf.x = pos.x + kl->off[i].x;
    in_buf.y = pos.y + kl->off[i].y;
    in_buf.z = pos.z + kl->off[i].z;
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
        addInToCTList(voxToCTListIn(in_buf), pos, vox, ct_l, ctl_i);
        ct_in[n] = ctl_i;
        n++;
      }
    }
  }
}

void initCTList(ComputeTensor *ct_l[CT_NUM], SparseTensor *st, Kernel *kl) {
  unsigned short m, p;
  // iterate over list of voxels to place in ct
  for (m = 0; m < st->num_vox; m++) {
    addVoxToCTList(st->in[m], &st->vox[m], ct_l, kl);
  }
  for (p = 0; p < CT_NUM; p++) {
    if (ct_l[p]) {
      printf("%d,%d: %d,%d,%d\n", p, ct_l[p]->n, ct_l[p]->loc.x, ct_l[p]->loc.y,
             ct_l[p]->loc.z);
    }
  }
}

#endif /* CONV3D_H_ */
