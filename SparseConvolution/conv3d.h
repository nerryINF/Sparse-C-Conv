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
  int x_l, y_l, z_l;       // xyz sizes of kernel
  int x_h, y_h, z_h;       // xyz half-lenghts
  int vol;                 // vol of kernel
  int x = 0, y = 0, z = 0; // xyz location of kernel
  int ***m;                // kernel matrix
};

int ***allocate3dMatrix(int n, int m, int l, int v) {
  int ***p;
  int i, j, k;

  p = (int ***)malloc(n * sizeof(int **));

  for (i = 0; i < n; i++)
    p[i] = (int **)malloc(m * sizeof(int *));

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      p[i][j] = (int *)malloc(l * sizeof(int));

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      for (k = 0; k < l; k++)
        p[i][j][k] = v;

  return p;
}

float ***allocate3dMatrix(int n, int m, int l, float v) {
  float ***p;
  int i, j, k;

  p = (float ***)malloc(n * sizeof(float **));

  for (i = 0; i < n; i++)
    p[i] = (float **)malloc(m * sizeof(float *));

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      p[i][j] = (float *)malloc(l * sizeof(float));

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      for (k = 0; k < l; k++)
        p[i][j][k] = v;

  return p;
}

void initKernel(Kernel *kl, int x_l, int y_l, int z_l) {
  kl->x_l = x_l;
  kl->y_l = y_l;
  kl->z_l = z_l;
  kl->x_h = (x_l - 1) / 2;
  kl->y_h = (y_l - 1) / 2;
  kl->z_h = (z_l - 1) / 2;
  kl->vol = x_l * y_l * z_l;
  kl->m = allocate3dMatrix(x_l, y_l, z_l, 1);
}

void addRule(int x, int y, int z, int i, int j, int k, int voxIdx, Kernel *kl,
             RuleBook *rb) {
  int p;
  // check that rule doesn't exist (according to vout coords)
  int ruleIdx = -1;
  for (p = 0; p <= rb->n; p++)
    if (x == rb->r[p].x && y == rb->r[p].y && z == rb->r[p].z)
      ruleIdx = p;
  // if rule doesn't exist, reallocate and add new rule
  if (ruleIdx < 0) {
    rb->n++;
    // allocate space for new rule
    rb->r = (Rule *)realloc(rb->r, (rb->n + 1) * sizeof(Rule));
    // assign vout coords for new rule
    rb->r[rb->n].x = x;
    rb->r[rb->n].y = y;
    rb->r[rb->n].z = z;
    // allocate space for rule matrix
    rb->r[rb->n].m = allocate3dMatrix(kl->x_l, kl->y_l, kl->z_l, -1);
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
        if (x >= ct->loc[0] && x < ct->loc[0] + ct->sh[0])
          if (y >= ct->loc[1] && y < ct->loc[1] + ct->sh[1])
            if (z >= ct->loc[2] &&
                z < ct->loc[2] + ct->sh[2] + 1) // if in bounds add rule
              addRule(x, y, z, i + kl->x_h, j + kl->y_h, k + kl->z_h, voxIdx,
                      kl, rb);
      }
}

void genRuleBook(ComputeTensor *ct, Kernel *kl) {
  int l;
  RuleBook *rb = (RuleBook *)malloc(sizeof(RuleBook));
  // counter for Vouts (Vout key/index)
  rb->n = 0;
  // allocate
  rb->r = (Rule *)malloc((rb->n + 1) * sizeof(Rule));
  // allocate space for first rule matrix
  rb->r[rb->n].m = allocate3dMatrix(kl->x_l, kl->y_l, kl->z_l, -1);
  // init first rule
  rb->r[rb->n].x = ct->in[0].x;
  rb->r[rb->n].y = ct->in[0].y;
  rb->r[rb->n].z = ct->in[0].z;
  // fill rulebook
  for (l = 0; l < ct->n + 1; l++) {
    // center kernel at voxel
    kl->x = ct->in[l].x;
    kl->y = ct->in[l].y;
    kl->z = ct->in[l].z;
    voxToRuleBook(l, ct, kl, rb);
  }
  printf("rb->n : %d\n", rb->n);
  ct->rb = rb;
}

void conv4(ComputeTensor *ct, Kernel *kl, const char *f) {
  int i, j, k, p;
  float res;
  std::ofstream of(f, std::ios_base::app);

  for (p = 0; p <= ct->rb->n; p++) { // iterate rules
    res = 0;
    for (i = 0; i < kl->x_l; i++)
      for (j = 0; j < kl->y_l; j++)
        for (k = 0; k < kl->z_l; k++) {
          // check that idx isn't -1
          if (ct->rb->r[p].m[i][j][k] >= 0) {
            // calculate result
            res += ((ct->vox[ct->rb->r[p].m[i][j][k]].r_m +
                     ct->vox[ct->rb->r[p].m[i][j][k]].x_m +
                     ct->vox[ct->rb->r[p].m[i][j][k]].y_m +
                     ct->vox[ct->rb->r[p].m[i][j][k]].z_m) *
                    kl->m[i][j][k]);
          }
        }
    of << ct->rb->r[p].x << ',' << ct->rb->r[p].y << ',' << ct->rb->r[p].z
       << ',' << res << '\n';
  }
}

ComputeTensor **genCTList(SparseTensor *st, Kernel *kl, unsigned char div) {
  // shape for all ct's
  unsigned short ct_sh[3] = {(unsigned short)(st->sh[0] / div),
                             (unsigned short)(st->sh[1] / div),
                             (unsigned short)(st->sh[2] / div)};
  unsigned short i, j, k, m, p = 0;
  _Bool hasVoxel = 0;
  // allocate list to ct pointers (div^3 ct's)
  ComputeTensor **ct_l =
      (ComputeTensor **)calloc(div * div * div, sizeof(ComputeTensor *));

  for (i = 0; i < st->sh[0]; i += ct_sh[0])
    for (j = 0; j < st->sh[1]; j += ct_sh[1])
      for (k = 0; k < st->sh[2]; k += ct_sh[2]) {
        // init voxel count
        hasVoxel = 0;
        // iterate voxels
        for (m = 0; m < st->num_vox; m++) {
          if (st->in[m][0] >= i - kl->x_h &&
              st->in[m][0] <= i + ct_sh[0] + kl->x_h)
            if (st->in[m][1] >= j - kl->y_h &&
                st->in[m][1] <= j + ct_sh[1] + kl->y_h)
              if (st->in[m][2] >= k - kl->z_h &&
                  st->in[m][2] <= k + ct_sh[2] + kl->z_h) {
                // init tensor
                if (!hasVoxel) {
                  printf("%d\n", p);
                  ct_l[p] = (ComputeTensor *)malloc(sizeof(ComputeTensor));
                  hasVoxel = 1;
                  // set voxel index
                  ct_l[p]->n = 0;
                  // init shape
                  memcpy(ct_l[p]->sh, ct_sh, sizeof(ct_sh));
                  // set tensor location
                  ct_l[p]->loc[0] = i;
                  ct_l[p]->loc[1] = j;
                  ct_l[p]->loc[2] = k;
                  // init indice and feature lists
                  ct_l[p]->in = (Indice *)malloc(sizeof(Indice));
                  ct_l[p]->vox = (Voxel *)malloc(sizeof(Voxel));
                  //                  ct_l[p]->e = (_Bool
                  //                  *)malloc(sizeof(_Bool));
                } else { // tensor already init
                  // increment voxel index
                  ct_l[p]->n++;
                  // reallocate indice and feature lists
                  ct_l[p]->in = (Indice *)realloc(
                      ct_l[p]->in, (ct_l[p]->n + 1) * sizeof(Indice));
                  ct_l[p]->vox = (Voxel *)realloc(
                      ct_l[p]->vox, (ct_l[p]->n + 1) * sizeof(Voxel));
                  //                  ct_l[p]->e = (_Bool *)realloc(ct_l[p]->e,
                  //                  (ct_l[p]->n + 1) *
                  //                                                                sizeof(_Bool));
                }
                // set indices
                ct_l[p]->in[ct_l[p]->n].x = st->in[m][0];
                ct_l[p]->in[ct_l[p]->n].y = st->in[m][1];
                ct_l[p]->in[ct_l[p]->n].z = st->in[m][2];
                // set features
                ct_l[p]->vox[ct_l[p]->n].n = st->vox[m].n;
                ct_l[p]->vox[ct_l[p]->n].x_m = st->vox[m].x_m;
                ct_l[p]->vox[ct_l[p]->n].y_m = st->vox[m].y_m;
                ct_l[p]->vox[ct_l[p]->n].z_m = st->vox[m].z_m;
                ct_l[p]->vox[ct_l[p]->n].r_m = st->vox[m].r_m;
                // (TODO) FIX ELSEIF
                //                ct_l[p]->e[ct_l[p]->n] = 0;
                //                if (st->in[m][0] < i || st->in[m][0] > i +
                //                ct_sh[0] ||
                //                    st->in[m][1] < j || st->in[m][1] > j +
                //                    ct_sh[1] || st->in[m][2] < k ||
                //                    st->in[m][2] > k + ct_sh[2])
                //                  ct_l[p]->e[ct_l[p]->n] = 1; // EDGE
              }
        }

        /////////////////// CONVOLUTION //////////////////////
        if (hasVoxel) {
          genRuleBook(ct_l[p], kl);
          conv4(ct_l[p], kl, "../SparseConvolution/data/conv4.csv");
        }
        //////////////////////////////////////////////////////

        // increment tensor count
        p++;
      }

  return ct_l;
}

#endif /* CONV3D_H_ */
