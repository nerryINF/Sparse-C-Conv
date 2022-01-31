/*
 * conv3d.h
 *
 *  Created on: Nov 26, 2021
 *      Author: user
 */

#ifndef CONV3D_H_
#define CONV3D_H_
#include <omp.h>

#include "preproc.h"

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

Voxel ***allocate3dMatrix(int n, int m, int l, Voxel v) {
  Voxel ***p;
  int i, j, k;

  p = (Voxel ***)malloc(n * sizeof(Voxel **));

  for (i = 0; i < n; i++)
    p[i] = (Voxel **)malloc(m * sizeof(Voxel *));

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      p[i][j] = (Voxel *)malloc(l * sizeof(Voxel));

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      for (k = 0; k < l; k++)
        p[i][j][k] = v;

  return p;
}

void genIdxMatrix(SparseTensor *st) {
  int n, m, l, i, x, y, z;
  n = st->sh[0];
  m = st->sh[1];
  l = st->sh[2];

  // allocate and init matrix
  int ***IdxMatrix = allocate3dMatrix(n, m, l, -1);

  // fill matrix
  for (i = 0; i < st->num_vox; i++) {
    x = st->in[i][0];
    y = st->in[i][1];
    z = st->in[i][2];
    IdxMatrix[x][y][z] = i;
  }

  st->idxMat = IdxMatrix;
}

/*

Voxel convp(SparseTensor *st, Kernel *kl) { // single point convolution
  Voxel v = {0, 0, 0, 0, 0};
  int i, j, k, x_abs, y_abs, z_abs; // coords in kernel/absolute
  int idx = -1;
  int wt = 1;        // weight at point
  int hlx, hly, hlz; // kernel half lengths
  // int vol = (kl->x_l * kl->y_l * kl->z_l); // kernel total volume
  hlx = (kl->x_l - 1) / 2;
  hly = (kl->y_l - 1) / 2;
  hlz = (kl->z_l - 1) / 2;
  for (i = -hlx; i <= hlx; i++)
    for (j = -hly; j <= hly; j++)
      for (k = -hlz; k <= hlz; k++) {
        x_abs = kl->x + i;
        y_abs = kl->y + j;
        z_abs = kl->z + k;
        if (x_abs >= 0 && x_abs < st->sh[0])
          if (y_abs >= 0 && y_abs < st->sh[1])
            if (z_abs >= 0 && z_abs < st->sh[2]) {
              idx = st->idxMat[x_abs][y_abs][z_abs];
              if (idx >= 0) {
                wt = kl->m[i + hlx][j + hly][k + hlz];
                v.n += st->vox[idx].n * wt;
                v.r_m += st->vox[idx].r_m * wt;
                v.x_m += st->vox[idx].x_m * wt;
                v.y_m += st->vox[idx].y_m * wt;
                v.z_m += st->vox[idx].z_m * wt;
              }
            }
      }
  return v;
}

Voxel ***conv1(SparseTensor *st, Kernel *kl) { // convolution v1
  int n, m, l;
  n = st->sh[0];
  m = st->sh[1];
  l = st->sh[2];

  // allocate convolution matrix
  Voxel ***c_m = allocate3dMatrix(n, m, l, {0, 0, 0, 0, 0});

  // convolution
  for (kl->x = 0; kl->x < n; kl->x++)
    for (kl->y = 0; kl->y < m; kl->y++)
      for (kl->z = 0; kl->z < l; kl->z++)
        c_m[kl->x][kl->y][kl->z] = convp(st, kl);

  return c_m;
}

Voxel ***submconv1(SparseTensor *st, Kernel *kl) {
  int i, n, m, l;
  n = st->sh[0];
  m = st->sh[1];
  l = st->sh[2];

  // allocate convolution matrix
  Voxel ***c_m = allocate3dMatrix(n, m, l, {0, 0, 0, 0, 0});

  // subm convolution
  for (i = 0; i < st->num_vox; i++) {
    kl->x = st->in[i][0];
    kl->y = st->in[i][1];
    kl->z = st->in[i][2];
    c_m[kl->x][kl->y][kl->z] = convp(st, kl);
  }

  return c_m;
}
// first method for generating indice pairs, assumes that memory access is less
// costly, requires the sparse tensor index matrix to be generated
void genIdxPairs1(SparseTensor *st, Kernel *kl) {
  int idx, i, j, k, n;
  int x_abs, y_abs, z_abs;               // absolute position for kernel weights
  int vol = kl->x_l * kl->y_l * kl->z_l; // kernel volume
  int hvol = (vol - 1) / 2;              // half kernel volume
  int hlx, hly, hlz;                     // kernel half lengths xyz
  hlx = (kl->x_l - 1) / 2;
  hly = (kl->y_l - 1) / 2;
  hlz = (kl->z_l - 1) / 2;
  // allocate Indice Pair matrix
  int ***ip = allocate3dMatrix(st->num_vox, 2, vol, -1);
  for (idx = 0; idx < st->num_vox; idx++) {
    // index for advancing in third dimension
    n = 0;
    // set kernel position
    kl->x = st->in[idx][0];
    kl->y = st->in[idx][1];
    kl->z = st->in[idx][2];
    // middle pair, equal indices
    ip[idx][0][hvol] = idx;
    ip[idx][1][hvol] = idx;
    // first (and second) half of pairs
    for (i = -hlx; i <= hlx; i++)
      for (j = -hly; j <= hly; j++)
        for (k = -hlz; k <= hlz; k++) {
          if (i == 0 && j == 0 && k == 0)
            goto idxdone; // only calculate half the indice pairs
          x_abs = kl->x + i;
          y_abs = kl->y + j;
          z_abs = kl->z + k;
          // check out of bounds
          if (x_abs >= 0 && x_abs < st->sh[0])
            if (y_abs >= 0 && y_abs < st->sh[1])
              if (z_abs >= 0 && z_abs < st->sh[2]) {
                ip[idx][0][n] = idx;
                ip[idx][1][vol - n - 1] = idx;
                ip[idx][0][vol - n - 1] = st->idxMat[x_abs][y_abs][z_abs];
                ip[idx][1][n] = st->idxMat[x_abs][y_abs][z_abs];
              }
          n++;
        }
  idxdone:;
    ;
  }
  st->idxPairs = ip;
}

Voxel ***submconv2(SparseTensor *st, Kernel *kl) {
  int idx, n, m, l;
  int i, j, k, x_abs, y_abs, z_abs; // coords in kernel/absolute
  int hlx, hly, hlz;                // kernel half lengths xyz
  hlx = (kl->x_l - 1) / 2;
  hly = (kl->y_l - 1) / 2;
  hlz = (kl->z_l - 1) / 2;
  // main and secondary indices (in index pairs)
  int m_idx, s_idx;
  // kernel weight
  int m_wt, s_wt = 1; // main and secondary weights
  // set main weight (invariable)
  m_wt = kl->m[hlx][hly][hlz];
  Voxel m_v, s_v; // main and secondary voxels
  // convmat size
  n = st->sh[0];
  m = st->sh[1];
  l = st->sh[2];
  // allocate convolution matrix
  Voxel ***c_m = allocate3dMatrix(n, m, l, {0, 0, 0, 0, 0});
  for (idx = 0; idx < st->num_vox; idx++) {
    // index for advancing in third dimension
    n = 0;
    // iterate all kernel weights
    for (i = -hlx; i <= hlx; i++)
      for (j = -hly; j <= hly; j++)
        for (k = -hlz; k <= hlz; k++) {
          if (st->idxPairs[idx][0][n] >= 0 && st->idxPairs[idx][1][n] >= 0) {
            // CORRIGER VERIFICATION POUR NON SUBM

            // set secondary kernel weight
            s_wt = kl->m[i + hlx][j + hly][k + hlz];
            // set main and secondary indexes from pairs
            m_idx = st->idxPairs[idx][0][n];
            s_idx = st->idxPairs[idx][1][n];
            // set kernel position at main index
            kl->x = st->in[m_idx][0];
            kl->y = st->in[m_idx][1];
            kl->z = st->in[m_idx][2];
            // absolute location of current convolution operation
            x_abs = kl->x + i;
            y_abs = kl->y + j;
            z_abs = kl->z + k;
            // main voxel
            m_v = st->vox[m_idx];
            s_v = st->vox[s_idx];
            // calculate convolution
            c_m[x_abs][y_abs][z_abs].n += (m_v.n * m_wt + s_v.n * s_wt);
            c_m[x_abs][y_abs][z_abs].r_m += (m_v.r_m * m_wt + s_v.r_m * s_wt);
            c_m[x_abs][y_abs][z_abs].x_m += (m_v.x_m * m_wt + s_v.x_m * s_wt);
            c_m[x_abs][y_abs][z_abs].y_m += (m_v.y_m * m_wt + s_v.y_m * s_wt);
            c_m[x_abs][y_abs][z_abs].z_m += (m_v.z_m * m_wt + s_v.z_m * s_wt);
            // TODO: USE LABEL TO HALF MEMORY ACCESSES, INVERSE MAIN AND
            // SECONDARY WEIGHTS
          }
          n++;
        }
    // doneconvp:;;
  }
  return c_m;
}
*/

void matToCsv(Voxel ***v, int n, int m, int l, const char *f) {
  std::ofstream of(f);
  int i, j, k;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      for (k = 0; k < l; k++)
        of << i << ',' << j << ',' << k << ',' << v[i][j][k].x_m << ','
           << v[i][j][k].y_m << ',' << v[i][j][k].z_m << ',' << v[i][j][k].r_m
           << ',' << v[i][j][k].n << '\n';
}

void matToCsv(float ***c_m, int n, int m, int l, const char *f) {
  std::ofstream of(f);
  int i, j, k;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      for (k = 0; k < l; k++)
        of << i << ',' << j << ',' << k << ',' << c_m[i][j][k] << '\n';
}

/* conp2, final sparse non-subm sub-function
 * c_m  : result matrix
 * in   : index of current voxel in st->vox
 * st   : sparse tensor
 * kl   : convolution kernel
 */
void convp2(float ***c_m, int in, SparseTensor *st, Kernel *kl) {
  int i, j, k; // coords relative to kernel of the output point
  int x, y, z; // buffer for absolute coordinates
  int wt, res; // buffer for weight
  for (i = -kl->x_h; i <= kl->x_h; i++)
    for (j = -kl->y_h; j <= kl->y_h; j++)
      for (k = -kl->z_h; k <= kl->z_h; k++) {
        // ALL ABSOLUTE COORDINATES COME FROM INVERTED KERNEL COORDINATES
        x = kl->x - i;
        y = kl->y - j;
        z = kl->z - k;
        // get weight (SAME FOR ALL FEATURES FOR TESTING PURPOSES)
        wt = kl->m[i + kl->x_h][j + kl->y_h][k + kl->z_h];
        // check out of bounds
        if (x >= 0 && x < st->sh[0])
          if (y >= 0 && y < st->sh[1])
            if (z >= 0 && z < st->sh[2]) {
              res = wt * (st->vox[in].r_m + st->vox[in].x_m + st->vox[in].y_m +
                          st->vox[in].z_m);
              c_m[x][y][z] += res;
            }
      }
}

/*
void convp3(float *** c_m, int in, SparseTensor *st, Kernel *kl) { //
PARALLELIZATION int i, j, k; // coords relative to kernel of the output point
  int x, y, z; // buffer for absolute coordinates
  int wt,res;      // buffer for weight and result
  #pragma omp parallel for private(j,k)
  for (i = -kl->x_h; i <= kl->x_h; i++)
    for (j = -kl->y_h; j <= kl->y_h; j++)
      for (k = -kl->z_h; k <= kl->z_h; k++) {
        // ALL ABSOLUTE COORDINATES COME FROM INVERTED KERNEL COORDINATES
        x = kl->x - i;
        y = kl->y - j;
        z = kl->z - k;
        // get weight (SAME FOR ALL FEATURES FOR TESTING PURPOSES)
        wt = kl->m[i + kl->x_h][j + kl->y_h][k + kl->z_h];
        // check out of bounds
        if (x >= 0 && x < st->sh[0])
          if (y >= 0 && y < st->sh[1])
            if (z >= 0 && z < st->sh[2]) {
                res = wt*(st->vox[in].r_m
                          + st->vox[in].x_m
                          + st->vox[in].y_m
                          + st->vox[in].z_m);
                #pragma omp critical
                {
                c_m[x][y][z] += res;
                }
            }
      }
}
*/

float ***conv2(SparseTensor *st, Kernel *kl) { // final sparse convolution
  int i;

  // allocate result matrix
  float ***c_m = allocate3dMatrix(st->sh[0], st->sh[1], st->sh[2], (float)0);

  // convolution
  for (i = 0; i < st->num_vox; i++) {
    // center kernel at voxel
    kl->x = st->in[i][0];
    kl->y = st->in[i][1];
    kl->z = st->in[i][2];
    convp2(c_m, i, st, kl);
  }

  return c_m;
}

struct Rule {
  // the xyz Vout coordinates
  int x, y, z;
  // [3][3][3] index of Vin affected by each Kernel weight
  // index ties back to voxel list in sparsetensor
  int ***m;
};

struct RuleBook {
  // rules
  Rule *r;
  // number of rules
  int n;
};

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

void voxToRuleBook(int voxIdx, SparseTensor *st, Kernel *kl, RuleBook *rb) {
  int x, y, z, i, j, k;
  for (i = -kl->x_h; i <= kl->x_h; i++)
    for (j = -kl->y_h; j <= kl->y_h; j++)
      for (k = -kl->z_h; k <= kl->z_h; k++) {
        // ALL ABSOLUTE COORDINATES COME FROM INVERTED KERNEL COORDINATES
        x = kl->x - i;
        y = kl->y - j;
        z = kl->z - k;
        // check out of bounds
        if (x >= 0 && x < st->sh[0])
          if (y >= 0 && y < st->sh[1])
            if (z >= 0 && z < st->sh[2]) // if in bounds add rule
              addRule(x, y, z, i + kl->x_h, j + kl->y_h, k + kl->z_h, voxIdx,
                      kl, rb);
      }
}
RuleBook *genRuleBook(SparseTensor *st, Kernel *kl) {
  int l;
  RuleBook *rb = new RuleBook;
  // counter for Vouts (Vout key/index)
  rb->n = 0;
  // allocate
  rb->r = (Rule *)malloc((rb->n + 1) * sizeof(Rule));
  // allocate space for first rule matrix
  rb->r[rb->n].m = allocate3dMatrix(kl->x_l, kl->y_l, kl->z_l, -1);
  // init first rule
  rb->r[rb->n].x = st->in[0][0];
  rb->r[rb->n].y = st->in[0][1];
  rb->r[rb->n].z = st->in[0][2];
  // fill rulebook
  for (l = 0; l < st->num_vox; l++) {
    // center kernel at voxel
    kl->x = st->in[l][0];
    kl->y = st->in[l][1];
    kl->z = st->in[l][2];
    voxToRuleBook(l, st, kl, rb);
  }
  printf("%d\n", rb->n);
  return rb;
}

void conv3(SparseTensor *st, Kernel *kl, RuleBook *rb, const char *f) {
  int i, j, k, p;
  float res;
  std::ofstream of(f);
  for (p = 0; p <= rb->n; p++) { // iterate rules
    res = 0;
    for (i = 0; i < kl->x_l; i++)
      for (j = 0; j < kl->y_l; j++)
        for (k = 0; k < kl->z_l; k++) {
          // check that idx isn't -1
          if (rb->r[p].m[i][j][k] >= 0) {
            // calculate result
            res += ((st->vox[rb->r[p].m[i][j][k]].r_m +
                     st->vox[rb->r[p].m[i][j][k]].x_m +
                     st->vox[rb->r[p].m[i][j][k]].y_m +
                     st->vox[rb->r[p].m[i][j][k]].z_m) *
                    kl->m[i][j][k]);
          }
        }
    of << rb->r[p].x << ',' << rb->r[p].y << ',' << rb->r[p].z << ',' << res
       << '\n';
  }
}

#endif /* CONV3D_H_ */
