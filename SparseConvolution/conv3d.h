/*
 * conv3d.h
 *
 *  Created on: Nov 26, 2021
 *      Author: user
 */

#ifndef CONV3D_H_
#define CONV3D_H_
#include "preproc.h"

struct Kernel {
  int x_l, y_l, z_l;       // xyz sizes of kernel
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

#endif /* CONV3D_H_ */
