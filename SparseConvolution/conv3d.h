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
  int wt = 1;                              // weight at point
  int hlx, hly, hlz;                       // kernel half lengths
  int vol = (kl->x_l * kl->y_l * kl->z_l); // kernel total volume
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
