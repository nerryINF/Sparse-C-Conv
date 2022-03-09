/*
 * preproc.h
 *
 *  Created on: 21 nov. 2021
 *      Author: narimane
 */

#ifndef PREPROC_H_
#define PREPROC_H_
#include <fstream>
#include <stdbool.h>

#define SH_X 41
#define SH_Y 800
#define SH_Z 1408

#define MAX_KL 3 // for allocation
#define KL_X 3
#define KL_Y 3
#define KL_Z 3
#define KL_VOL (KL_X * KL_Y * KL_Z)

#define CT_DIV_X 2
#define CT_DIV_Y 2
#define CT_DIV_Z 2
#define CT_NUM (CT_DIV_X * CT_DIV_Y * CT_DIV_Z)

// shape for all ct's
#define CT_SH_X (SH_X / CT_DIV_X)
#define CT_SH_Y (SH_Y / CT_DIV_Y)
#define CT_SH_Z (SH_Z / CT_DIV_Z)
#define CT_VOL (CT_SH_X * CT_SH_Y * CT_SH_Z)
#define MAX_VOX 100000
#define MAX_VOX_CT 15000               // arbitrary
#define MAX_R_CT (MAX_VOX_CT * KL_VOL) // max number of rules

#include <iostream>
class Position {
public:
  Position() {
    x = 0;
    y = 0;
    z = 0;
  }
  Position(short _x, short _y, short _z) {
    x = _x;
    y = _y;
    z = _z;
  }

  Position operator+(Position _pos) {
    Position ret(x + _pos.x, y + _pos.y, z + _pos.z);
    return ret;
  }
  Position operator-(Position _pos) {
    Position ret(x - _pos.x, y - _pos.y, z - _pos.z);
    return ret;
  }
  void print() { std::cout << x << "," << y << "," << z; }
  short x, y, z;
};

class Voxel {
public:
  Position pos = {0, 0, 0};
  float x_m = 0.0, y_m = 0.0, z_m = 0.0,
        r_m = 0.0; // xyz medians and reflectance means
  int n = 0;       // number of points
};

struct SparseTensor {
  int num_vox;          // num_voxels (data_shape)
  int sh[3];            // sparse_shape
  unsigned char b_sz;   // batch_size
  Voxel vox[MAX_VOX];   // features_numpoints (voxels)
  Position in[MAX_VOX]; // indices, dynamically allocated
};

class DataImporter {
public:
  DataImporter() {
    std::ifstream fin("../SparseConvolution/data/data_shape.csv");
    fin >> st.num_vox;
    fin.close();
    read_indices();
    // read_sparse_shape();
    read_batch_size();
    read_features();
  }
  SparseTensor getTensor() { return st; }

private:
  void read_indices() {
    std::ifstream fin("../SparseConvolution/data/indices.csv");
    char cbuf;
    // iterate lines
    for (int i = 0; i < st.num_vox; i++) {
      // read line
      fin >> st.in[i].x >> cbuf >> st.in[i].y >> cbuf >> st.in[i].z;
      st.vox[i].pos = st.in[i];
    }
    fin.close();
  }
  void read_sparse_shape() {
    std::ifstream fin("../SparseConvolution/data/sparse_shape.csv");
    char cbuf;
    fin >> st.sh[0] >> cbuf >> st.sh[1] >> cbuf >> st.sh[2];
    fin.close();
  }
  void read_batch_size() {
    std::ifstream fin("../SparseConvolution/data/batch_size.csv");
    fin >> st.b_sz;
    fin.close();
  }
  void read_features() {
    std::ifstream fin("../SparseConvolution/data/features_numpoints.csv");
    char cbuf;
    // iterate lines
    for (int i = 0; i < st.num_vox; i++) {
      // read line
      fin >> st.vox[i].x_m >> cbuf >> st.vox[i].y_m >> cbuf >> st.vox[i].z_m >>
          cbuf >> st.vox[i].r_m >> cbuf >> st.vox[i].n;
    }
    fin.close();
  }
  SparseTensor st;
};

#endif /* PREPROC_H_ */
