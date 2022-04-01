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

#define SH_X 8
#define SH_Y 8
#define SH_Z 8

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
#define MAX_VOX_CU 1500                // arbitrary
#define MAX_R_CT (MAX_VOX_CT * KL_VOL) // max number of rules
// table containing CU shapes
int CU_SH[3] = {CT_SH_X, CT_SH_Y, CT_SH_Z};

#include <iostream>
class Position {
public:
  Position() {
    x = 0;
    y = 0;
    z = 0;
  }
  Position(int _x, int _y, int _z) {
    x = _x;
    y = _y;
    z = _z;
  }
  Position toCURelativePosition() {
    int bufx, bufy, bufz;
    bufx = x % CT_SH_X;
    bufy = y % CT_SH_Y;
    bufz = z % CT_SH_Z;
    return Position(bufx, bufy, bufz);
  }
  Position toCURelativePosition(Position cu_pos) {
    int bufx, bufy, bufz;
    bufx = (x - cu_pos.x);
    bufy = (y - cu_pos.y);
    bufz = (z - cu_pos.z);
    return Position(bufx, bufy, bufz);
  }
  int toIdx(int sh_x, int sh_y, int sh_z) {
    return x + (y * sh_x) + (z * sh_x * sh_y);
  }
  int toCUIdx() { return this->toIdx(CT_SH_X, CT_SH_Y, CT_SH_Z); }
  int toCUListIdx() { // assuming this is a CU position
    Position bufp(x / CT_SH_X, y / CT_SH_Y, z / CT_SH_Z);
    return bufp.toIdx(CT_DIV_X, CT_DIV_Y, CT_DIV_Z);
  }
  Position toCUPosition() { return *this - toCURelativePosition(); }
  Position operator+(Position _pos) {
    Position ret(x + _pos.x, y + _pos.y, z + _pos.z);
    return ret;
  }
  Position operator-(Position _pos) {
    Position ret(x - _pos.x, y - _pos.y, z - _pos.z);
    return ret;
  }
  void print() { std::cout << x << "," << y << "," << z; }
  int x, y, z;
};

class Voxel {
public:
  void init(Position pos) {
    _pos = pos;
    boundCheck();
  }
  void setFeatures(float x_m, float y_m, float z_m, float r_m, int n) {
    _x_m = x_m;
    _y_m = y_m;
    _z_m = z_m;
    _r_m = r_m;
    _n = n;
  }
  Position getPos() { return _pos; }
  int getBoundCheck() { return scen; }
  // should be made private after cleanup
  float _x_m = 0.0, _y_m = 0.0, _z_m = 0.0,
        _r_m = 0.0; // xyz medians and reflectance means
  int _n = 0;       // number of points
  // int vout[KL_VOL] = {CT_VOL, CT_VOL, CT_VOL, CT_VOL, CT_VOL, CT_VOL, CT_VOL,
  //                     CT_VOL, CT_VOL, CT_VOL, CT_VOL, CT_VOL, CT_VOL, CT_VOL,
  //                     CT_VOL, CT_VOL, CT_VOL, CT_VOL, CT_VOL, CT_VOL, CT_VOL,
  //                     CT_VOL, CT_VOL, CT_VOL, CT_VOL, CT_VOL, CT_VOL};
  //********************************************

private:
  Position _pos;
  int scen; // boundcheck scenario
  // returns the position relative to its CU
  void boundCheck() {
    short px = 1, py = 1, pz = 1;
    // if <=0 or >=shape-1  we are at an edge
    switch (_pos.toCURelativePosition().x) {
    case 0:
      px = 0;
      break;
    case CT_SH_X - 1:
      px = 2;
      break;
    }
    switch (_pos.toCURelativePosition().y) {
    case 0:
      py = 0;
      break;
    case CT_SH_Y - 1:
      py = 2;
      break;
    }
    switch (_pos.toCURelativePosition().z) {
    case 0:
      pz = 0;
      break;
    case CT_SH_Z - 1:
      pz = 2;
      break;
    }
    scen = px * 9 + py * 3 + pz;
  }
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
    // TODO: wesh
    std::ifstream fin("../SparseConvolution/data_test/data_shape.csv");
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
    std::ifstream fin("../SparseConvolution/data_test/indices.csv");
    char cbuf;
    // iterate lines
    for (int i = 0; i < st.num_vox; i++) {
      // read line
      fin >> st.in[i].x >> cbuf >> st.in[i].y >> cbuf >> st.in[i].z;
      st.vox[i].init(st.in[i]);
    }
    fin.close();
  }
  void read_sparse_shape() {
    std::ifstream fin("../SparseConvolution/data_test/sparse_shape.csv");
    char cbuf;
    fin >> st.sh[0] >> cbuf >> st.sh[1] >> cbuf >> st.sh[2];
    fin.close();
  }
  void read_batch_size() {
    std::ifstream fin("../SparseConvolution/data_test/batch_size.csv");
    fin >> st.b_sz;
    fin.close();
  }
  void read_features() {
    std::ifstream fin("../SparseConvolution/data_test/features_numpoints.csv");
    char cbuf;
    float x_m, y_m, z_m, r_m;
    int n;
    // iterate lines
    for (int i = 0; i < st.num_vox; i++) {
      // read line
      fin >> x_m >> cbuf >> y_m >> cbuf >> z_m >> cbuf >> r_m >> cbuf >> n;
      st.vox[i].setFeatures(x_m, y_m, z_m, r_m, n);
    }
    fin.close();
  }
  SparseTensor st;
};

#endif /* PREPROC_H_ */
