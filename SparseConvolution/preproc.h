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
#include <string.h>

#define CEILING(x, y) (((x) + (y)-1) / (y))

#define SH_X 41
#define SH_Y 800
#define SH_Z 1408

int KL_VOL = 0;

#define CU_SH_X 10
#define CU_SH_Y 100
#define CU_SH_Z 100

#define CU_DIV_X CEILING(SH_X, CU_SH_X)
#define CU_DIV_Y CEILING(SH_Y, CU_SH_Y)
#define CU_DIV_Z CEILING(SH_Z, CU_SH_Z)
#define CU_NUM (CU_DIV_X * CU_DIV_Y * CU_DIV_Z)

#define CU_VOL (CU_SH_X * CU_SH_Y * CU_SH_Z)
#define MAX_VOX 100000
#define MAX_VOX_CU 1500 // arbitrary
// table containing CU shapes
int CU_SH[3] = {CU_SH_X, CU_SH_Y, CU_SH_Z};

#include <iostream>
/**
 * @brief The Position class
 */
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
  /**
   * @brief getVoxRelativePosition
   * @return Relative position of voxel % unit to which it belongs
   */
  Position getVoxRelativePosition() {
    int bufx, bufy, bufz;
    bufx = x % CU_SH_X;
    bufy = y % CU_SH_Y;
    bufz = z % CU_SH_Z;
    return Position(bufx, bufy, bufz);
  }
  /**
   * @brief getCURelativePosition
   * @param cu_pos position of unit
   * @return Relative position of voxel % unit passed in param
   */
  Position getCURelativePosition(Position cu_pos) {
    int bufx, bufy, bufz;
    bufx = (x - cu_pos.x);
    bufy = (y - cu_pos.y);
    bufz = (z - cu_pos.z);
    return Position(bufx, bufy, bufz);
  }
  /**
   * @brief toIdx : 3d to 1d index
   * @param sh_x : spatial shape in x
   * @param sh_y : spatial shape in y
   * @return  1d index (in a list)
   */
  int toIdx(int sh_x, int sh_y) { return x + (y * sh_x) + (z * sh_x * sh_y); }
  /**
   * @brief toCUIdx
   * @return Index of voxel in unit (uses toIx with sh_x=CU_SH_X and
   * sh_y=CU_SH_Y)
   */
  int toCUIdx() { return this->toIdx(CU_SH_X, CU_SH_Y); }
  /**
   * @brief toCUListIdx
   * @return index of unit in list of units
   */
  int toCUListIdx() { // assuming this is a CU position
    Position bufp(x / CU_SH_X, y / CU_SH_Y, z / CU_SH_Z);
    return bufp.toIdx(CU_DIV_X, CU_DIV_Y);
  }
  /**
   * @brief getCUPosition
   * @return _loc of unit to wich this(voxel) belongs
   */
  Position getCUPosition() { return *this - getVoxRelativePosition(); }
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
  /**
   * @brief init : initialization of voxel with setup of scen value
   * @param pos : position of voxel in the grid
   */
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
  float _x_m = 0.0, _y_m = 0.0, _z_m = 0.0,
        _r_m = 0.0; // xyz medians and reflectance means
  int _n = 0;       // number of points

private:
  /**
   * @brief _pos : Position of voxel
   */
  Position _pos;
  /**
   * @brief scen boundcheck scenario
   */
  int scen;
  /**
   * @brief boundCheck returns the position relative to its CU (value from 0 to
   * 26)
   */
  void boundCheck() {
    // scen = XYZ, 0 is in bounds, 1 is at edge below, 2 above
    short px = 1, py = 1, pz = 1;
    // if <=0 or >=shape-1  we are at an edge
    switch (_pos.getVoxRelativePosition().x) {
    case 0:
      px = 0;
      break;
    case CU_SH_X - 1:
      px = 2;
      break;
    }
    switch (_pos.getVoxRelativePosition().y) {
    case 0:
      py = 0;
      break;
    case CU_SH_Y - 1:
      py = 2;
      break;
    }
    switch (_pos.getVoxRelativePosition().z) {
    case 0:
      pz = 0;
      break;
    case CU_SH_Z - 1:
      pz = 2;
      break;
    }
    scen = px * 9 + py * 3 + pz;
  }
};

struct SparseTensor {
  int num_vox;        // num_voxels (data_shape)
  int sh[3];          // sparse_shape
  unsigned char b_sz; // batch_size
  Voxel *vox;         // features_numpoints (voxels)
  Position *in;       // indices, dynamically allocated
};

class DataImporter {
public:
  DataImporter(std::string folder) {
    _f = folder;
    std::ifstream fin(_f + "/data_shape.csv");
    fin >> st.num_vox;
    st.in = new Position[st.num_vox];
    st.vox = new Voxel[st.num_vox];
    //
    fin.close();
    read_indices();
    // read_sparse_shape();
    read_batch_size();
    read_features();
  }
  SparseTensor getTensor() { return st; }

private:
  void read_indices() {
    std::ifstream fin(_f + "/indices.csv");
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
    std::ifstream fin(_f + "/sparse_shape.csv");
    char cbuf;
    fin >> st.sh[0] >> cbuf >> st.sh[1] >> cbuf >> st.sh[2];
    fin.close();
  }
  void read_batch_size() {
    std::ifstream fin(_f + "/batch_size.csv");
    fin >> st.b_sz;
    fin.close();
  }
  void read_features() {
    std::ifstream fin(_f + "/features_numpoints.csv");
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

private:
  std::string _f;
};

#endif /* PREPROC_H_ */
