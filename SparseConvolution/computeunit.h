#ifndef COMPUTEUNIT_H
#define COMPUTEUNIT_H
#include "preproc.h"

class Kernel {
public:
  Kernel(const int _x_l, const int _y_l, const int _z_l) {
    int i, j, k, n = 0;
    x_l = _x_l;
    y_l = _y_l;
    z_l = _z_l;
    x_h = (x_l - 1) / 2;
    y_h = (y_l - 1) / 2;
    z_h = (z_l - 1) / 2;
    vol = x_l * y_l * z_l;
    // default kernel
    for (i = 0; i < KL_VOL; i++) {
      m[i] = 1;
    }

    for (i = -x_h; i <= x_h; i++)
      for (j = -y_h; j <= y_h; j++)
        for (k = -z_h; k <= z_h; k++) {
          Position pbuf(i, j, k);
          off[n] = pbuf;
          pbuf = pbuf + Position(1, 1, 1);
          // pbuf.print();
          // printf(": %d\n", n);
          n++;
        }
  }

  int x_l, y_l, z_l;       // xyz sizes of kernel
  int x_h, y_h, z_h;       // xyz half-lenghts
  int vol;                 // vol of kernel
  int x = 0, y = 0, z = 0; // xyz location of kernel
  float m[KL_VOL];         // kernel matrix
  Position off[KL_VOL];

private:
};

class ComputeUnit {
public:
  ComputeUnit(Position loc) { _loc = loc; }

  void addVox(Voxel vox) {
    _vox[_n] = vox;
    _n++;
  }
  void conv(Kernel kl) {
    // iterate over voxels
    for (int i = 0; i <= _n; i++) {
      // _vox[i] current voxel
      // the relative voxel position is
      Position rel = _vox[i].pos - _loc;
      // iterate over kernel offsets
      for (int j = 0; j < KL_VOL; j++) {
        auto wt = kl.m[j];
        // inverse for position
        Position pout = rel - kl.off[j];
        // check relative oob
        if (!isOOB(pout)) {
          // flatten
          int iout = flatten(pout);
          _out[iout] +=
              wt * (_vox[i].r_m + _vox[i].x_m + _vox[i].y_m + _vox[i].z_m);
        }
      }
    }
  }
  void write(const char *f) {
    std::ofstream of(f, std::ios_base::app);
    for (int i = 0; i < CT_VOL; i++) {
      if (_out[i] != 0) {
        Position pos = deflat(i);
        pos = pos + _loc;
        of << pos.x << ',' << pos.y << ',' << pos.z << ',' << _out[i] << '\n';
      }
    }
  }
  Position getCULoc() { return _loc; }
  unsigned short getVoxN() { return _n; }

private:
  int flatten(Position pos) {
    int ret;
    // x + (y*x_max) + (z*x_max*y_max)
    ret = pos.x + (pos.y * CT_SH_X) + (pos.z * CT_SH_X * CT_SH_Y);
    return ret;
  }
  Position deflat(int idx) {
    Position pos;
    pos.z = idx / (CT_SH_X * CT_SH_Y);
    idx -= (pos.z * CT_SH_X * CT_SH_Y);
    pos.y = idx / CT_SH_X;
    pos.x = idx % CT_SH_X;
    return pos;
  }
  bool isOOB(Position out) {
    bool isOOB = 1;
    if (out.x >= 0 && out.x < CT_SH_X)
      if (out.y >= 0 && out.y < CT_SH_Y)
        if (out.z >= 0 && out.z < CT_SH_Z)
          isOOB = 0;
    return isOOB;
  }
  Position _loc;          // absolute location of the unit
  unsigned short _n = 0;  // counter-1 for voxels
  Voxel _vox[MAX_VOX_CT]; // voxels inside shape and at the edge
  float _out[CT_VOL] = {0};
};

#endif // COMPUTEUNIT_H
