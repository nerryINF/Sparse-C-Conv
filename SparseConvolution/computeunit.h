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
          //          pbuf = pbuf + Position(1, 1, 1);
          //          pbuf.print();
          //          printf("::::::::::::::::::::::::::: %d\n", n);
          //          comb(pbuf);
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
  ComputeUnit(Position loc, Kernel *kl) {
    _loc = loc;
    _kl = kl;
  }

  void procVox(Voxel vox) {
    // convolution
    // the relative voxel position is
    Position rel = vox.pos - _loc;
    // iterate over kernel offsets
    for (int j = 0; j < KL_VOL; j++) {
      auto wt = _kl->m[j];
      // inverse for position
      Position pout = rel - _kl->off[j];
      // check relative oob
      if (!isOOB(pout)) {
        // flatten
        int iout = flatten(pout);
        _out[iout] += wt * (vox.r_m + vox.x_m + vox.y_m + vox.z_m);
      }
    }
    _n++;
  }
  void conv() {}
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
  void reset() {
    for (int i = 0; i < CT_VOL; i++) {
      _out[i] = {0};
    }
    _n = 0;
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
    out = out + Position(1, 1, 1);
    int modx = (out.x % (CT_SH_X + 1));
    int mody = (out.y % (CT_SH_Y + 1));
    int modz = (out.z % (CT_SH_Z + 1));
    bool isIB = modx && mody && modz;
    return !isIB;
  }
  Kernel *_kl;
  Position _loc;         // absolute location of the unit
  unsigned short _n = 0; // counter-1 for voxels
  float _out[CT_VOL] = {0};
};

#endif // COMPUTEUNIT_H
