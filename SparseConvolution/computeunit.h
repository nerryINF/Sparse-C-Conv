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
      m[i] = 1.0;
    }

    for (i = -x_h; i <= x_h; i++)
      for (j = -y_h; j <= y_h; j++)
        for (k = -z_h; k <= z_h; k++) {
          Position pbuf(i, j, k);
          off[n] = pbuf;
          n++;
        }
  }
  float posToWt(Position pos) {
    pos = pos + Position(1, 1, 1);
    return m[pos.toIdx(x_l, y_l, z_l)];
  }
  Position *getoff() { return off; }
  int x_l, y_l, z_l;       // xyz sizes of kernel
  int x_h, y_h, z_h;       // xyz half-lenghts
  int vol;                 // vol of kernel
  int x = 0, y = 0, z = 0; // xyz location of kernel
  float m[KL_VOL];         // kernel matrix
private:
  Position off[KL_VOL];
};

class ComputeUnit {
public:
  ComputeUnit(Kernel *kl, Position pos) {
    _kl = kl;
    _loc = pos;
  }

  void procVoxel(Voxel vox) {
    int i, j, k;
    float wt;
    Position pbuf;
    // add to list
    vlist[_n] = vox;
    // initial checks
    int x_min_rel = vox.getPos().toCURelativePosition(_loc).x - _kl->x_h;
    int y_min_rel = vox.getPos().toCURelativePosition(_loc).y - _kl->y_h;
    int z_min_rel = vox.getPos().toCURelativePosition(_loc).z - _kl->z_h;
    int x_max_rel = vox.getPos().toCURelativePosition(_loc).x + _kl->x_h;
    int y_max_rel = vox.getPos().toCURelativePosition(_loc).y + _kl->y_h;
    int z_max_rel = vox.getPos().toCURelativePosition(_loc).z + _kl->z_h;
    // calculate for ranges
    forRange(x_min_rel, x_max_rel, 0);
    forRange(y_min_rel, y_max_rel, 1);
    forRange(z_min_rel, z_max_rel, 2);
    // for loops
    for (i = min[0]; i <= max[0]; i++)
      for (j = min[1]; j <= max[1]; j++)
        for (k = min[2]; k <= max[2]; k++) {
          //
          pbuf = Position(i, j, k) - vox.getPos().toCURelativePosition(_loc);
          wt = _kl->posToWt(pbuf);
          pbuf = Position(i, j, k).toCURelativePosition();
          int outIdx = pbuf.toCUIdx();
          _out[outIdx] += wt * (vox._r_m + vox._x_m + vox._y_m + vox._z_m);
        }
    _n++;
  }
  void write(const char *f) {
    std::ofstream of(f, std::ios_base::app);
    for (int i = 0; i < CT_VOL; i++) {
      if (_out[i] != 0) {
        nonempty++;
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
  unsigned int getVoxN() { return _n; }

private:
  void forRange(int min_rel, int max_rel, int dim) {
    min[dim] = min_rel;
    max[dim] = max_rel;
    if (min_rel < 0)
      min[dim] = 0;
    if (max_rel >= CU_SH[dim])
      max[dim] = CU_SH[dim] - 1;
  }
  // TODO: move to Position class
  Position deflat(int idx) {
    Position pos;
    pos.z = idx / (CT_SH_X * CT_SH_Y);
    idx -= (pos.z * CT_SH_X * CT_SH_Y);
    pos.y = idx / CT_SH_X;
    pos.x = idx % CT_SH_X;
    return pos;
  }

  Kernel *_kl;
  Position _loc;       // absolute location of the unit
  unsigned int _n = 0; // counter-1 for voxels
  Voxel vlist[MAX_VOX_CU];
  float _out[CT_VOL + 1] = {0};
  int min[3], max[3];
  // debugging
  int nonempty = 0;
};

#endif // COMPUTEUNIT_H
