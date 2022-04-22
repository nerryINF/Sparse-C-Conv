#ifndef COMPUTEUNIT_H
#define COMPUTEUNIT_H
#include "preproc.h"
#include <bits/stdc++.h>
/**
 * @brief The Kernel class
 */
class Kernel {
public:
  /**
   * @brief Kernel : loading kernel
   * @param _f : input file
   */
  Kernel(std::string _f) {
    std::ifstream fin(_f + "/kernel.csv");
    int i, j, k, n = 0;
    // read kernel vol in first line
    fin >> KL_VOL;
    // decalre
    m = new float[KL_VOL];
    off = new Position[KL_VOL];
    for (i = 0; i < KL_VOL; i++) {
      fin >> m[i];
    }

    x_l = std::cbrt(KL_VOL);
    y_l = x_l;
    z_l = x_l;
    x_h = (x_l - 1) / 2;
    y_h = (y_l - 1) / 2;
    z_h = (z_l - 1) / 2;
    vol = x_l * y_l * z_l;
    // gen offsets
    for (i = -x_h; i <= x_h; i++)
      for (j = -y_h; j <= y_h; j++)
        for (k = -z_h; k <= z_h; k++) {
          Position pbuf(i, j, k);
          off[n] = pbuf;
          n++;
        }
  }
  /**
   * @brief posToWt : from position to corresponding weight
   * @param pos : position from which to retrieve the weight
   * @return : weight at a position
   */
  float posToWt(Position pos) {
    pos = pos + Position(x_h, y_h, z_h);
    return m[pos.toIdx(x_l, y_l)];
  }
  int x_l, y_l, z_l;       // xyz sizes of kernel
  int x_h, y_h, z_h;       // xyz half-lenghts
  int vol;                 // vol of kernel
  int x = 0, y = 0, z = 0; // xyz location of kernel
  float *m;                // kernel matrix
  Position *off;
};

class ComputeUnit {
public:
  ComputeUnit(Kernel *kl, Position pos) {
    _kl = kl;
    _loc = pos;
  }
  /**
   * @brief procVoxel : process added voxel by calculation the convolution
   * @param vox : addede voxel to unit
   */
  void procVoxel(Voxel vox) {
    int i, j, k, outIdx;
    float wt;
    Position pbuf;
    // add to list
    vlist[_n] = vox;
    // initial checks
    int x_min_rel = vox.getPos().getCURelativePosition(_loc).x - _kl->x_h;
    int y_min_rel = vox.getPos().getCURelativePosition(_loc).y - _kl->y_h;
    int z_min_rel = vox.getPos().getCURelativePosition(_loc).z - _kl->z_h;
    int x_max_rel = vox.getPos().getCURelativePosition(_loc).x + _kl->x_h;
    int y_max_rel = vox.getPos().getCURelativePosition(_loc).y + _kl->y_h;
    int z_max_rel = vox.getPos().getCURelativePosition(_loc).z + _kl->z_h;
    // calculate for ranges
    forRange(x_min_rel, x_max_rel, 0);
    forRange(y_min_rel, y_max_rel, 1);
    forRange(z_min_rel, z_max_rel, 2);
    // for loops
    for (i = min[0]; i <= max[0]; i++)
      for (j = min[1]; j <= max[1]; j++)
        for (k = min[2]; k <= max[2]; k++) {
          //
          pbuf = Position(i, j, k) - vox.getPos().getCURelativePosition(_loc);
          wt = _kl->posToWt(pbuf); // corresponding weight
          pbuf = Position(i, j, k).getVoxRelativePosition();
          outIdx = pbuf.toCUIdx();
          _out[outIdx] += wt * (vox._r_m + vox._x_m + vox._y_m +
                                vox._z_m); // TO be changed for binary kernel
        }
    _n++;
  }
  /**
   * @brief write : write output file, result of conv
   * @param f : name of output file
   */
  void write(const char *f) {
    std::ofstream of(f, std::ios_base::app);
    for (int i = 0; i < CU_VOL; i++) {
      if (_out[i] != 0) {
        Position pos = toPos(i);
        pos = pos + _loc;
        of << pos.x << ',' << pos.y << ',' << pos.z << ',' << _out[i] << '\n';
      }
    }
  }
  /**
   * @brief reset : reset output to 0
   */
  void reset() {
    for (int i = 0; i < CU_VOL; i++) {
      _out[i] = {0};
    }
    _n = 0;
  }
  /**
   * @brief getCULoc
   * @return unit location
   */
  Position getCULoc() { return _loc; }
  /**
   * @brief getVoxN
   * @return number pf voxel in this unit
   */
  unsigned int getVoxN() { return _n; }

private:
  /**
   * @brief forRange : find range of for loop in procVox
   * @param min_rel : minimum relative position reached by the voxel
   * @param max_rel : maximum relative position reached by voxel
   * @param dim : 0 for x, 1 for y and 2 for z
   */
  void forRange(int min_rel, int max_rel, int dim) {
    min[dim] = min_rel;
    max[dim] = max_rel;
    if (min_rel < 0)
      min[dim] = 0;
    if (max_rel >= CU_SH[dim])
      max[dim] = CU_SH[dim] - 1;
  }
  /**
   * @brief toPos : from 1d list to 3d matrix
   * @param idx
   * @return a 3d position
   */
  Position toPos(int idx) {
    Position pos;
    pos.z = idx / (CU_SH_X * CU_SH_Y);
    idx -= (pos.z * CU_SH_X * CU_SH_Y);
    pos.y = idx / CU_SH_X;
    pos.x = idx % CU_SH_X;
    return pos;
  }
  Kernel *_kl;
  Position _loc;            // absolute location of the unit
  unsigned int _n = 0;      // counter-1 for voxels
  Voxel vlist[MAX_VOX_CU];  // list of voxels in unit
  float _out[CU_VOL] = {0}; // output of conv
  int min[3], max[3];       // ranges for loops
};

#endif // COMPUTEUNIT_H
