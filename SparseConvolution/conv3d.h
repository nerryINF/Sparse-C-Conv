/*
 * conv3d.h
 *
 *  Created on: Nov 26, 2021
 *      Author: user
 */

#ifndef CONV3D_H_
#define CONV3D_H_
#include "computeunit.h"
#include <string.h>

class ComputeUnitList {
public:
  ComputeUnitList(Kernel *kl) { _kl = kl; }
  /**
   * @brief addVox : add a voxel to one or many units depending on its impact
   * @param vox
   */
  void addVox(Voxel *vox) {
    _vox = vox;
    // add to main CU
    addToCuList(13);
    // XYZ, 0 is in bounds, 1 is at edge below, 2 above
    int scen = vox->getBoundCheck();
    switch (scen) {
    case 0: {
      multiAddVox(0, 1, 3, 4, 9, 10, 12);
      break;
    }
    case 1: {
      multiAddVox(1, 4, 10);
      break;
    }
    case 2: {
      multiAddVox(1, 2, 4, 5, 10, 11, 14);
      break;
    }
    case 3: {
      multiAddVox(3, 4, 12);
      break;
    }
    case 4: {
      multiAddVox(4);
      break;
    }
    case 5: {
      multiAddVox(4, 5, 14);
      break;
    }
    case 6: {
      multiAddVox(4, 3, 6, 7, 12, 15, 16);
      break;
    }
    case 7: {
      multiAddVox(4, 7, 16);
      break;
    }
    case 8: {
      multiAddVox(8, 17, 14, 16, 5, 4, 7);
      break;
    }
    case 9: {
      multiAddVox(9, 12, 10);
      break;
    }
    case 10: {
      multiAddVox(10);
      break;
    }
    case 11: {
      multiAddVox(11, 14, 10);
      break;
    }
    case 12: {
      multiAddVox(12);
      break;
    }
    case 14: {
      multiAddVox(14);
      break;
    }
    case 15: {
      multiAddVox(15, 12, 16);
      break;
    }
    case 16: {
      multiAddVox(16);
      break;
    }
    case 17: {
      multiAddVox(17, 14, 16);
      break;
    }
    case 18: {
      multiAddVox(18, 9, 12, 10, 21, 19, 22);
      break;
    }
    case 19: {
      multiAddVox(19, 10, 22);
      break;
    }
    case 20: {
      multiAddVox(20, 11, 14, 10, 23, 22, 19);
      break;
    }
    case 21: {
      multiAddVox(21, 12, 22);
      break;
    }
    case 22: {
      multiAddVox(22);
      break;
    }
    case 23: {
      multiAddVox(23, 14, 22);
      break;
    }
    case 24: {
      multiAddVox(24, 15, 12, 16, 21, 25, 22);
      break;
    }
    case 25: {
      multiAddVox(25, 16, 22);
      break;
    }
    case 26: {
      multiAddVox(26, 17, 14, 16, 23, 22, 25);
      break;
    }
    }
  }
  /**
   * @brief display : print CU location
   */
  void display() {
    int m;
    for (m = 0; m < CU_NUM; m++) {
      // check if ct is initialized
      if (CUs[m]) {
        CUs[m]->getCULoc().print();
        printf("\n%d\n", CUs[m]->getVoxN());
      }
    }
  }
  /**
   * @brief conv : convolution on non-empty units
   * @param f : name of output file
   */
  void conv(const char *f) {
    int m;
    for (m = 0; m < CU_NUM; m++) {
      // check if ct is initialized
      if (CUs[m]) {
        CUs[m]->write(f);
      }
    }
  }

private:
  /**
   * @brief addToCuList : add voxel to list of units
   * @param p : param of addVox (between 0 and 26) cf README
   */
  void addToCuList(short p) {
    // position of CU associated with offsetted position
    Position cu_pos = (_vox->getPos() + _kl->off[p]).getCUPosition();
    int ctl_i = cu_pos.toCUListIdx();
    // check if ct is initialized
    if (!CUs[ctl_i]) {
      // init ct
      CUs[ctl_i] = new ComputeUnit(_kl, cu_pos);
    }
    CUs[ctl_i]->procVoxel(*_vox);
  }
  /**
   * @brief multiAddVox : calls addToCuList with p0 (voxel is in a face)
   * @param p0
   */
  void multiAddVox(short p0) { addToCuList(p0); }
  /**
   * @brief multiAddVox calls addToCuList 3 times (voxel is in an endge)
   * @param p0
   * @param p1
   * @param p2
   */
  void multiAddVox(short p0, short p1, short p2) {
    int i;
    short posb[3] = {p0, p1, p2};
    for (i = 0; i < 3; i++)
      addToCuList(posb[i]);
  }
  /**
   * @brief multiAddVox : calls addToCuList 7 times (voxel is in a vertex)
   * @param p0
   * @param p1
   * @param p2
   * @param p3
   * @param p4
   * @param p5
   * @param p6
   */
  void multiAddVox(short p0, short p1, short p2, short p3, short p4, short p5,
                   short p6) {
    int i;
    short posb[7] = {p0, p1, p2, p3, p4, p5, p6};
    for (i = 0; i < 7; i++)
      addToCuList(posb[i]);
  }
  //**** variables
  ComputeUnit *CUs[CU_NUM] = {
      0};      // List of pointers to units, empty units points to 0
  Kernel *_kl; // Kernel
  Voxel *_vox; // voxel currently being computed
};

#endif /* CONV3D_H_ */
