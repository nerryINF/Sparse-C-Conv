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

#define CT_DIV_X 8
#define CT_DIV_Y 8
#define CT_DIV_Z 8
#define CT_NUM (CT_DIV_X * CT_DIV_Y * CT_DIV_Z)

// shape for all ct's
#define CT_SH_X (SH_X / CT_DIV_X)
#define CT_SH_Y (SH_Y / CT_DIV_Y)
#define CT_SH_Z (SH_Z / CT_DIV_Z)

#define MAX_VOX 100000
#define MAX_VOX_CT 1500                // arbitrary
#define MAX_R_CT (MAX_VOX_CT * KL_VOL) // max number of rules
struct Voxel {
  float x_m, y_m, z_m, r_m; // xyz medians and reflectance means
  int n;                    // number of points
};

struct Position {
  short x, y, z;
};

struct Rule {
  // the xyz Vout coordinates
  int x, y, z;
  // [3][3][3] index of Vin affected by each Kernel weight
  // index ties back to voxel list in sparsetensor
  // int m[KL_VOL];
  int m[3][3][3];
};

struct RuleBook {
  // rules
  Rule r[MAX_R_CT];
  // number of rules
  int n;
};

struct SparseTensor {
  int num_vox;        // num_voxels (data_shape)
  int sh[3];          // sparse_shape
  unsigned char b_sz; // batch_size
  Voxel vox[MAX_VOX]; // features_numpoints (voxels)
  Position in[MAX_VOX]; // indices, dynamically allocated
};

struct ComputeTensor {
  Position loc;            // absolute location of the tensor
  unsigned short n = 0;  // index for voxels
  Voxel vox[MAX_VOX_CT]; // voxels inside shape and at the edge (within kernel
                         // half-lenght)
  Position in[MAX_VOX_CT]; // indices for above voxels (NEW STRUCT)
                         //  _Bool *e;   // 1 if outside ct shape (edge), 0 if
                         //  inside
  RuleBook rb;           // rulebook
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
