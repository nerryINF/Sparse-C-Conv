/*
 * preproc.h
 *
 *  Created on: 21 nov. 2021
 *      Author: narimane
 */

#ifndef PREPROC_H_
#define PREPROC_H_
#include <fstream>

struct Voxel {
  float x_m, y_m, z_m, r_m; // xyz medians and reflectance means
  int n;                    // number of points
};

struct SparseTensor {
  int num_vox;     // num_voxels (data_shape)
  int sh[3];       // sparse_shape
  int b_sz;        // batch_size
  Voxel *vox;      // features_numpoints (voxels), dynamically allocated
  int ***idxMat;   // dynamically allocated index matrix
  int ***idxPairs; // dynamically allocated index pair matrix
  int **in;        // indices, dynamically allocated
};

class DataImporter {
public:
  DataImporter() {
    std::ifstream fin("../SparseConvolution/data/data_shape.csv");
    fin >> st.num_vox;
    fin.close();
    read_indices();
    read_sparse_shape();
    read_batch_size();
    read_features();
  }
  SparseTensor getTensor() { return st; }

private:
  void read_indices() {
    std::ifstream fin("../SparseConvolution/data/indices.csv");
    char cbuf;
    int **in_buf = (int **)malloc(st.num_vox * sizeof(int *));
    // iterate lines
    for (int i = 0; i < st.num_vox; i++) {
      // allocate memory
      in_buf[i] = (int *)malloc(3 * sizeof(int));
      // read line
      fin >> in_buf[i][0] >> cbuf >> in_buf[i][1] >> cbuf >> in_buf[i][2];
    }
    fin.close();
    st.in = in_buf;
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
    Voxel *vbuf = (Voxel *)malloc(st.num_vox * sizeof(Voxel));
    // iterate lines
    for (int i = 0; i < st.num_vox; i++) {
      // read line
      fin >> vbuf[i].x_m >> cbuf >> vbuf[i].y_m >> cbuf >> vbuf[i].z_m >>
          cbuf >> vbuf[i].r_m >> cbuf >> vbuf[i].n;
    }
    fin.close();
    st.vox = vbuf;
  }
  SparseTensor st;
};

#endif /* PREPROC_H_ */
