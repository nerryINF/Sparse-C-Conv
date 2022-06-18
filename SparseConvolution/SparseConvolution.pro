TEMPLATE = app
CONFIG += console c++17 optimize_full
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -fopenmp -funroll-loops
LIBS += -fopenmp

SOURCES += \
        main.cpp

HEADERS += \
    computeunit.h \
    conv3d.h \
    conv3d.h \
    data.h \
    preproc.h

DISTFILES += \
    ../../SparseConvolution/data/batch_size.csv \
    ../../SparseConvolution/data/data_shape.csv \
    ../../SparseConvolution/data/features_numpoints.csv \
    ../../SparseConvolution/data/indices.csv \
    ../../SparseConvolution/data/sparse_shape.csv
