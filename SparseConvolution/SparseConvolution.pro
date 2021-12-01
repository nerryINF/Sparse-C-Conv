TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp

HEADERS += \
    conv3d.h \
    conv3d.h \
    preproc.h

DISTFILES += \
    ../../SparseConvolution/data/batch_size.csv \
    ../../SparseConvolution/data/data_shape.csv \
    ../../SparseConvolution/data/features_numpoints.csv \
    ../../SparseConvolution/data/indices.csv \
    ../../SparseConvolution/data/sparse_shape.csv
