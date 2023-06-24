/*****************************************************************************
 * C++ functions to support MKDE computations in R via Rcpp
 * Author: Jeff A. Tracey, PhD; USGS-WERC San Diego Field Station
 *     jatracey@usgs.gov OR jatracey2005@gmail.com
 * Created: 23 August 2011
 * Revised: December 2013 - January 2014
*****************************************************************************/

// #define MATHLIB_STANDALONE // not needed with Rcpp


#ifndef _MKDE_H
#define _MKDE_H

#include <Rcpp.h>

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 1
    #define omp_get_thread_num() 0
#endif

#include <math.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>

// #include "Rmath.h" // not needed with Rcpp, use R namespace

const double MY_EPS = 0.00000001;
const double INT_EPS = 1.0e-6;
const double MY_PI = 3.141592653589793;

// for integrateNormal, phi
const double A1 = 0.31938153;
const double A2 = -0.356563782;
const double A3 = 1.781477937;
const double A4 = -1.821255978;
const double A5 = 1.330274429;
const double RSQRT2PI = 0.39894228040143267793994605993438; // 1/sqrt(2*pi)

const int JMAX = 20;

bool isMachineBigEndian(void);

// translation between coordinates and cell/voxel indexes
long getLinearIndex(long row, long col, long level, long nR, long nC);
int coordToIndex(double x, double minGridCoord, double cellSz);
double indexToCellCenterCoord(int i, double minGridCoord, double cellSz);
int getLowerCellIndex(double minZ, double minGridCoord, double cellSz);
int getUpperCellIndex(double maxZ, double minGridCoord, double cellSz);
double getLowerZ(double zGridMin, double zGridMax, double zMin, double cellSz);
double getUpperZ(double zGridMin, double zGridMax, double zMax, double cellSz);

// truncated normal
double integrateNormal(double x0, double x1, double mu, double sigma);
double dAtNormalDensityThreshold(double p, double sigma2);
double univariateNormalDensityThreshold(double p, double sigma2);
double reflectedNormalDensity(double x, double mu, double sigma2,
                              double xMin, double xMax, double thresh);

// for interaction model
double kernelBC(double x, double mu1, double sigma1sq, double mu2, double sigma2sq);
double trapzdKernelBC(double x0, double x1, double mu1, double sigma1sq, double mu2, double sigma2sq, int n);
double integrateKernelBC(double x0, double x1, double mu1, double sigma1, double mu2, double sigma2, double pThresh);


inline bool doubleEquals(double a, double b) {
  return fabs(a - b) < MY_EPS;
}

RcppExport SEXP mkde2dGrid02(SEXP obsT, SEXP obsX, SEXP obsY, SEXP useObs,
                             SEXP xgrid, SEXP ygrid, SEXP mvSig2xy,
                             SEXP obsSig2, SEXP tStep, SEXP pdfThresh);

RcppExport SEXP mkde2dGridv02interact(SEXP obsT, SEXP obsX0, SEXP obsY0,
                                      SEXP obsX1, SEXP obsY1, SEXP useObs,
                                      SEXP xgrid, SEXP ygrid,
                                      SEXP mvSig2xy0, SEXP mvSig2xy1,
                                      SEXP obsSig2xy0, SEXP obsSig2xy1,
                                      SEXP tStep, SEXP pdfThresh);

RcppExport SEXP mkde3dGridv02(SEXP obsT, SEXP obsX, SEXP obsY, SEXP obsZ, SEXP useObs,
                              SEXP xgrid, SEXP ygrid, SEXP zgrid,
                              SEXP zMinMatrix, SEXP zMaxMatrix,
                              SEXP mvSig2xy, SEXP mvSig2z, SEXP obsSig2, SEXP obsSig2z,
                              SEXP tStep, SEXP pdfThresh);

RcppExport SEXP mkde3dGridv02interact(SEXP obsT, SEXP obsX0, SEXP obsY0, SEXP obsZ0,
                                      SEXP obsX1, SEXP obsY1, SEXP obsZ1, SEXP useObs,
                                      SEXP xgrid, SEXP ygrid, SEXP zgrid,
                                      SEXP zMinMatrix, SEXP zMaxMatrix,
                                      SEXP mvSig2xy0, SEXP mvSig2xy1, SEXP mvSig2z0,
                                      SEXP mvSig2z1, SEXP obsSig2xy0, SEXP obsSig2xy1,
                                      SEXP obsSig2z0, SEXP obsSig2z1,
                                      SEXP tStep, SEXP pdfThresh);

RcppExport SEXP computeCellSurfaceArea(SEXP elevGrid, SEXP cellSize);

RcppExport SEXP writeMKDE3DtoVTK(SEXP xgrid, SEXP ygrid, SEXP zgrid, SEXP density,
                                 SEXP filename, SEXP description);

RcppExport SEXP writeMKDE3DtoGRASS(SEXP xgrid, SEXP ygrid, SEXP zgrid, SEXP density,
                                   SEXP filename, SEXP nodata);

RcppExport SEXP writeMKDE3DtoXDMF(SEXP xgrid, SEXP ygrid, SEXP zgrid, SEXP density,
                                  SEXP filenameXDMF, SEXP filenameDAT);

RcppExport SEXP writeRasterToXDMF02(SEXP xgrid, SEXP ygrid, SEXP rast,
                                  SEXP filenameXDMF, SEXP filenameDAT);

RcppExport SEXP writeRasterToVTK02(SEXP xgrid, SEXP ygrid, SEXP elev, SEXP rd,
                                 SEXP gr, SEXP bl, SEXP description, SEXP filenameVTK);

#endif
