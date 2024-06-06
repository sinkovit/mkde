/*****************************************************************************
 * C++ functions to support MKDE computations in R via Rcpp
 * Author: Jeff A. Tracey, PhD; USGS-WERC San Diego Field Station
 *     jatracey@usgs.gov OR jatracey2005@gmail.com
 * Created: 23 August 2011
 * Revised: December 2013 - January 2014
 *****************************************************************************/

#include "mkdeFunctions.h"

using namespace std;

bool isMachineBigEndian(void) {
	union {
		long num;
		unsigned char uc[sizeof(long)];
	} u;
	u.num = 1;
	if (u.uc[0] == 1) {
		return false;
	} else {
		return true;
	}
}


/*****************************************************************************
 * Helper functions to translate between cell or voxel indexes and coordinates
 *****************************************************************************/

long getLinearIndex(long row, long col, long level, long nR, long nC) {
    return row + col*nR + level * nR * nC;
}

int coordToIndex(double x, double minGridCoord, double cellSz) {
    return (int)floor((x - minGridCoord - 0.5*cellSz)/cellSz);
}

double indexToCellCenterCoord(int i, double minGridCoord, double cellSz) {
    return minGridCoord + ((double)i)*cellSz;
}

// can switch floor and ceiling between the two functions to be more restrictive

int getLowerCellIndex(double minZ, double minGridCoord, double cellSz) {
    int i = coordToIndex(minZ, minGridCoord, cellSz);
    double cellZ = indexToCellCenterCoord(i, minGridCoord, cellSz);
    if (cellZ < minZ) {
        i++;
    }
    return i;
}

int getUpperCellIndex(double maxZ, double minGridCoord, double cellSz) {
    int i = coordToIndex(maxZ, minGridCoord, cellSz);
    double cellZ = indexToCellCenterCoord(i, minGridCoord, cellSz);
    if (cellZ > maxZ) {
        i--;
    }
    return i;
}

double getLowerZ(double zGridMin, double zGridMax, double zMin, double cellSz) {
    double zg0 = zGridMin - 0.5 * cellSz;
    if (zg0 >= zMin) {
        return zg0;
    } else {
        return zg0 + cellSz * floor((zMin - zg0) / cellSz);
    }
}

double getUpperZ(double zGridMin, double zGridMax, double zMax, double cellSz) {
    double zg0 = zGridMin - 0.5 * cellSz;
    double zg1 = zGridMax + 0.5 * cellSz;
    if (zg1 <= zMax) {
        return zg1;
    } else {
        return  zg0 + cellSz * (floor((zMax - zg0) / cellSz) + 1.0);
    }
}


/*****************************************************************************
 * Helper functions for truncated normal kernel
 *****************************************************************************/

// Uses R standalone math library
// get area under normal PDF from x0 to x1 with mean mu and standard deviation sigma
double integrateNormal(double x0, double x1, double mu, double sigma) {
    double p0 = R::pnorm(x0, mu, sigma, 1, 0);
    double p1 = R::pnorm(x1, mu, sigma, 1, 0);
    return p1 - p0;
}

/*****************************************************************************
 * Helper functions for spatial-temporal interaction
 *****************************************************************************/

/* Evaluate the square root of the product of the two kernels at x
 * We will very likely want to optimize the hell out of this function
 */
double kernelBC(double x, double mu1, double sigma1sq, double mu2, double sigma2sq) {
    double res = exp(-0.25*(x - mu1)*(x - mu1)/sigma1sq - 0.25*(x - mu2)*(x - mu2)/sigma2sq);
    return(res);
}

/* Based on trapzd() from Press et al.
 * Meant to be used only in a function like integrateKernelBC
 */
double trapzdKernelBC(double x0, double x1, double mu1, double sigma1sq, double mu2, double sigma2sq, int n, double old_s) {
    double x, tnm, sum, del;
    double s; // this was static in Press et al, but I made it so you have to pass the old value as an arg
    int it, j;
    if (n == 1) {
        s = 0.5*(x1 - x0)*(kernelBC(x1, mu1, sigma1sq, mu2, sigma2sq) + kernelBC(x0, mu1, sigma1sq, mu2, sigma2sq));
    } else {
        for (it = 1, j = 1; j < n - 1; j++) {
            it <<= 1;
        }
        tnm = it;
        del = (x1 - x0)/tnm;
        x = x0 + 0.5*del;
        for (sum = 0.0, j = 1; j <= it; j++, x += del) {
            sum += kernelBC(x, mu1, sigma1sq, mu2, sigma2sq);
        }
        s = 0.5*(old_s + (x1 - x0)*sum/tnm); // updates s if n > 1
    }
    return s;
}

/* Based on qsimp() from Press et al. */
double integrateKernelBC(double x0, double x1, double mu1, double sigma1, double mu2, double sigma2, double pThresh) {
    double sigma1sq = sigma1*sigma1, sigma2sq = sigma2*sigma2;
    double s, st, ost = 0.0, os = 0.0; // is this a local variable s?
    for (int j = 1; j <= JMAX; j++) {
        st = trapzdKernelBC(x0, x1, mu1, sigma1sq, mu2, sigma2sq, j, ost);
        s = (4.0*st - ost)/3.0;
        if (j > 5) {
            if (fabs(s - os) < INT_EPS*fabs(os) || (s == 0.0 && os == 0.0) || (st < pThresh)) {
                return s;
            } // also bail if st is less than minimum threshold so don't waste time on negligible associations
        }
        os = s;
        ost = st;
    }
    return 0.0;
}

/*****************************************************************************
 * Helper functions to determine range of cells or voxels over which to
 * update density.
 *****************************************************************************/

/*
* This function computes the absolute value of the distance from
* the mean at which the normal PDF achieves a specified density d
*/
double dAtNormalDensityThreshold(double p, double sigma2) {
    double y = -2.0*sigma2*log(p*2.0*sigma2*MY_PI);
    double res = NAN;
    if (y >= 0.0) {
        res = sqrt(y);
    }
    return res;
}

double univariateNormalDensityThreshold(double p, double sigma2) {
    double z = -2.0*sigma2*log(p*sqrt(2.0*sigma2*MY_PI));
    double res = NAN;
    if (z >= 0.0) {
        res = sqrt(z);
    }
    return res;
}


/*****************************************************************************
 * Helper functions for reflected normal kernel
 *****************************************************************************/

double reflectedNormalDensity(double x, double mu, double sigma,
                              double xMin, double xMax, double thresh) {
    double res = 0.0;
    int w0 = 0, w1 = 0, i = 1;
    double tmp = 2.0*thresh;
    if ((x >= xMin) && (x <= xMax) && (xMin < xMax)) {
        double d0 = x - xMin, d1 = xMax - x, x0 = 0.0, x1 = 0.0;
        res = R::dnorm(x, mu, sigma, 0);
        while (tmp > thresh) {
            i++; // start at 2
            w0 = i - i%2;
            w1 = (i - 1) - (i - 1)%2;
            x0 = x - ((double)w0)*d0 - ((double)w1)*d1;
            x1 = x + ((double)w0)*d1 + ((double)w1)*d0;
            tmp = R::dnorm(x0, mu, sigma, 0) + R::dnorm(x1, mu, sigma, 0);
            res += tmp;
        }
    } // else out of bounds and density is 0.0
    return res;
}

/*****************************************************************************
 * Functions that calculate MKDEs for a set of grid cells
 *****************************************************************************/

// 2D CASE
// obsT, obsX, obsY: observed data
// xyTable: a 2d array with pairs of (x,y) coordinates of cell centers
// tMax, tStep, mvSig2xy, obsSig2: parameters
RcppExport SEXP mkde2dGrid02(SEXP obsT, SEXP obsX, SEXP obsY, SEXP useObs,
                             SEXP xgrid, SEXP ygrid, SEXP mvSig2xy,
                             SEXP obsSig2, SEXP tStep, SEXP pdfThresh) {
    // observations
    Rcpp::NumericVector T(obsT); // observed times
    Rcpp::NumericVector X(obsX); // observed x-coordinates
    Rcpp::NumericVector Y(obsY); // observed y-coordinates
    Rcpp::IntegerVector isValid(useObs); // move step is flagged for use
    long nObs = T.length();
    int nObsStep;
    // grid speces
    Rcpp::NumericVector xGrid(xgrid); // cell center coordinates in the x-dimension
    Rcpp::NumericVector yGrid(ygrid); // cell center coordinates in the y-dimension
    int nX = xGrid.length();
    int nY = yGrid.length();
    double xSz = xGrid[1] - xGrid[0];
    double ySz = yGrid[1] - yGrid[0];
    // variance parameters (vectors)
    Rcpp::NumericVector msig2xy(mvSig2xy); // move variance in the(x,y)-dimensions
    Rcpp::NumericVector osig2xy(obsSig2); // observation error variance in the(x,y)-dimensions
    // control parameters
    Rcpp::NumericVector stepT(tStep); // integration time step
    Rcpp::NumericVector pdfMin(pdfThresh); // integration time step

    // arrays for MKDE computation
    double * ydens = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    std::vector<double> mkde(nX*nY);
    for (int i = 0; i < nX*nY; i++) {
        mkde[i] = 0.0;
    }

    // set up time variables
    double t0, t1, t, tOld, dt, alpha, totalT;
    // set up tmp variables
    double eX, eY;
    double sig2xy;
    double distMaxXY;
    int halo1, halo2, i1k, i2k;

    // start computing MKDE
    Rcpp::Rcout << "2D MKDE Computation: STARTING" << std::endl;
    nObsStep = nObs/50;
    if (nObsStep == 0) {
      nObsStep = 1;
    }
    totalT = 0.0;
    for (int j = 0; j < (nObs - 1); j++) {
      if (j % nObsStep == 0) {
	Rcpp::Rcout << "\tProcessing move step " << (j + 1) << " of " << (nObs - 1) << std::endl;
      }
        // report percent complete after each observation
        t0 = T[j];
        t1 = T[j + 1];
        dt = t1 - t0;
        t = t0;
        if (isValid[j] == 1) {
            bool exitLoop = false;
            bool finalLoop = false;
            totalT += dt; // add to total time iff move step used
            while (!exitLoop) { // iterate over integration time steps
                // Calculate fractional distance between t0 and current time
                alpha = (t - t0) / dt;
                // Calculate parameters for bilinear interpolation
                sig2xy = dt * alpha * (1.0 - alpha) * msig2xy[j] +
                    osig2xy[j] * (1.0 - alpha) * (1.0 - alpha) +
                    osig2xy[j+1] * alpha * alpha;
                // Get (x,y,z) coordinates of kernel origin using linear interpolation
                eX = X[j] + alpha * (X[j + 1] - X[j]);
                eY = Y[j] + alpha * (Y[j + 1] - Y[j]);
                // Convert to grid indices
                i1k = coordToIndex(eX, xGrid[0], xSz);
                i2k = coordToIndex(eY, yGrid[0], ySz);
                //
                distMaxXY = univariateNormalDensityThreshold(pdfMin[0], sig2xy);
                halo1 = (int)ceil(distMaxXY/xSz);
                halo2 = (int)ceil(distMaxXY/ySz);
                // Precompute voxel density in y dimension
                for(int i2 = 0; i2 < nY; i2++) {
                    ydens[i2] =  integrateNormal(yGrid[i2] - 0.5*ySz, yGrid[i2] + 0.5*ySz, eY, sqrt(sig2xy));
                }

                for (int i1 = std::max(0, i1k-halo1); i1 < std::min(nX, i1k+halo1); i1++) { // x-dimension
                    double voxx = xGrid[i1]; // voxel x
                    double xdens = integrateNormal(voxx - 0.5*xSz, voxx + 0.5*xSz, eX, sqrt(sig2xy));
                    for (int i2 = std::max(0, i2k-halo2); i2 < std::min(nY, i2k+halo2); i2++) { // y-dimension
                        double pXY = xdens * ydens[i2];
                        double tmpDens;
                        if (doubleEquals(t, t0)) { // first term
                            tmpDens = stepT[0] * pXY / 2.0;
                        } else if (doubleEquals(t, t1)) { // last term
                            tmpDens = (t - tOld) * pXY / 2.0;
                        } else {   // intermediate terms
                            tmpDens = stepT[0] * pXY;
                        }
                        long kk = getLinearIndex(i1, i2, 0, nX, nY);
                        mkde[kk] = mkde[kk] + tmpDens;
                    }
                }
                // update the eval time (t) and stopping conditions
                if (finalLoop) {
                    exitLoop = true;
                } else {
                    tOld = t;
                    t += stepT[0];
                    if (t >= t1) {
                        t = t1;
                        finalLoop = true;
                    }
                }
            }
        }
    }
    double maxDens = 0.0, sumDens = 0.0;
    long kk;
    for (int i1 = 0; i1 < nX; i1++) {
        for (int i2 = 0; i2 < nY; i2++) {
            kk = getLinearIndex(i1, i2, 0, nX, nY);
            mkde[kk] = mkde[kk]/totalT;
            if (mkde[kk] > maxDens) {
                maxDens = mkde[kk];
            }
            sumDens += mkde[kk];
        }
    }
    Rcpp::Rcout << "\tMaximum voxel density = " << maxDens << std::endl;
    Rcpp::Rcout << "\tSum of voxel densities = " << sumDens << std::endl;
    Rcpp::Rcout << "2D MKDE Computation: DONE" << std::endl;
    // RETURN DENSITY HERE
    return Rcpp::wrap(mkde);
}


/* this version computes the probability that the two individuals will be in the
   same cell, not at the exact same location.  To compute the probability of
   occurrence at the same point within a cell, we will have to integrate the
   product of the two kernels over the area of the cell */
RcppExport SEXP mkde2dGridv02interact(SEXP obsT, SEXP obsX0, SEXP obsY0,
                                      SEXP obsX1, SEXP obsY1, SEXP useObs,
                                      SEXP xgrid, SEXP ygrid,
                                      SEXP mvSig2xy0, SEXP mvSig2xy1,
                                      SEXP obsSig2xy0, SEXP obsSig2xy1,
                                      SEXP tStep, SEXP pdfThresh) {
    //
    // observations
    Rcpp::NumericVector T(obsT); // observed times
    Rcpp::NumericVector X0(obsX0); // observed x-coordinates
    Rcpp::NumericVector Y0(obsY0); // observed y-coordinates
    Rcpp::NumericVector X1(obsX1); // observed x-coordinates
    Rcpp::NumericVector Y1(obsY1); // observed y-coordinates
    Rcpp::IntegerVector isValid(useObs); // move step is flagged for use
    int nObs = T.length();
    // grid speces
    Rcpp::NumericVector xGrid(xgrid); // cell center coordinates in the x-dimension
    Rcpp::NumericVector yGrid(ygrid); // cell center coordinates in the y-dimension
    int nX = xGrid.length();
    int nY = yGrid.length();
    long nCells = (long)nX*nY;
    double xSz = xGrid[1] - xGrid[0];
    double ySz = yGrid[1] - yGrid[0];
    // variance parameters (vectors)
    Rcpp::NumericVector msig2xy0(mvSig2xy0); // move variance in the(x,y)-dimensions
    Rcpp::NumericVector osig2xy0(obsSig2xy0); // observation error variance in the(x,y)-dimensions
    Rcpp::NumericVector msig2xy1(mvSig2xy1); // move variance in the(x,y)-dimensions
    Rcpp::NumericVector osig2xy1(obsSig2xy1); // observation error variance in the(x,y)-dimensions
    // control parameters
    Rcpp::NumericVector stepT(tStep); // integration time step
    Rcpp::NumericVector pdfMin(pdfThresh); // integration time step
    // arrays for MKDE computation
    double * yprob0 = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double * yprob1 = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double * ybhattdist = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    std::vector<double> mkde(nCells);
    for (long i = 0; i < nCells; i++) {
        mkde[i] = 0.0;
    }
    // set up time variables
    double t0, t1, t, tOld, dt, alpha;
    // set up tmp variables
    double eX0, eY0, eX1, eY1;
    double W0 = 0.0, W1 = 0.0;
    // Save for possible later use double totalT; // T, Ttotal
    double sig2xy0, sig2xy1;
    double distMaxXY0, distMaxXY1, haloMinX, haloMaxX, haloMinY, haloMaxY, bhattaFact;
    int halo1min, halo1max, halo2min, halo2max;

    // FOLLOWING FOR DEBUGGING
    // int i1min = nX, i1max = 0, i2min = nY, i2max = 0;

    // start computing MKDE
    Rcpp::Rcout << "2D MKDE Interaction Computation: STARTING" << std::endl;
    // Save for possible later use totalT = 0.0;
    for (int j = 0; j < (nObs - 1); j++) {
        Rcpp::Rcout << "\tProcessing move step " << (j + 1) << " of " << (nObs - 1) << std::endl;
        // report percent complete after each observation
        t0 = T[j];
        t1 = T[j + 1];
        dt = t1 - t0;
        t = t0;
        if (isValid[j] == 1) {
	  // Save for possible later use totalT += dt;
            bool exitLoop = false;
            bool finalLoop = false;
            while (!exitLoop) { // iterate over integration time steps
                // Calculate fractional distance between t0 and current time
                alpha = (t - t0) / dt;
                // Calculate parameters for bilinear interpolation
                sig2xy0 = dt * alpha * (1.0 - alpha) * msig2xy0[j] +
                    osig2xy0[j] * (1.0 - alpha) * (1.0 - alpha) +
                    osig2xy0[j+1] * alpha * alpha;
                sig2xy1 = dt * alpha * (1.0 - alpha) * msig2xy1[j] +
                    osig2xy1[j] * (1.0 - alpha) * (1.0 - alpha) +
                    osig2xy1[j+1] * alpha * alpha;
                // Get (x,y,z) coordinates of kernel origin using linear interpolation
                eX0 = X0[j] + alpha * (X0[j + 1] - X0[j]);
                eY0 = Y0[j] + alpha * (Y0[j + 1] - Y0[j]);
                eX1 = X1[j] + alpha * (X1[j + 1] - X1[j]);
                eY1 = Y1[j] + alpha * (Y1[j + 1] - Y1[j]);
                // halo
                distMaxXY0 = univariateNormalDensityThreshold(pdfMin[0], sig2xy0); // ADD
                distMaxXY1 = univariateNormalDensityThreshold(pdfMin[0], sig2xy1); // ADD
                // x
                haloMinX = std::min(eX0 - distMaxXY0, eX1 - distMaxXY1); // ADD
                haloMaxX = std::max(eX0 + distMaxXY0, eX1 + distMaxXY1); // ADD
                halo1min = coordToIndex(haloMinX, xGrid[0], xSz); // ADD
                halo1max = coordToIndex(haloMaxX, xGrid[0], xSz); // ADD
                // y
                haloMinY = std::min(eY0 - distMaxXY0, eY1 - distMaxXY1); // ADD
                haloMaxY = std::max(eY0 + distMaxXY0, eY1 + distMaxXY1); // ADD
                halo2min = coordToIndex(haloMinY, yGrid[0], ySz); // ADD
                halo2max = coordToIndex(haloMaxY, yGrid[0], ySz); // ADD
                //
                /*
                Rcpp::Rcout << "\t\thalo1 = (" << halo1min << ", " << halo1max << 
                "), halo2 = (" << halo2min << ", " << halo2max << ")" << std::endl; // FOR TESTING
                */
                
                // Precompute voxel density in y dimension
                for(int i2 = 0; i2 < nY; i2++) {
                    yprob0[i2] =  integrateNormal(yGrid[i2] - 0.5*ySz, yGrid[i2] + 0.5*ySz, eY0, sqrt(sig2xy0));
                    yprob1[i2] =  integrateNormal(yGrid[i2] - 0.5*ySz, yGrid[i2] + 0.5*ySz, eY1, sqrt(sig2xy1));
                    ybhattdist[i2] = integrateKernelBC(yGrid[i2] - 0.5*ySz, yGrid[i2] + 0.5*ySz, eY0, sqrt(sig2xy0),eY1, sqrt(sig2xy1), pdfMin[0]);
                }
                bhattaFact = (RSQRT2PI * RSQRT2PI) / (sqrt(sig2xy0) * sqrt(sig2xy1));

                for (int i1 = std::max(0, halo1min); i1 < std::min(nX, halo1max); i1++) { // x-dimension
                    double voxx = xGrid[i1]; // voxel x
                    double xprob0 = integrateNormal(voxx - 0.5*xSz, voxx + 0.5*xSz, eX0, sqrt(sig2xy0));
                    double xprob1 = integrateNormal(voxx - 0.5*xSz, voxx + 0.5*xSz, eX1, sqrt(sig2xy1));
                    double xbhattdist = integrateKernelBC(voxx - 0.5*xSz, voxx + 0.5*xSz, eX0, sqrt(sig2xy0),eX1, sqrt(sig2xy1), pdfMin[0]);
                    for (int i2 = std::max(0, halo2min); i2 < std::min(nY, halo2max); i2++) { // y-dimension
                        // Calculate contribution of kernel to voxel
                        double pXY0 = xprob0*yprob0[i2];
                        double pXY1 = xprob1*yprob1[i2];
                        double pXY = xbhattdist*ybhattdist[i2]*bhattaFact;
                        // update (trapezoid rule)
                        double tmpDens, tmpDens0, tmpDens1;
                        if (doubleEquals(t, t0)) { // first term
                            tmpDens = stepT[0] * pXY / 2.0;
                            tmpDens0 = stepT[0] * pXY0 / 2.0;
                            tmpDens1 = stepT[0] * pXY1 / 2.0;
                        } else if (doubleEquals(t, t1)) { // last term
                            tmpDens = (t - tOld) * pXY / 2.0;
                            tmpDens0 = (t - tOld) * pXY0 / 2.0;
                            tmpDens1 = (t - tOld) * pXY1 / 2.0;
                        } else { // intermediate terms
                            tmpDens = stepT[0] * pXY;
                            tmpDens0 = stepT[0] * pXY0;
                            tmpDens1 = stepT[0] * pXY1;
                        }
                        // Add contribution to voxel (removed Kahan summation for now)
                        long kk = getLinearIndex(i1, i2, 0, nX, nY);
                        mkde[kk] += tmpDens;
                        W0 += tmpDens0;
                        W1 += tmpDens1;
                    }
                }
                // update the eval time (t) and stopping conditions
                if (finalLoop) {
                    exitLoop = true;
                } else {
                    tOld = t;
                    t += stepT[0];
                    if (t >= t1) {
                        t = t1;
                        finalLoop = true;
                    }
                }
            }
        }
    }
    // divide by totalT
    double maxDist = 0.0, sumDens = 0.0;
    double normConst = 1.0 / sqrt(W0*W1); // 1.0 / totalT
    long kk;
    for (int i1 = 0; i1 < nX; i1++) {
        for (int i2 = 0; i2 < nY; i2++) {
            kk = getLinearIndex(i1, i2, 0, nX, nY);
            mkde[kk] *= normConst;
            if (mkde[kk] > maxDist) {
                maxDist = mkde[kk];
            }
            sumDens += mkde[kk];
        }
    }
    Rcpp::Rcout << "\tNormalizing Constant = " << normConst << "(" << 1.0/normConst << ")" << std::endl;
    Rcpp::Rcout << "\tMaximum cell Bhattacharyya coeff = " << maxDist << std::endl;
    Rcpp::Rcout << "\tOverall Bhattacharyya coeff = " << sumDens << std::endl;
    Rcpp::Rcout << "2D MKDE Interaction Computation: DONE" << std::endl;
    // RETURN DENSITY HERE
    return Rcpp::wrap(mkde);
}


//
RcppExport SEXP mkde3dGridv02(SEXP obsT, SEXP obsX, SEXP obsY, SEXP obsZ, SEXP useObs,
                              SEXP xgrid, SEXP ygrid, SEXP zgrid,
                              SEXP zMinMatrix, SEXP zMaxMatrix,
                              SEXP mvSig2xy, SEXP mvSig2z, SEXP obsSig2, SEXP obsSig2z,
                              SEXP tStep, SEXP pdfThresh) {
    //
    // observations
    Rcpp::NumericVector T(obsT); // observed times
    Rcpp::NumericVector X(obsX); // observed x-coordinates
    Rcpp::NumericVector Y(obsY); // observed y-coordinates
    Rcpp::NumericVector Z(obsZ); // observed z-coordinates
    Rcpp::IntegerVector isValid(useObs); // move step is flagged for use
    int nObs = T.length();
    // grid speces
    Rcpp::NumericVector xGrid(xgrid); // cell center coordinates in the x-dimension
    Rcpp::NumericVector yGrid(ygrid); // cell center coordinates in the y-dimension
    Rcpp::NumericVector zGrid(zgrid); // cell center coordinates in the z-dimension
    int nX = xGrid.length();
    int nY = yGrid.length();
    int nZ = zGrid.length();
    long nVoxels = (long)nX*nY*nZ;
    double xSz = xGrid[1] - xGrid[0];
    double ySz = yGrid[1] - yGrid[0];
    double zSz = zGrid[1] - zGrid[0];
    // pysical constraints in z-dimension
    Rcpp::NumericMatrix zMin(zMinMatrix);  // lower z-limit at each (x,y)
    Rcpp::NumericMatrix zMax(zMaxMatrix);  // upper z-limit at each (x,y)
    // variance parameters (vectors)
    Rcpp::NumericVector msig2xy(mvSig2xy); // move variance in the(x,y)-dimensions
    Rcpp::NumericVector msig2z(mvSig2z);   // move variance in the z-dimension
    Rcpp::NumericVector osig2xy(obsSig2);  // observation error variance in the(x,y)-dimensions
    Rcpp::NumericVector osig2z(obsSig2z);  // observation error variance in the z-dimension
    // control parameters
    Rcpp::NumericVector stepT(tStep); // integration time step
    Rcpp::NumericVector pdfMin(pdfThresh); // integration time step
    // arrays for MKDE computation
    double * ydens = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double * zdens = (double *) malloc(nZ * sizeof(double)); // To Precompute z
    std::vector<double> mkde(nVoxels);
    for (long i = 0; i < nVoxels; i++) {
        mkde[i] = 0.0;
    }
    // set up time variables
    double t0, t1, t, tOld, dt, alpha;
    // set up tmp variables
    double eX, eY, eZ;
    double W;
    double sig2xy, sig2z;
    double distMaxXY, distMaxZ;
    int halo1, halo2, halo3, i1k, i2k, i3k;

    // FOLLOWING FOR DEBUGGING
    // int i1min = nX, i1max = 0, i2min = nY, i2max = 0, i3min = nZ, i3max = 0;

    // start computing MKDE
    Rcpp::Rcout << "3D MKDE Computation: STARTING" << std::endl;
    W = 0.0;
    for (int j = 0; j < (nObs - 1); j++) {
        Rcpp::Rcout << "\tProcessing move step " << (j + 1) << " of " << (nObs - 1) << std::endl;
        // report percent complete after each observation
        t0 = T[j];
        t1 = T[j + 1];
        dt = t1 - t0;
        t = t0;
        if (isValid[j] == 1) {
            bool exitLoop = false;
            bool finalLoop = false;
            while (!exitLoop) { // iterate over integration time steps
                // Calculate fractional distance between t0 and current time
                alpha = (t - t0) / dt;
                // Calculate parameters for bilinear interpolation
                sig2xy = dt * alpha * (1.0 - alpha) * msig2xy[j] +
                    osig2xy[j] * (1.0 - alpha) * (1.0 - alpha) +
                    osig2xy[j+1] * alpha * alpha;
                sig2z = dt * alpha * (1.0 - alpha) * msig2z[j] +
                    osig2z[j] * (1.0 - alpha) * (1.0 - alpha) +
                    osig2z[j+1] * alpha * alpha;
                // Get (x,y,z) coordinates of kernel origin using linear interpolation
                eX = X[j] + alpha * (X[j + 1] - X[j]);
                eY = Y[j] + alpha * (Y[j + 1] - Y[j]);
                eZ = Z[j] + alpha * (Z[j + 1] - Z[j]);
                // Convert to grid indices
                i1k = coordToIndex(eX, xGrid[0], xSz);
                i2k = coordToIndex(eY, yGrid[0], ySz);
                i3k = coordToIndex(eZ, zGrid[0], zSz);
                //
                distMaxXY = univariateNormalDensityThreshold(pdfMin[0], sig2xy);
                halo1 = ceil(distMaxXY/xSz);
                halo2 = ceil(distMaxXY/ySz);
                distMaxZ = univariateNormalDensityThreshold(pdfMin[0], sig2z);
                halo3 = ceil(distMaxZ/zSz);
                // Precompute voxel density in y dimension
                for(int i2 = 0; i2 < nY; i2++) {
                    ydens[i2] =  integrateNormal(yGrid[i2] - 0.5*ySz, yGrid[i2] + 0.5*ySz, eY, sqrt(sig2xy));
                }
                // Precompute voxel density in z dimension
                for(int i3 = 0; i3 < nZ; i3++) {
                    zdens[i3] = integrateNormal(zGrid[i3] - 0.5*zSz, zGrid[i3] + 0.5*zSz, eZ, sqrt(sig2z));
                }

                for (int i1 = std::max(0, i1k-halo1); i1 < std::min(nX, i1k+halo1); i1++) { // x-dimension
                    double voxx = xGrid[i1]; // voxel x
                    double xdens = integrateNormal(voxx - 0.5*xSz, voxx + 0.5*xSz, eX, sqrt(sig2xy));
                    for (int i2 = std::max(0, i2k-halo2); i2 < std::min(nY, i2k+halo2); i2++) { // y-dimension
                        // get the range of indexes and coordinates based on the physical boundaries
                        int i3lo = std::max(0, getLowerCellIndex(zMin(i1, i2), zGrid[0], zSz));
                        int i3hi = std::min(nZ, getUpperCellIndex(zMax(i1, i2), zGrid[0], zSz) + 1); // add 1 because less than
                        double loZ = indexToCellCenterCoord(i3lo, zGrid[0], zSz) - 0.5 * zSz;
                        double hiZ = indexToCellCenterCoord(i3hi - 1, zGrid[0], zSz) + 0.5 * zSz;
                        // Reflect E[Z] about the lower and upper boundaries
                        double loReflZ = 2.0 * loZ - eZ;
                        double hiReflZ = 2.0 * hiZ - eZ;
                        //
                        for (int i3 = std::max(i3lo, i3k-halo3); i3 < std::min(i3hi, i3k+halo3); i3++) { // z-dimension
                            //
                            // set up for reflection
                            // only compute if the expected location is within distance of boundary
                            double voxz = zGrid[i3];
                            double loDZ = 0.0;
                            if (std::fabs(eZ - loZ) <= distMaxZ) {
                                loDZ = integrateNormal(voxz - 0.5*zSz, voxz + 0.5*zSz, loReflZ, sqrt(sig2z));
                            }
                            double hiDZ = 0.0;
                            if (std::fabs(hiZ - eZ) <= distMaxZ) {
                                hiDZ = integrateNormal(voxz - 0.5*zSz, voxz + 0.5*zSz, hiReflZ, sqrt(sig2z));
                            }

                            // Calculate contribution of kernel to voxel
                            double pXYZ = xdens*ydens[i2]*(zdens[i3] + loDZ + hiDZ);
                            // update density (trapezoid rule)
                            double tmpDens;
                            if (doubleEquals(t, t0)) { // first term
                                tmpDens = stepT[0] * pXYZ / 2.0;
                            } else if (doubleEquals(t, t1)) { // last term
                                tmpDens = (t - tOld) * pXYZ / 2.0;
                            } else { // intermediate terms
                                tmpDens = stepT[0] * pXYZ;
                            }
                            // Add contribution to voxel (removed Kahan summation for now)
                            long kk = getLinearIndex(i1, i2, i3, nX, nY);
                            mkde[kk] += tmpDens;
                            W += tmpDens;
                        }
                    }
                }
                // end parallel for

                // update the eval time (t) and stopping conditions
                if (finalLoop) {
                    exitLoop = true;
                } else {
                    tOld = t;
                    t += stepT[0];
                    if (t >= t1) {
                        t = t1;
                        finalLoop = true;
                    }
                }
            }
        }
    }
    // divide by totalT
    double maxDens = 0.0, sumDens = 0.0;
    long kk;
    for (int i1 = 0; i1 < nX; i1++) {
        for (int i2 = 0; i2 < nY; i2++) {
            for (int i3 = 0; i3 < nZ; i3++) {
                kk = getLinearIndex(i1, i2, i3, nX, nY);
                mkde[kk] = mkde[kk] / W;
                if (mkde[kk] > maxDens) {
                    maxDens = mkde[kk];
                }
                sumDens += mkde[kk];
            }
        }
    }
    Rcpp::Rcout << "\tMaximum voxel density = " << maxDens << std::endl;
    Rcpp::Rcout << "\tSum of voxel densities = " << sumDens << std::endl;
    Rcpp::Rcout << "3D MKDE Computation: DONE" << std::endl;
    // RETURN DENSITY HERE
    return Rcpp::wrap(mkde);
}

/* this version computes the probability that the two individuals will be in the
   same voxel, not at the exact same location.  To compute the probability of
   occurrence at the same point within a voxel, we will have to integrate the
   product of the two kernels over the volume of the voxel */
RcppExport SEXP mkde3dGridv02interact(SEXP obsT, SEXP obsX0, SEXP obsY0, SEXP obsZ0,
                                      SEXP obsX1, SEXP obsY1, SEXP obsZ1, SEXP useObs,
                                      SEXP xgrid, SEXP ygrid, SEXP zgrid,
                                      SEXP zMinMatrix, SEXP zMaxMatrix,
                                      SEXP mvSig2xy0, SEXP mvSig2xy1, SEXP mvSig2z0,
                                      SEXP mvSig2z1, SEXP obsSig2xy0, SEXP obsSig2xy1,
                                      SEXP obsSig2z0, SEXP obsSig2z1,
                                      SEXP tStep, SEXP pdfThresh) {
    //
    // observations
    Rcpp::NumericVector T(obsT); // observed times
    Rcpp::NumericVector X0(obsX0); // observed x-coordinates
    Rcpp::NumericVector Y0(obsY0); // observed y-coordinates
    Rcpp::NumericVector Z0(obsZ0); // observed z-coordinates
    Rcpp::NumericVector X1(obsX1); // observed x-coordinates
    Rcpp::NumericVector Y1(obsY1); // observed y-coordinates
    Rcpp::NumericVector Z1(obsZ1); // observed z-coordinates
    Rcpp::IntegerVector isValid(useObs); // move step is flagged for use
    int nObs = T.length();
    // grid speces
    Rcpp::NumericVector xGrid(xgrid); // cell center coordinates in the x-dimension
    Rcpp::NumericVector yGrid(ygrid); // cell center coordinates in the y-dimension
    Rcpp::NumericVector zGrid(zgrid); // cell center coordinates in the z-dimension
    int nX = xGrid.length();
    int nY = yGrid.length();
    int nZ = zGrid.length();
    long nVoxels = (long)nX*nY*nZ;
    double xSz = xGrid[1] - xGrid[0];
    double ySz = yGrid[1] - yGrid[0];
    double zSz = zGrid[1] - zGrid[0];
    // pysical constraints in z-dimension
    Rcpp::NumericMatrix zMin(zMinMatrix); // lower z-limit at each (x,y)
    Rcpp::NumericMatrix zMax(zMaxMatrix); // upper z-limit at each (x,y)
    // variance parameters (vectors)
    Rcpp::NumericVector msig2xy0(mvSig2xy0); // move variance in the(x,y)-dimensions
    Rcpp::NumericVector msig2z0(mvSig2z0); // move variance in the z-dimension
    Rcpp::NumericVector osig2xy0(obsSig2xy0); // observation error variance in the(x,y)-dimensions
    Rcpp::NumericVector osig2z0(obsSig2z0); // observation error variance in the z-dimension
    Rcpp::NumericVector msig2xy1(mvSig2xy1); // move variance in the(x,y)-dimensions
    Rcpp::NumericVector msig2z1(mvSig2z1); // move variance in the z-dimension
    Rcpp::NumericVector osig2xy1(obsSig2xy1); // observation error variance in the(x,y)-dimensions
    Rcpp::NumericVector osig2z1(obsSig2z1); // observation error variance in the z-dimension
    // control parameters
    Rcpp::NumericVector stepT(tStep); // integration time step
    Rcpp::NumericVector pdfMin(pdfThresh); // integration time step
    // arrays for MKDE computation
    double * yprob0 = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double * yprob1 = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double * ybhattdist = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double * zprob0 = (double *) malloc(nZ * sizeof(double)); // To Precompute z
    double * zprob1 = (double *) malloc(nZ * sizeof(double)); // To Precompute z
    double * zbhattdist = (double *) malloc(nZ * sizeof(double)); // To PRECOMPUTE Z
    std::vector<double> mkde(nVoxels);
    for (long i = 0; i < nVoxels; i++) {
        mkde[i] = 0.0;
    }
    // set up time variables
    double t0, t1, t, tOld, dt, alpha;
    // set up tmp variables
    double eX0, eY0, eZ0, eX1, eY1, eZ1;
    double W0 = 0.0, W1 = 0.0;
    double sig2xy0, sig2z0, sig2xy1, sig2z1;
    double distMaxXY0, distMaxXY1, distMaxZ0, distMaxZ1, bhattaFact;
    double haloMinX, haloMaxX, haloMinY, haloMaxY, haloMinZ, haloMaxZ;
    int halo1min, halo1max, halo2min, halo2max, halo3min, halo3max; //, i1k, i2k, i3k;

    // FOLLOWING FOR DEBUGGING
    // int i1min = nX, i1max = 0, i2min = nY, i2max = 0, i3min = nZ, i3max = 0;

    // start computing MKDE
    Rcpp::Rcout << "3D MKDE Interaction Computation: STARTING" << std::endl;
    for (int j = 0; j < (nObs - 1); j++) {
        Rcpp::Rcout << "\tProcessing move step " << (j + 1) << " of " << (nObs - 1) << std::endl;
        // report percent complete after each observation
        t0 = T[j];
        t1 = T[j + 1];
        dt = t1 - t0;
        t = t0;
        if (isValid[j] == 1) {
            bool exitLoop = false;
            bool finalLoop = false;
            while (!exitLoop) { // iterate over integration time steps
                // Calculate fractional distance between t0 and current time
                alpha = (t - t0) / dt;
                // Calculate parameters for bilinear interpolation
                sig2xy0 = dt * alpha * (1.0 - alpha) * msig2xy0[j] +
                    osig2xy0[j] * (1.0 - alpha) * (1.0 - alpha) +
                    osig2xy0[j+1] * alpha * alpha;
                sig2z0 = dt * alpha * (1.0 - alpha) * msig2z0[j] +
                    osig2z0[j] * (1.0 - alpha) * (1.0 - alpha) +
                    osig2z0[j+1] * alpha * alpha;
                sig2xy1 = dt * alpha * (1.0 - alpha) * msig2xy1[j] +
                    osig2xy1[j] * (1.0 - alpha) * (1.0 - alpha) +
                    osig2xy1[j+1] * alpha * alpha;
                sig2z1 = dt * alpha * (1.0 - alpha) * msig2z1[j] +
                    osig2z1[j] * (1.0 - alpha) * (1.0 - alpha) +
                    osig2z1[j+1] * alpha * alpha;
                // Get (x,y,z) coordinates of kernel origin using linear interpolation
                eX0 = X0[j] + alpha * (X0[j + 1] - X0[j]);
                eY0 = Y0[j] + alpha * (Y0[j + 1] - Y0[j]);
                eZ0 = Z0[j] + alpha * (Z0[j + 1] - Z0[j]);
                eX1 = X1[j] + alpha * (X1[j + 1] - X1[j]);
                eY1 = Y1[j] + alpha * (Y1[j + 1] - Y1[j]);
                eZ1 = Z1[j] + alpha * (Z1[j + 1] - Z1[j]);
                // halo
                distMaxXY0 = univariateNormalDensityThreshold(pdfMin[0], sig2xy0);
                distMaxXY1 = univariateNormalDensityThreshold(pdfMin[0], sig2xy1);
                distMaxZ0 = univariateNormalDensityThreshold(pdfMin[0], sig2z0);
                distMaxZ1 = univariateNormalDensityThreshold(pdfMin[0], sig2z1);
                // x
                haloMinX = std::min(eX0 - distMaxXY0, eX1 - distMaxXY1);
                haloMaxX = std::max(eX0 + distMaxXY0, eX1 + distMaxXY1);
                halo1min = coordToIndex(haloMinX, xGrid[0], xSz);
                halo1max = coordToIndex(haloMaxX, xGrid[0], xSz);
                // y
                haloMinY = std::min(eY0 - distMaxXY0, eY1 - distMaxXY1);
                haloMaxY = std::max(eY0 + distMaxXY0, eY1 + distMaxXY1);
                halo2min = coordToIndex(haloMinY, yGrid[0], ySz);
                halo2max = coordToIndex(haloMaxY, yGrid[0], ySz);
                // z
                haloMinZ = std::min(eZ0 - distMaxZ0, eZ1 - distMaxZ1);
                haloMaxZ = std::max(eZ0 + distMaxZ0, eZ1 + distMaxZ1);
                halo3min = coordToIndex(haloMinZ, zGrid[0], zSz);
                halo3max = coordToIndex(haloMaxZ, zGrid[0], zSz);
                //
                /*
                Rcpp::Rcout << "\t\thalo1 = (" << halo1min << ", " << 
                  halo1max << "), halo2 = (" << halo2min << ", " << 
                  halo2max << "), halo3 = (" << halo3min << ", " << 
                  halo3max << ")" << std::endl; // FOR TESTING
                */
                
                // Precompute voxel density in y dimension
                for(int i2 = 0; i2 < nY; i2++) {
                    yprob0[i2] =  integrateNormal(yGrid[i2] - 0.5*ySz, yGrid[i2] + 0.5*ySz, eY0, sqrt(sig2xy0));
                    yprob1[i2] =  integrateNormal(yGrid[i2] - 0.5*ySz, yGrid[i2] + 0.5*ySz, eY1, sqrt(sig2xy1));
                    ybhattdist[i2] = integrateKernelBC(yGrid[i2] - 0.5*ySz, yGrid[i2] + 0.5*ySz, eY0, sqrt(sig2xy0),eY1, sqrt(sig2xy1), pdfMin[0]);
                }
                // Precompute voxel density in z dimension
                for(int i3 = 0; i3 < nZ; i3++) {
                    zprob0[i3] = integrateNormal(zGrid[i3] - 0.5*zSz, zGrid[i3] + 0.5*zSz, eZ0, sqrt(sig2z0));
                    zprob1[i3] = integrateNormal(zGrid[i3] - 0.5*zSz, zGrid[i3] + 0.5*zSz, eZ1, sqrt(sig2z1));
                    zbhattdist[i3] = integrateKernelBC(zGrid[i3] - 0.5*zSz, zGrid[i3] + 0.5*zSz, eZ0, sqrt(sig2z0),eZ1, sqrt(sig2z1), pdfMin[0]);
                }
                bhattaFact = (RSQRT2PI * RSQRT2PI * RSQRT2PI) / (sqrt(sig2xy0) * sqrt(sig2xy1) * sqrt(sqrt(sig2z0) * sqrt(sig2z1)));

                for (int i1 = std::max(0, halo1min); i1 < std::min(nX, halo1max); i1++) { // x-dimension
                    double voxx = xGrid[i1]; // voxel x
                    double xprob0 = integrateNormal(voxx - 0.5*xSz, voxx + 0.5*xSz, eX0, sqrt(sig2xy0));
                    double xprob1 = integrateNormal(voxx - 0.5*xSz, voxx + 0.5*xSz, eX1, sqrt(sig2xy1));
                    double xbhattdist = integrateKernelBC(voxx - 0.5*xSz, voxx + 0.5*xSz, eX0, sqrt(sig2xy0),eX1, sqrt(sig2xy1), pdfMin[0]);
                    for (int i2 = std::max(0, halo2min); i2 < std::min(nY, halo2max); i2++) { // y-dimension
                        // get the range of indexes and coordinates based on the physical boundaries
                        int i3lo = std::max(0, getLowerCellIndex(zMin(i1, i2), zGrid[0], zSz));
                        int i3hi = std::min(nZ, getUpperCellIndex(zMax(i1, i2), zGrid[0], zSz) + 1); // add 1 because less than
                        double loZ = indexToCellCenterCoord(i3lo, zGrid[0], zSz) - 0.5 * zSz;
                        double hiZ = indexToCellCenterCoord(i3hi - 1, zGrid[0], zSz) + 0.5 * zSz;
                        // Reflect E[Z] about the lower and upper boundaries
                        double loReflZ0 = 2.0 * loZ - eZ0;
                        double hiReflZ0 = 2.0 * hiZ - eZ0;
                        double loReflZ1 = 2.0 * loZ - eZ1;
                        double hiReflZ1 = 2.0 * hiZ - eZ1;
                        //
                        for (int i3 = std::max(i3lo, halo3min); i3 < std::min(i3hi, halo3max); i3++) { // z-dimension
                            // set up for reflection: HOW RELEVANT IS THE REFLECTION METHOD FOR INTERACTION?
                            // only compute if the expected location is within distance of boundary
                            double voxz = zGrid[i3];
                            double loDZ0 = 0.0, hiDZ0 = 0.0, loDZ1 = 0.0, hiDZ1 = 0.0, loBD = 0.0, hiBD = 0.0;
                            if ((std::fabs(eZ0 - loZ) <= distMaxZ0) || (std::fabs(eZ1 - loZ) <= distMaxZ1)) {
                                loDZ0 = integrateNormal(voxz - 0.5*zSz, voxz + 0.5*zSz, loReflZ0, sqrt(sig2z0));
                                loDZ1 = integrateNormal(voxz - 0.5*zSz, voxz + 0.5*zSz, loReflZ1, sqrt(sig2z1));
                                loBD = integrateKernelBC(voxz - 0.5*zSz, voxz + 0.5*zSz, loReflZ0, sqrt(sig2z0), loReflZ1, sqrt(sig2z1), pdfMin[0]);
                            }
                            if ((std::fabs(hiZ - eZ0) <= distMaxZ0) || (std::fabs(hiZ - eZ1) <= distMaxZ1)) {
                                hiDZ0 = integrateNormal(voxz - 0.5*zSz, voxz + 0.5*zSz, hiReflZ0, sqrt(sig2z0));
                                hiDZ1 = integrateNormal(voxz - 0.5*zSz, voxz + 0.5*zSz, hiReflZ1, sqrt(sig2z1));
                                hiBD = integrateKernelBC(voxz - 0.5*zSz, voxz + 0.5*zSz, hiReflZ0, sqrt(sig2z0), hiReflZ1, sqrt(sig2z1), pdfMin[0]);
                            }

                            // Calculate contribution of kernel to voxel
                            double pXYZ0 = xprob0*yprob0[i2]*(zprob0[i3] + loDZ0 + hiDZ0);
                            double pXYZ1 = xprob1*yprob1[i2]*(zprob1[i3] + loDZ1 + hiDZ1);
                            double pXYZ = xbhattdist*ybhattdist[i2]*(zbhattdist[i3] + loBD + hiBD) * bhattaFact;
                            // update density (trapezoid rule)
                            double tmpDens, tmpDens0, tmpDens1;
                            if (doubleEquals(t, t0)) { // first term
                                tmpDens0 = stepT[0] * pXYZ0 / 2.0;
                                tmpDens1 = stepT[0] * pXYZ1 / 2.0;
                                tmpDens = stepT[0] * pXYZ / 2.0;
                            } else if (doubleEquals(t, t1)) { // last term
                                tmpDens0 = (t - tOld) * pXYZ0 / 2.0;
                                tmpDens1 = (t - tOld) * pXYZ1 / 2.0;
                                tmpDens = (t - tOld) * pXYZ / 2.0;
                            } else { // intermediate terms
                                tmpDens0 = stepT[0] * pXYZ0;
                                tmpDens1 = stepT[0] * pXYZ1;
                                tmpDens = stepT[0] * pXYZ;
                            }
                            // Add contribution to voxel (removed Kahan summation for now)
                            long kk = getLinearIndex(i1, i2, i3, nX, nY);
                            mkde[kk] += tmpDens; // mkde[kk] + tmpDens;
                            W0 += tmpDens0;
                            W1 += tmpDens1;
                        }
                    }
                }
                // update the eval time (t) and stopping conditions
                if (finalLoop) {
                    exitLoop = true;
                } else {
                    tOld = t;
                    t += stepT[0];
                    if (t >= t1) {
                        t = t1;
                        finalLoop = true;
                    }
                }
            }
        }
        // Rcpp::Rcout << "\t\tW0 = " << W0 << ", W1 = " << W1 << std::endl;
    }
    // divide by totalT
    double maxDist = 0.0, sumDens = 0.0;
    double normConst = 1.0 / sqrt(W0*W1);
    long kk;
    for (int i1 = 0; i1 < nX; i1++) {
        for (int i2 = 0; i2 < nY; i2++) {
            for (int i3 = 0; i3 < nZ; i3++) {
                kk = getLinearIndex(i1, i2, i3, nX, nY);
                mkde[kk] = normConst * mkde[kk];
                if (mkde[kk] > maxDist) {
                    maxDist = mkde[kk];
                }
                sumDens += mkde[kk];
            }
        }
    }
    Rcpp::Rcout << "\tNormalizing Constant = " << normConst << "(" << 1.0/normConst << ")" << std::endl;
    Rcpp::Rcout << "\tMaximum voxel Bhattacharyya coeff = " << maxDist << std::endl;
    Rcpp::Rcout << "\tOverall Bhattacharyya coeff = " << sumDens << std::endl;
    Rcpp::Rcout << "3D MKDE Interaction Computation: DONE" << std::endl;
    // RETURN DENSITY HERE
    return Rcpp::wrap(mkde);
}

/* Compute Elevation Raster Cell Surface Area
   Based on algorithm by Jeff Jenness, modified to use bilinear
   interpolation in diagonal directions
*/
RcppExport SEXP computeCellSurfaceArea(SEXP elevGrid, SEXP cellSize) {
    Rcpp::NumericMatrix elev(elevGrid);
    Rcpp::NumericVector sz(cellSize);
    int nRows = elev.nrow();
    int nCols = elev.ncol();
    long totalCells = nRows * nCols;
    std::vector<double> area(totalCells);
    // Notes:
    //   elev arranged so that elev[0][0] is LL corner
    //   Triangle Corner Order: E, NE, N, NW, W, SW, S, SE
    double halfW = 0.5 * sz[0];
    double dx[] = {halfW, halfW, 0.0, -halfW, -halfW, -halfW, 0.0, halfW};
    double dy[] = {0.0, halfW, halfW, halfW, 0.0, -halfW, -halfW, -halfW};
    int nInterpPoints[] = {2, 4, 2, 4, 2, 4, 2, 4};
    int dr[8][4] = {
        {0, 0, 2, 2},
        {1, 1, 0, 0},
        {1, 0, 2, 2},
        {1, 1, 0, 0},
        {0, 0, 2, 2},
        {0, 0, -1, -1},
        {0, -1, 2, 2},
        {0, 0, -1, -1}
    }; // 2 means NA
    int dc[8][4] = {
        {0, 1, 2, 2},
        {0, 1, 0, 1},
        {0, 0, 2, 2},
        {-1, 0, -1, 0},
        {-1, 0, 2, 2},
        {-1, 0, -1, 0},
        {0, 0, 2, 2},
        {0, 1, 0, 1}
    }; // 2 means NA
    //
    long areaIndex = 0;
    int kk;
    double focalElev, x2, y2, z2, side1, side2, side3, abc;
    std::vector<double> pElev(8);
    //
    for (int r = 0; r < nRows; r++) {
        for (int c = 0; c < nCols; c++) {
            area[areaIndex] = 0.0;
            // check if cell is on edge
            if ((r > 0) && (r < (nRows - 1)) && (c > 0) && (c < (nCols - 1))) {
                // compute area of cell
                focalElev = elev(r, c); // focal cell elevation
                // interpolate elev at triangle corners on focal cell sides
                for (int k = 0; k < 8; k++) {
                    pElev[k] = 0.0;
                    for (int p = 0; p < nInterpPoints[k]; p++) {
                        pElev[k] += elev(r + dr[k][p], c + dc[k][p]);
                    }
                    pElev[k] /= ((double)nInterpPoints[k]);
                }
                // compute area
                for (int k = 0; k < 8; k++) {
                    kk = (k + 1)%8;
                    // focal to p0
                    x2 = dx[k] * dx[k];
                    y2 = dy[k] * dy[k];
                    z2 = (pElev[k] - focalElev) * (pElev[k] - focalElev);
                    side1 = sqrt(x2 + y2 + z2);
                    // p0 to p1
                    x2 = (dx[kk] - dx[k]) * (dx[kk] - dx[k]);
                    y2 = (dy[kk] - dy[k]) * (dy[kk] - dy[k]);
                    z2 = (pElev[kk] - pElev[k]) * (pElev[kk] - pElev[k]);
                    side2 = sqrt(x2 + y2 + z2);
                    // focal to p1
                    x2 = dx[kk] * dx[kk];
                    y2 = dy[kk] * dy[kk];
                    z2 = (pElev[kk] - focalElev) * (pElev[kk] - focalElev);
                    side3 = sqrt(x2 + y2 + z2);
                    //
                    abc = 0.5 * (side1 + side2 + side3);
                    area[areaIndex] += sqrt(abc * (abc - side1) * (abc - side2)*(abc - side3));
                }
            }
            areaIndex++;
        }
    }
    // return result (will have to be put into 2D area in R, filled by row)
    return Rcpp::wrap(area); //
}

RcppExport SEXP writeMKDE3DtoVTK(SEXP xgrid, SEXP ygrid, SEXP zgrid, SEXP density, SEXP filename, SEXP description) {
    Rcpp::NumericVector xGrid(xgrid); // cell centers in the x-dimension
    Rcpp::NumericVector yGrid(ygrid); // cell centers in the y-dimension
    Rcpp::NumericVector zGrid(zgrid); // cell centers in the z-dimension
    int nX = (long)xGrid.length();
    int nY = (long)yGrid.length();
    int nZ = (long)zGrid.length();
    long nPoints = (long)nX*nY*nZ, ijk = 0;
    std::vector<double> d = Rcpp::as<std::vector<double> >(density);
    std::string descr = Rcpp::as<std::string>(description);
    std::string fname = Rcpp::as<std::string>(filename);
    char * fnm = new char[fname.size()+1];
    fnm[fname.size()] = 0;
    memcpy(fnm, fname.c_str(), fname.size());


    double densTmp;

    std::ofstream vtkfile;
    vtkfile.open (fnm);
    // write header
    vtkfile << "# vtk DataFile Version 3.0" << std::endl;
    vtkfile << descr << std::endl;
    vtkfile << "ASCII" << std::endl;
    vtkfile << "DATASET STRUCTURED_GRID" << std::endl;
    vtkfile << "DIMENSIONS " << nX << " " << nY << " " << nZ << std::endl;
    // write points
    vtkfile << "POINTS " << nPoints << " float" << std::endl;
    vtkfile << std::scientific;
    for (int k = 0; k < nZ; k++) {
        for (int j = 0; j < nY; j++) {
            for (int i = 0; i < nX; i++) {
                vtkfile << xGrid[i] << " " << yGrid[j] << " " << zGrid[k] << std::endl;
            }
        }
    }
    // write data
    vtkfile << std::endl << "POINT_DATA " << nPoints << std::endl;
    vtkfile << "SCALARS density float 1" << std::endl;
    vtkfile << "LOOKUP_TABLE densityTable" << std::endl;
    for (int k = 0; k < nZ; k++) {
        for (int j = 0; j < nY; j++) {
            for (int i = 0; i < nX; i++) {
                ijk = getLinearIndex(i, j, k, nX, nY);
                densTmp = d[ijk];
                if (std::isnan(densTmp)) {
                    vtkfile << "0.0000000000000" << std::endl;
                } else {
                    // possibly use std::scientific or #include <iomanip> with std::setprecision(12)
                    vtkfile << densTmp << std::endl;
                }
            }
        }
    }
    // write lookup table
    vtkfile << std::endl << "LOOKUP_TABLE densityTable 9" << std::endl;
    vtkfile << "255 255 204 0.3" << std::endl;
    vtkfile << "255 237 160 0.4" << std::endl;
    vtkfile << "254 217 118 0.5" << std::endl;
    vtkfile << "254 178 76 0.6" << std::endl;
    vtkfile << "253 141 60 0.7" << std::endl;
    vtkfile << "252 78 42 0.8" << std::endl;
    vtkfile << "227 26 28 0.9" << std::endl;
    vtkfile << "189 0 38 0.9" << std::endl;
    vtkfile << "128 0 38 1.0" << std::endl;
    // done!
    vtkfile.close();
    return Rcpp::wrap(1);
}

// write to GRASS GIS 3D ASCII raster file
RcppExport SEXP writeMKDE3DtoGRASS(SEXP xgrid, SEXP ygrid, SEXP zgrid, SEXP density, SEXP filename, SEXP nodata) {
    Rcpp::NumericVector xGrid(xgrid); // cell centers in the x-dimension
    Rcpp::NumericVector yGrid(ygrid); // cell centers in the y-dimension
    Rcpp::NumericVector zGrid(zgrid); // cell centers in the z-dimension
    int nX = (long)xGrid.length();
    int nY = (long)yGrid.length();
    int nZ = (long)zGrid.length();
    long ijk = 0;
    double xSz = xGrid[1] - xGrid[0];
    double ySz = yGrid[1] - yGrid[0];
    double zSz = zGrid[1] - zGrid[0];
    double densTmp;
    std::vector<double> d = Rcpp::as<std::vector<double> >(density);
    std::string fname = Rcpp::as<std::string>(filename);
    char * fnm = new char[fname.size()+1];
    fnm[fname.size()] = 0;
    memcpy(fnm, fname.c_str(), fname.size());
    std::string nv = Rcpp::as<std::string>(nodata);

    // open file
    std::ofstream r3file;
    r3file.open (fnm);
    r3file << std::setprecision(12);

    // may need to set precision first
    // header
    r3file << "north: " << yGrid[nY - 1] + 0.5 * ySz << std::endl; // north: y.max
    r3file << "south: " << yGrid[0] - 0.5 * ySz << std::endl; // south: y.min
    r3file << "east: " << xGrid[nX - 1] + 0.5 * xSz << std::endl; // east: x.max
    r3file << "west: " << xGrid[0] - 0.5 * xSz << std::endl; // west: x.min
    r3file << "top: " << zGrid[nZ - 1] + 0.5 * zSz << std::endl; // top: z.max
    r3file << "bottom: " << zGrid[0] - 0.5 * zSz << std::endl;// bottom: z.min
    r3file << "rows: " << nY << std::endl; // rows
    r3file << "cols: " << nX << std::endl; // cols
    r3file << "levels: " << nZ << std::endl; // levels: integer

    // data
    // possibly use std::scientific or #include <iomanip> with std::setprecision(12)
    for (int k = 0; k < nZ; k++) {
        for (int j = 0; j < nY; j++) {
            for (int i = 0; i < nX; i++) {
                ijk = getLinearIndex(i, j, k, nX, nY);
                densTmp = d[ijk];
                if (std::isnan(densTmp)) {
                    r3file << nv;
                } else {
                    r3file << densTmp;
                }
                if (i == (nX - 1)) {
                    r3file << std::endl;
                } else {
                    r3file << " ";
                }
            }
        }
    }

    // done...
    r3file.close();
    return Rcpp::wrap(1);
}

// write to XDMF file
RcppExport SEXP writeMKDE3DtoXDMF(SEXP xgrid, SEXP ygrid, SEXP zgrid, SEXP density, SEXP filenameXDMF, SEXP filenameDAT) {
    Rcpp::NumericVector xGrid(xgrid); // cell centers in the x-dimension
    Rcpp::NumericVector yGrid(ygrid); // cell centers in the y-dimension
    Rcpp::NumericVector zGrid(zgrid); // cell centers in the z-dimension
    int nX = (long)xGrid.length();
    int nY = (long)yGrid.length();
    int nZ = (long)zGrid.length();
    long ijk = 0;
    double xSz = xGrid[1] - xGrid[0];
    double ySz = yGrid[1] - yGrid[0];
    double zSz = zGrid[1] - zGrid[0];
    double densTmp;
    std::vector<double> d = Rcpp::as<std::vector<double> >(density);
    std::string strXDMF = Rcpp::as<std::string>(filenameXDMF);
    std::string strDAT = Rcpp::as<std::string>(filenameDAT);
    char * fnmXDMF = new char[strXDMF.size()+1];
    fnmXDMF[strXDMF.size()] = 0;
    memcpy(fnmXDMF, strXDMF.c_str(), strXDMF.size());
    char * fnmDAT =new char[strDAT.size()+1];
    fnmDAT[strDAT.size()] = 0;
    memcpy(fnmDAT, strDAT.c_str(), strDAT.size());
    int nmSize = 0;
    // scan backward through char array, count chars, break when hit "/"; should count \0
    for (int i = strDAT.size(); i >= 0; i--) {
        if (fnmDAT[i] == '/') {
            break;
        } else {
            nmSize++;
        }
    }
    char * binName = new char[nmSize];
    int j = 0;
    for (int i = 0; i < nmSize; i++) {
        j = strDAT.size() + 1 - nmSize + i; // CHECK THIS!!!!
        binName[i] = fnmDAT[j];
    }
    // now copy name to string

    // write XML wrapper
    std::ofstream xmffile;
    xmffile.open (fnmXDMF);
    xmffile << std::setprecision(12);

    xmffile << "<?xml version=\"1.0\" ?>" << std::endl;
    xmffile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
    xmffile << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">" << std::endl;
    xmffile << "<Domain>" << std::endl;
    xmffile << "    <Grid Name=\"Mesh\" GridType=\"Uniform\">" << std::endl;
    xmffile << "        <Topology name=\"topo\" TopologyType=\"3DCoRectMesh\"" << std::endl;
    xmffile << "            Dimensions=\"" << (nZ + 1) << " " << (nY + 1) << " " << (nX + 1) << "\">" << std::endl; // levels, rows, cols?
    xmffile << "        </Topology>" << std::endl;
    xmffile << "        <Geometry name=\"geo\" Type=\"ORIGIN_DXDYDZ\">" << std::endl;
    xmffile << "            <!-- Origin -->" << std::endl;
    xmffile << "            <DataItem Format=\"XML\" Dimensions=\"3\">" << std::endl;
    xmffile << "             " << " " << (xGrid[0] - 0.5*xSz) << " " << (yGrid[0] - 0.5*ySz) << " " << (zGrid[0] - 0.5*zSz) << std::endl;
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "            <!-- DxDyDz -->" << std::endl;
    xmffile << "            <DataItem Format=\"XML\" Dimensions=\"3\">" << std::endl;
    xmffile << "             " << xSz << " " << ySz << " " << zSz <<  std::endl;
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "        </Geometry>" << std::endl;
    xmffile << "        <Attribute Name=\"Density\" Center=\"Cell\">" << std::endl; // need AttributeType="Scalar" or Type="Scalar" ?
    xmffile << "            <DataItem Format=\"Binary\"" << std::endl;
    xmffile << "             DataType=\"Double\"" << std::endl;
    xmffile << "             Precision=\"8\"" << std::endl;
    if (isMachineBigEndian()) {
        xmffile << "             Endian=\"Big\"" << std::endl;
    } else {
        xmffile << "             Endian=\"Little\"" << std::endl;
    }
    xmffile << "             Dimensions=\"" << nZ << " " << nY << " " << nX << "\">" << std::endl;
    xmffile << "               " << binName << std::endl;
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "        </Attribute>" << std::endl;
    xmffile << "    </Grid>" << std::endl;
    xmffile << "</Domain>" << std::endl;
    xmffile << "</Xdmf>" << std::endl;

    // close XML file
    xmffile.close();

    // write binary data (kji order)
    std::ofstream datfile(fnmDAT, std::ios::out | std::ios::trunc | std::ios::binary); // std::ios::out | std::ios::app | std::ios::binary
    if (!datfile.is_open()) {
        Rcpp::Rcout << "Error in writeMKDE3DtoXDMF(): Output file "<< fnmDAT << " could not be opened." << std::endl;
    } else {
        for (int k = 0; k < nZ; k++) {
            for (int j = 0; j < nY; j++) {
                for (int i = 0; i < nX; i++) {
                    ijk = getLinearIndex(i, j, k, nX, nY); // i, k, k, ...
                    densTmp = d[ijk];
                    if (std::isnan(densTmp)) {
                        densTmp = 0.0;
                    }
                    datfile.write((char *)(&densTmp), sizeof(densTmp));
                }
            }
        }
        datfile.close();
    }

    // done...
    return Rcpp::wrap(1);
}

RcppExport SEXP writeRasterToXDMF02(SEXP xgrid, SEXP ygrid, SEXP rast, SEXP filenameXDMF, SEXP filenameDAT) {
    Rcpp::NumericVector xGrid(xgrid); // cell centers in the x-dimension
    Rcpp::NumericVector yGrid(ygrid); // cell centers in the y-dimension
    int nX = (long)xGrid.length();
    int nY = (long)yGrid.length();
    double xSz = xGrid[1] - xGrid[0];
    double ySz = yGrid[1] - yGrid[0];
    double densTmp;
    std::vector<double> r = Rcpp::as<std::vector<double> >(rast);
    std::string strXDMF = Rcpp::as<std::string>(filenameXDMF);
    std::string strDAT = Rcpp::as<std::string>(filenameDAT);
    char * fnmXDMF = new char[strXDMF.size()+1];
    fnmXDMF[strXDMF.size()] = 0;
    memcpy(fnmXDMF, strXDMF.c_str(), strXDMF.size());
    char * fnmDAT =new char[strDAT.size()+1];
    fnmDAT[strDAT.size()] = 0;
    memcpy(fnmDAT, strDAT.c_str(), strDAT.size());
    int nmSize = 0;
    // scan backward through char array, count chars, break when hit "/"; should count \0
    for (int i = strDAT.size(); i >= 0; i--) {
        if (fnmDAT[i] == '/') {
            break;
        } else {
            nmSize++;
        }
    }
    char * binName = new char[nmSize];
    int j = 0;
    for (int i = 0; i < nmSize; i++) {
        j = strDAT.size() + 1 - nmSize + i; // CHECK THIS!!!!
        binName[i] = fnmDAT[j];
    }
    // now copy name to string

    // write XML wrapper
    std::ofstream xmffile;
    xmffile.open (fnmXDMF);
    xmffile << std::setprecision(12);

    xmffile << "<?xml version=\"1.0\" ?>" << std::endl;
    xmffile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
    xmffile << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">" << std::endl;
    xmffile << "<Domain>" << std::endl;
    xmffile << "    <Grid Name=\"Mesh\" GridType=\"Uniform\">" << std::endl;
    xmffile << "        <Topology name=\"topo\" TopologyType=\"2DCoRectMesh\"" << std::endl;
    xmffile << "            Dimensions=\"" << (nY + 1) << " " << (nX + 1) << "\">" << std::endl; // y, x or x, y ?
    xmffile << "        </Topology>" << std::endl;
    xmffile << "        <Geometry name=\"geo\" Type=\"ORIGIN_DXDY\">" << std::endl;
    xmffile << "            <!-- Origin -->" << std::endl;
    xmffile << "            <DataItem Format=\"XML\" Dimensions=\"2\">" << std::endl;
    xmffile << "             " << (xGrid[0] - 0.5*xSz) << " " << (yGrid[0] - 0.5*ySz) << std::endl; // y, x or x, y ?
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "            <!-- DxDy -->" << std::endl;
    xmffile << "            <DataItem Format=\"XML\" Dimensions=\"2\">" << std::endl;
    xmffile << "             " << xSz << " " << ySz <<  std::endl; // y, x or x, y ?
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "        </Geometry>" << std::endl;
    xmffile << "        <Attribute Name=\"Raster\" Center=\"Cell\">" << std::endl; // need AttributeType="Scalar" or Type="Scalar" ?
    xmffile << "            <DataItem Format=\"Binary\"" << std::endl;
    xmffile << "             DataType=\"Double\"" << std::endl;
    xmffile << "             Precision=\"8\"" << std::endl;
    if (isMachineBigEndian()) {
        xmffile << "             Endian=\"Big\"" << std::endl;
    } else {
        xmffile << "             Endian=\"Little\"" << std::endl;
    }
    xmffile << "             Dimensions=\"" << nY << " " << nX << "\">" << std::endl; // y, x or x, y ?
    xmffile << "               " << binName << std::endl;
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "        </Attribute>" << std::endl;
    xmffile << "    </Grid>" << std::endl;
    xmffile << "</Domain>" << std::endl;
    xmffile << "</Xdmf>" << std::endl;

    // close XML file
    xmffile.close();

    // write binary data (kji order)
    std::ofstream datfile(fnmDAT, std::ios::out | std::ios::trunc | std::ios::binary); // std::ios::out | std::ios::app | std::ios::binary
    if (!datfile.is_open()) {
        Rcpp::Rcout << "Error in writeMKDE3DtoXDMF(): Output file "<< fnmDAT << " could not be opened." << std::endl;
    } else {
        // this is effectively the same as the above
        for (int i = 0; i < r.size(); i++) {
            if (i%100000 == 0) {
                Rcpp::Rcout << "writing raster cell " << (i +1 ) << " of " << r.size() << " to file " << binName << std::endl;
            }
            densTmp = r[i];
            if (std::isnan(densTmp)) {
                densTmp = 0.0;
            }
            datfile.write((char *)(&densTmp), sizeof(densTmp));
        }
        datfile.close();
    }

    // done...
    return Rcpp::wrap(1);
}


RcppExport SEXP writeRasterToVTK02(SEXP xgrid, SEXP ygrid, SEXP elev, SEXP rd, SEXP gr, SEXP bl, SEXP description, SEXP filenameVTK) {
    Rcpp::NumericVector xGrid(xgrid); // cell centers in the x-dimension
    Rcpp::NumericVector yGrid(ygrid); // cell centers in the y-dimension
    int nX = (long)xGrid.length();
    int nY = (long)yGrid.length();
    std::string descr = Rcpp::as<std::string>(description);
    std::vector<double> r = Rcpp::as<std::vector<double> >(elev);
    std::vector<double> red = Rcpp::as<std::vector<double> >(rd);
    std::vector<double> grn = Rcpp::as<std::vector<double> >(gr);
    std::vector<double> blu = Rcpp::as<std::vector<double> >(bl);
    std::string strVTK = Rcpp::as<std::string>(filenameVTK);
    char * fnmVTK = new char[strVTK.size()+1];
    fnmVTK[strVTK.size()] = 0;
    memcpy(fnmVTK, strVTK.c_str(), strVTK.size());
    //
    long nPoints = nX*nY, ij, ul, ur, ll, lr;
    double tmpX, tmpY, tmpElev;

    // write VTK

    std::ofstream vtkfile;
    vtkfile.open (fnmVTK);
    // write header
    vtkfile << "# vtk DataFile Version 3.0" << std::endl;
    vtkfile << descr << std::endl;
    vtkfile << "ASCII" << std::endl;
    vtkfile << "DATASET POLYDATA" << std::endl;
    // write points
    Rcpp::Rcout << "Writing points to file " << fnmVTK << std::endl;
    vtkfile << "POINTS " << nPoints << " float" << std::endl;
    vtkfile << std::scientific;
    for (int i = 0; i < nY; i++) {
        for (int j = 0; j < nX; j++) {
            ij = i*nX + j; // i, k, k, ...
            tmpX = xGrid[j];
            tmpY = yGrid[i];
            tmpElev = r[ij];
            if (std::isnan(tmpElev)) {
                tmpElev = 0.0;
            }
            vtkfile << tmpX << " " << tmpY << " " << tmpElev << std::endl;
        }
    }
    // write data
    Rcpp::Rcout << "Writing polygons to file " << fnmVTK << std::endl;
    vtkfile << std::endl << "POLYGONS " << (nX - 1)*(nY - 1) << " " << (nX - 1)*(nY - 1)*5 << std::endl;
    for (int i = 0; i < (nY - 1); i++) {
        for (int j = 0; j < (nX - 1); j++) {
            ul = i*nX + j;
            ur = ul + 1;
            ll = (i + 1)*nX + j;
            lr = ll + 1;
            vtkfile << "4 " << ul << " " << ur << " " << lr << " " << ll << std::endl;
        }
    }
    // write lookup table
    Rcpp::Rcout << "Writing lookup table to file " << fnmVTK << std::endl;
    vtkfile << std::endl << "POINT_DATA " << nX*nY << std::endl;
    vtkfile << "SCALARS value float 1" << std::endl;
    vtkfile << "LOOKUP_TABLE default " << std::endl;
    for (int i = 0; i < nY; i++) {
        for (int j = 0; j < nX; j++) {
            ij = i*nX + j;
            tmpElev = r[ij];
            if (std::isnan(tmpElev)) {
                tmpElev = 0.0;
            }
            vtkfile << tmpElev << std::endl;
        }
    }

    // colors
    Rcpp::Rcout << "Writing colors to file " << fnmVTK << std::endl;
    vtkfile << std::endl << "COLOR_SCALARS RGB_Image 3" << std::endl;
    for (int i = 0; i < nY; i++) {
        for (int j = 0; j < nX; j++) {
            ij = i*nX + j;
            vtkfile << red[ij]/255.0 << " " << grn[ij]/255.0 << " " << blu[ij]/255.0 << std::endl;
        }
    }
    // done!
    vtkfile.close();
    return Rcpp::wrap(1);
}


