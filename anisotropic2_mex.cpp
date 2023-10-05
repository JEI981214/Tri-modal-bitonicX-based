#include "mex.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctype.h>
#include <vector>
#include <list>
using namespace std;

// Anisotropic Gaussian filter for 2D double images
void anisotropic_filter_2d(double *iarray, int w, int h, double f, double alpha, double *an, double *dir, double *co, double t, double *oarray)
{
  long i, j, x, y, filter, xmin, xmax, ymax, l, yi, li, xyi, iji, iw;
  double *theta, *ancosdir, *ansindir, *dsintheta, *dcostheta, *gauss, dx, dy, *d, *total, w1, w2, w3, *weight, ttotal, wweight, thresh;
  bool use_min_lambda = (co != NULL);

  // Calculate some constants
  thresh = t;
  if (alpha < 0.1) alpha = 0.1;
  if (thresh == 0.0) thresh = 1e10;
  filter = (long)floor(f);
  l = 2 * filter + 1;

  // Some additional arrays we need at this stage
  theta = new double[l*l];
  d = new double[l*l];
  dsintheta = new double[l*l];
  dcostheta = new double[l*l];
  gauss = new double[l*l];
  total = new double[w*h];
  weight = new double[w*h];
  ancosdir = new double[w*h];
  ansindir = new double[w*h];

  // precalculate some values
  for (i=0; i < l; i++) {
    for (j=0; j < l; j++) {
      dy = (double)(i - filter);
      dx = (double)(j - filter);
      d[i*l + j] = sqrt(dx*dx + dy * dy);
      theta[i*l + j] = atan2(dy, dx);
      if (theta[i*l + j] < 0.0) theta[i*l + j] += M_PI;
    }
  }

  // Transform some values for faster computation
  // after this:
  // theta = d x sin(theta)
  //     d = d x cos(theta)
  // ansindir = an x sin(dir)
  // ancosdir = an x cos(dir)
  // Also ensure non-gaussian for very short filters, otherwise they are not effective
  w3 = 0.65;
  for (i=l * filter; i < l*l; i++) {
    w2 = 2.0*d[i] / f;
    if ((d[i] > 0.0) && (filter == 1)) {
      w1 = 0.75 / d[i];
      gauss[i] = w1 + (1.0 - w1)*exp(-w2 * w2);
    } else if ((d[i] > 0.0) && (filter == 2)) {
      w1 = 0.25 / d[i];
      gauss[i] = w1 + (1.0 - w1)*exp(-w2 * w2);
    } else {
      gauss[i] = exp(-w2 * w2);
    }
    //d[i] = d[i] / (alpha * alpha);
    dsintheta[i] = d[i] * sin(theta[i]);
    dcostheta[i] = d[i] * cos(theta[i]);
    if (use_min_lambda) {
      d[i] = d[i] / w3;
    }
  }
  w1 = alpha * alpha;
  for (i=0; i < w*h; i++) {
    if ((an[i]*an[i]) > w1 * 1.2) {
      ansindir[i] = 1.2 * sin(dir[i]);
      ancosdir[i] = 1.2 * cos(dir[i]);
    } else {
      ansindir[i] = an[i] * an[i] * sin(dir[i]) / w1;
      ancosdir[i] = an[i] * an[i] * cos(dir[i]) / w1;
    }
    if (use_min_lambda) {
      if (co[i] > w3) {
        co[i] = w3;
      }
    }
    total[i] = 0.0;
    weight[i] = 0.0;
  }

  // Now filter actual image using these values
  // really want gauss / (1 + (an^2 / alpha^2) x d x sin(theta - dir)) OR
  // gauss / (1 + (an^2 / alpha^2) x ((d/2) x sin(theta - dir))^2) which is replaced using
  // trigonometric formulae with above pre-calculated values
  // Because the weights are symmetric, we can store them for both pixels being compared
  // which means we only have to search forwards: all the backwards calculations are done
  for (y=0; y < h; y++) {
    if (y > (h - filter - 1)) ymax = h - 1;
    else ymax = y + filter;
    for (x=0; x < w; x++) {
      xyi = y * w + x;
      if (x < filter) xmin = 0;
      else xmin = x - filter;
      if (x > (w - filter - 1)) xmax = w - 1;
      else xmax = x + filter;
      for (i=y; i <= ymax; i++) {
        yi = (i - y + filter)*l - x + filter;
        iw = i * w;
        if (i == y) {
          ttotal = total[xyi] + iarray[xyi];
          wweight = weight[xyi] + 1.0;
          for (j=(x + 1); j <= xmax; j++) {
            iji = iw + j;
            if (fabs(iarray[xyi] - iarray[iji]) < thresh) {
              li = yi + j;
              if (use_min_lambda) {
                w1 = 1.0 + d[li] * co[xyi] + fabs(dsintheta[li] * ancosdir[xyi] - dcostheta[li] * ansindir[xyi]);
                w1 *= (1.0 + d[li] * co[iji] + fabs(dsintheta[li] * ancosdir[iji] - dcostheta[li] * ansindir[iji]));
              } else {
                w1 = 1.0 + fabs(dsintheta[li] * ancosdir[xyi] - dcostheta[li] * ansindir[xyi]);
                w1 *= (1.0 + fabs(dsintheta[li] * ancosdir[iji] - dcostheta[li] * ansindir[iji]));
              }
              w1 = gauss[li] / w1;
              ttotal += w1 * iarray[iji];
              wweight += w1;
              total[iji] += w1 * iarray[xyi];
              weight[iji] += w1;
            }
          }
        } else {
          for (j=xmin; j <= xmax; j++) {
            iji = iw + j;
            if (fabs(iarray[xyi] - iarray[iji]) < thresh) {
              li = yi + j;
              if (use_min_lambda) {
                w1 = 1.0 + d[li] * co[xyi] + fabs(dsintheta[li] * ancosdir[xyi] - dcostheta[li] * ansindir[xyi]);
                w1 *= (1.0 + d[li] * co[iji] + fabs(dsintheta[li] * ancosdir[iji] - dcostheta[li] * ansindir[iji]));
              } else {
                w1 = 1.0 + fabs(dsintheta[li] * ancosdir[xyi] - dcostheta[li] * ansindir[xyi]);
                w1 *= (1.0 + fabs(dsintheta[li] * ancosdir[iji] - dcostheta[li] * ansindir[iji]));
              }
              w1 = gauss[li] / w1;
              ttotal += w1 * iarray[iji];
              wweight += w1;
              total[iji] += w1 * iarray[xyi];
              weight[iji] += w1;
            }
          }
        }
      }
      total[xyi] = ttotal;
      weight[xyi] = wweight;
    }
  }
  for (i=0; i < w*h; i++) {
    if (weight[i] > 0.0) {
      oarray[i] = total[i] / weight[i];
    } else {
      oarray[i] = 0.0;
    }
  }

  // clear up allocated arrays
  delete[] theta;
  delete[] d;
  delete[] gauss;
  delete[] total;
  delete[] weight;
  delete[] ancosdir;
  delete[] ansindir;
  delete[] dsintheta;
  delete[] dcostheta;
}

// Call this function using [X] = anisotropic2_mex(A, F, ALPHA, T, ANIS, DIR, CO)
#define A_IN prhs[0]
#define F_IN prhs[1]
#define ALPHA_IN prhs[2]
#define T_IN prhs[3]
#define ANIS_IN prhs[4]
#define DIR_IN prhs[5]
#define CO_IN prhs[6]
#define X_OUT plhs[0]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *A, *Anis, *Dir, f, alpha, t, *X, *Co;
  bool do_corner = false;
  int w, h;
  
  // check we have consistent input and output arguments
  if (nrhs<6) {
    mexErrMsgTxt("A, F, ALPHA, ANIS, DIR and T must all be defined.");
  }
  if (nrhs > 6) {
    do_corner = true;
  }
  if (nlhs>1) {
    mexErrMsgTxt("Only one output argument can be defined.");
  }

  // check input data is of the right type
  if (mxIsComplex(A_IN) || mxIsSparse(A_IN) || (mxGetClassID(A_IN) != mxDOUBLE_CLASS) || (mxGetNumberOfDimensions(A_IN) != 2)) {
    mexErrMsgTxt("Input data must be a 2D matrix of double real values.");
  }
  if (mxIsComplex(F_IN) || (mxGetNumberOfElements(F_IN) != 1)) {
    mexErrMsgTxt("F must be a scalar.");
  } else {
    f = mxGetScalar(F_IN);
  }
  if (mxIsComplex(ALPHA_IN) || (mxGetNumberOfElements(ALPHA_IN) != 1)) {
    mexErrMsgTxt("ALPHA must be a scalar.");
  } else {
    alpha = mxGetScalar(ALPHA_IN);
  }
  if (mxIsComplex(T_IN) || (mxGetNumberOfElements(T_IN) != 1)) {
    mexErrMsgTxt("T must be a scalar.");
  } else {
    t = mxGetScalar(T_IN);
  }
  if (mxIsComplex(ANIS_IN) || mxIsSparse(ANIS_IN) || (mxGetClassID(ANIS_IN) != mxDOUBLE_CLASS) || (mxGetNumberOfDimensions(ANIS_IN) != 2)
    || (mxGetM(ANIS_IN) != mxGetM(A_IN)) || (mxGetN(ANIS_IN) != mxGetN(A_IN)) ) {
    mexErrMsgTxt("ANIS must be a 2D matrix of double real values, the same size as A.");
  }
  if (mxIsComplex(DIR_IN) || mxIsSparse(DIR_IN) || (mxGetClassID(DIR_IN) != mxDOUBLE_CLASS) || (mxGetNumberOfDimensions(DIR_IN) != 2)
    || (mxGetM(DIR_IN) != mxGetM(A_IN)) || (mxGetN(DIR_IN) != mxGetN(A_IN))) {
    mexErrMsgTxt("DIR must be a 2D matrix of double real values, the same size as A.");
  }
  if (do_corner) {
    if (mxIsComplex(CO_IN) || mxIsSparse(CO_IN) || (mxGetClassID(CO_IN) != mxDOUBLE_CLASS) || (mxGetNumberOfDimensions(CO_IN) != 2)
        || (mxGetM(CO_IN) != mxGetM(A_IN)) || (mxGetN(CO_IN) != mxGetN(A_IN))) {
      mexErrMsgTxt("MINE must be a 2D matrix of double real values, the same size as A.");
    }
  }
  
  // Get input data and domain
  w = (int)mxGetM(A_IN);
  h = (int)mxGetN(A_IN);

  // Get input and output data
  A = (double *)mxGetData(A_IN);
  Anis = (double *)mxGetData(ANIS_IN);
  Dir = (double *)mxGetData(DIR_IN);
  if (do_corner) {
    Co = (double *)mxGetData(CO_IN);
  } else {
    Co = NULL;
  }
  X_OUT = mxCreateDoubleMatrix(w, h, mxREAL);
  X = (double *)mxGetData(X_OUT);

  // Finally run actual filter on this data
  anisotropic_filter_2d(A, w, h, f, alpha, Anis, Dir, Co, t, X);
}
