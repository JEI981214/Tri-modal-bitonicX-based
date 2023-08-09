#include "mex.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctype.h>
#include <vector>
#include <list>
#include <algorithm>
using namespace std;

typedef enum { LEFT=0, RIGHT=1, UP=2, DOWN=3 } dir_t;

// Integer version using histograms
void xrankopen_2d(int *iarray, int w, int h, int c, int f, int t, double g1, int mlevel, int *oarray1, int *oarray2, int *imarray, int *igarray)
{
  long i, j, k, cc, p, mp, ii, jj, iii, jjj, l, l2;
  int m, e, ig1, ig1r, ig2, ig2r, maxval, minval, mii, mii2, v, vrgb, v1[4], v2[4];
  int mtotal, *mmask, *mask, bins, b, bstart;
  int *igrey, *oa1, *oa2, *pa, *forward1, *forward2;
  int mine[4], maxe[4], minv, maxv, thresh, mthresh, g, ming, maxg, minrgb[4], maxrgb[4];
  int *harray[4];
  double dm, dm2;
  dir_t dir;
  bool inverse = true;
  bool reduce_memory = true;
  bool use_raw_rgb_threshold = (c > 1) && !((mlevel == 0) || (mlevel == -1));
  bool update_whole_mask, use_local_threshold; // set once we know max and min values

  // return if incorrect inputs
  if ((f < 1) || (w < 1) || (h < 1) || (c < 1) || (c > 4) || (g1 < 0.0) || (g1 > 0.5)) return;
  if ((oarray1 == NULL) && (oarray2 == NULL)) return;
  if ((imarray == NULL) || ((c > 1) && (igarray == NULL))) return;

  // Just copy image across if it is too small
  if ((f > (h / 2)) || (f > (w / 2))) {
    if (oarray1 != NULL) {
      for (i = 0; i < w * h * c; i++) {
        oarray1[i] = iarray[i];
      }
    }
    if (oarray2 != NULL) {
      for (i = 0; i < w * h * c; i++) {
        oarray2[i] = iarray[i];
      }
    }
    return;
  }

  // Calculate mask and bin sizes, and centile location
  l = 2 * f + 1;
  l2 = l * l;
  mp = l2;
  ig1 = (int)floor(g1 * (double)(mp - 1) + 0.5);
  if ((ig1 == 0) && (mp > 3) && (g1 > 0.0)) ig1 = 1;
  bins = ig1 + 1;
  if ((g1 == 0.5) && (oarray2 == NULL)) {
    // special case for median
    inverse = false;
  } else {
    // can do the whole process in one pass
    inverse = true;
  }

  // set up starting and ending points based on rectangular mask
  mmask = new int[l2];
  mask = new int[l2];
  for (i = 0; i < l2; i++) {
    mask[i] = 1;
    mmask[i] = 1;
  }

  // Check for colour input, and adjust input image and threshold
  if (c > 1) {
    igrey = igarray;
  } else {
    igrey = iarray;
  }

  // Compensation due to prior filtering, which reduces the noise
  if (mlevel == 0) {
    dm = sqrt(M_PI) * ((double)f + 0.5) * 1.2 * 0.85;
    // Further compensation to reduce the proportion of pixels being used as the mask gets larger
    dm2 = dm * erf(((double)f + 0.5) * 1.2 / 24.0);
  } else {
    // After the first level, the noise is correlated so filtering reduces it less
    // The filtering is twice as long as the previous level, so the effective filter length is about 2.0
    dm = 2.0;
    dm2 = dm * erf(((double)f + 0.5) * 1.2 / 12.0);
  }
  thresh = t;
  if (c > 1) thresh /= sqrt(3.0);
  mthresh = (int)((double)thresh / dm2 + 0.5);

  // check for input array range, use one for all components
  minval = iarray[0];
  maxval = minval;
  for (i = 1; i < w * h * c; i++) {
    if (iarray[i] > maxval) maxval = iarray[i];
    if (iarray[i] < minval) minval = iarray[i];
  }

  // use range to work out if we are in a low noise scenario
  update_whole_mask = (t <= (maxval-minval)/20); // more efficient for low noise scenario
  use_local_threshold = (t > (maxval-minval)/20) && (mlevel == 0);

  // initialise output if doing inverse in same loop as forward
  oa1 = oarray1;
  oa2 = oarray2;
  if (inverse) {
    if (oa1 != NULL) forward1 = new int[w * h * c];
    if (oa2 != NULL) forward2 = new int[w * h * c];
    pa = new int[w * h];
    for (i = 0; i < w * h; i++) pa[i] = 0;
    if (bins == 1) reduce_memory = false;
    if (reduce_memory) {
      // this only stores w * l rather than w * h
      if (oarray1 != NULL) oa1 = new int[w * l * c * bins];
      if (oarray2 != NULL) oa2 = new int[w * l * c * bins];
      if (oa1 != NULL) {
        for (i = 0; i < w * l * c * bins; i++) oa1[i] = minval - 1;
      }
      if (oa2 != NULL) {
        for (i = 0; i < w * l * c * bins; i++) oa2[i] = maxval + 1;
      }
    } else {
      if (bins > 1) {
        if (oarray1 != NULL) oa1 = new int[w * h * c * bins];
        if (oarray2 != NULL) oa2 = new int[w * h * c * bins];
      }
      if (oa1 != NULL) {
        for (i = 0; i < w * h * c * bins; i++) oa1[i] = minval - 1;
      }
      if (oa2 != NULL) {
        for (i = 0; i < w * h * c * bins; i++) oa2[i] = maxval + 1;
      }
    }
  }

  // initialise histogram at first pixel
  mp = 0;
  for (cc = 0; cc < c; cc++) {

    // create space for histogram, and fill with centre at (0, -1)
    e = maxval - minval + 1;
    harray[cc] = new int[e];
    for (m = 0; m < e; m++) harray[cc][m] = 0;
    if (!update_whole_mask) {
      for (i = -f - 1; i < f; i++) {
        if (i >= 0) {
          for (j = -f; j <= f; j++) {
            if (j >= 0) {
              m = iarray[(i * w + j) * c + cc] - minval;
              harray[cc][m]++;
              if (cc == 0) mp++;
            }
          }
        }
      }
    }
    mine[cc] = 0;
    maxe[cc] = e - 1;
  }

  // start walking through image using forwards then backwards approach to make the addition
  // of extra values to the sorted array more efficient
  i = 0;
  j = 0;
  dir = DOWN;
  while ((i < h) && (j < w)) {

    // Work out the list of old and new values we need to correct the sorted array with
    if (!update_whole_mask) {
      switch (dir) {

      case DOWN:
        for (k = 0; k < l; k++) {
          jj = (j + k - f);
          if ((jj >= 0) && (jj < w)) {
            ii = (i + f);
            if ((ii >= 0) && (ii < h)) {
              ii = (ii * w + jj) * c;
              mp++;
              for (cc = 0; cc < c; cc++) {
                m = iarray[ii + cc] - minval;
                harray[cc][m]++;
                if (m > maxe[cc]) maxe[cc] = m;
                if (m < mine[cc]) mine[cc] = m;
              }
            }
            iii = (i - f - 1);
            if ((iii >= 0) && (iii < h)) {
              iii = (iii * w + jj) * c;
              mp--;
              for (cc = 0; cc < c; cc++) {
                m = iarray[iii + cc] - minval;
                harray[cc][m]--;
              }
            }
          }
        }
        break;

      case UP:
        for (k = 0; k < l; k++) {
          jj = j + k - f;
          if ((jj >= 0) && (jj < w)) {
            ii = (i - f);
            if ((ii >= 0) && (ii < h)) {
              ii = (ii * w + jj) * c;
              mp++;
              for (cc = 0; cc < c; cc++) {
                m = iarray[ii + cc] - minval;
                harray[cc][m]++;
                if (m > maxe[cc]) maxe[cc] = m;
                if (m < mine[cc]) mine[cc] = m;
              }
            }
            iii = (i + f + 1);
            if ((iii >= 0) && (iii < h)) {
              iii = (iii * w + jj) * c;
              mp--;
              for (cc = 0; cc < c; cc++) {
                m = iarray[iii + cc] - minval;
                harray[cc][m]--;
              }
            }
          }
        }
        break;

      case RIGHT:
        for (k = 0; k < l; k++) {
          ii = i + k - f;
          if ((ii >= 0) && (ii < h)) {
            jj = j + f;
            if ((jj >= 0) && (jj < w)) {
              jj = (ii * w + jj) * c;
              mp++;
              for (cc = 0; cc < c; cc++) {
                m = iarray[jj + cc] - minval;
                harray[cc][m]++;
                if (m > maxe[cc]) maxe[cc] = m;
                if (m < mine[cc]) mine[cc] = m;
              }
            }
            jjj = j - f - 1;
            if ((jjj >= 0) && (jjj < w)) {
              jjj = (ii * w + jjj) * c;
              mp--;
              for (cc = 0; cc < c; cc++) {
                m = iarray[jjj + cc] - minval;
                harray[cc][m]--;
              }
            }
          }
        }
        break;

      case LEFT:
        for (k = 0; k < l; k++) {
          ii = i + k - f;
          if ((ii >= 0) && (ii < h)) {
            jj = j - f;
            if ((jj >= 0) && (jj < w)) {
              jj = (ii * w + jj) * c;
              mp++;
              for (cc = 0; cc < c; cc++) {
                m = iarray[jj + cc] - minval;
                harray[cc][m]++;
                if (m > maxe[cc]) maxe[cc] = m;
                if (m < mine[cc]) mine[cc] = m;
              }
            }
            jjj = j + f + 1;
            if ((jjj >= 0) && (jjj < w)) {
              jjj = (ii * w + jjj) * c;
              mp--;
              for (cc = 0; cc < c; cc++) {
                m = iarray[jjj + cc] - minval;
                harray[cc][m]--;
              }
            }
          }
        }
        break;
      }

      // Might be able to shrink max and min of histogram
      for (cc = 0; cc < c; cc++) {
        while (harray[cc][mine[cc]] == 0) mine[cc]++;
        while (harray[cc][maxe[cc]] == 0) maxe[cc]--;
      }
    }

    // Construct threshold mask for this location
    int iimax, iimin, jjmax, jjmin;
    iimin = i - f;
    if (iimin < 0) iimin = 0;
    iimax = i + f;
    if (iimax >= h) iimax = h - 1;
    jjmin = j - f;
    if (jjmin < 0) jjmin = 0;
    jjmax = j + f;
    if (jjmax >= w) jjmax = w - 1;
    p = mp;

    // Check distance of current value from min and max in histogram
    double tfact = 1.0;

    // This only happens in the first full size frame if doing multiresolution
    if (use_local_threshold) {
      double mrange, mminv;
      v = igrey[i * w + j];
      if ((c > 1) || update_whole_mask) {
        m = (f - i) * l + (f - j);
        maxg = v;
        ming = v;
        for (ii = iimin; ii <= iimax; ii++) {
          for (jj = jjmin; jj <= jjmax; jj++) {
            if (mask[ii * l + m + jj] > 0) { // in overall mask shape
              mii = ii * w + jj;
              if (igrey[mii] > maxg) maxg = igrey[mii];
              if (igrey[mii] < ming) ming = igrey[mii];
            }
          }
        }
      } else {
        ming = mine[0];
        maxg = maxe[0];
        ming = ming + minval;
        maxg = maxg + minval;
      }
      mrange = (double)(maxg - ming);
      mminv = (double)abs(v - ming);
      if (abs(maxg - v) < mminv) mminv = (double)abs(maxg - v);
      mrange = mrange / (double)thresh;
      mminv = mminv / (double)thresh;
      // These limits presume that 'thresh' has been set to four times the std of the noise.
      double low_filter, high_filter, range_limit, low_range_limit;
      low_filter = 0.7;
      high_filter = 1.6;
      range_limit = 1.57;
      low_range_limit = 0.9;
      if (mrange > range_limit) {
        // Reduce the threshold limits if the data range is very high
        // (and therefore the noise range is a small proportion of this)
        tfact = low_filter;
      } else if (mrange < low_range_limit) {
        // Increase the threshold limits if we are pretty sure that the noise occupies the whole range
        tfact = high_filter;
      } else {
        // Otherwise a blend between the two scenarios above
        tfact = high_filter - (mrange - low_range_limit) * (high_filter - low_filter) / (range_limit - low_range_limit);
      }
      if ((tfact < 1.0) && (mminv < low_filter)) {
        // As long as we are in the relatively high range case, also increase threshold if we are really close to the range edges
        // (since then we know the distribution is skewed towards other direction)
        tfact = 1.0 - mminv * (1.0 - tfact) / low_filter;
      }
    }

    // compensate for higher chance to de-select mask pixel because considering several channels
    // also, this factor isn't variable like tfact
    double tfactrgb = sqrt((double)c)*sqrt(3.0);
    if (tfact > 1.0) tfactrgb *= tfact;

    mii = i * w + j;

    // Work out what range of pixel values to include from (probably filtered) grey mask image imarray
    v = imarray[mii];
    minv = v - (int)(mthresh * tfact + 0.5);
    maxv = v + (int)(mthresh * tfact + 0.5);

    // Work out what range of pixel values to include from (un-filtered) grey image
    g = igrey[mii];
    ming = g - (int)(thresh * tfact + 0.5);
    maxg = g + (int)(thresh * tfact + 0.5);

    // Use raw threshold / average comparison to adjust mthresh
    // which allows us to ignore filtered image if it isn't from the same distribution as un-filtered data
    // the 0.5 allows for half the threshold range in the un-filtered data,
    // the 0.5/dm allows for additional half threshold range in the filtered data
    int test = (int)((double)thresh * (0.5 + 0.5 / dm) + 0.5);
    int vg = abs(v - g);
    if (vg > test) {
      minv -= 2 * (vg - test);
      maxv += 2 * (vg - test);
    }

    if (use_raw_rgb_threshold) {
      for (cc = 0; cc < c; cc++) {
        vrgb = iarray[mii * c + cc];
        minrgb[cc] = vrgb - (int)(thresh * tfactrgb + 0.5);
        maxrgb[cc] = vrgb + (int)(thresh * tfactrgb + 0.5);
      }
    }

    // Threshold based on imarray (and possibly also imarray2 and igrey)
    for (m = 0; m < l2; m++) mmask[m] = 0;
    mmask[(l2 - 1) / 2] = 1; // always include central pixel
    p = 1;
    m = (f - i) * l + (f - j);
    for (ii = iimin; ii <= iimax; ii++) {
      for (jj = jjmin; jj <= jjmax; jj++) {
        mii2 = ii * l + m + jj;
        if (2 * mii2 == (l2 - 1)) continue;
        if (mask[mii2] > 0) { // in overall mask shape
          mii = ii * w + jj;
          if ((imarray[mii] >= minv) && (imarray[mii] <= maxv)) { // within threshold of imarray
            if ((igrey[mii] >= ming) && (igrey[mii] <= maxg)) { // within threshold of igrey
              mmask[mii2] = 1;
              p++;
            }
          }
          if (use_raw_rgb_threshold && (mmask[mii2] == 1)) {
            for (cc = 0; cc < c; cc++) {
              vrgb = iarray[mii * c + cc];
              if ((vrgb < minrgb[cc]) || (vrgb > maxrgb[cc])) {
                mmask[mii2] = 0;
                p--;
                break;
              }
            }
          }
          if (!update_whole_mask) {
            if (mmask[mii2] <= 0) {
              // Adjust histograms based on new mask
              for (cc = 0; cc < c; cc++) {
                v = iarray[mii * c + cc] - minval;
                harray[cc][v]--;
              }
            }
          }
        }
      }
    }

    // Update current histogram or sorted list from mmask
    if (update_whole_mask) {
      bool first_value = true;
      int m2 = 0;
      m = (f - i) * l + (f - j);
      for (ii = iimin; ii <= iimax; ii++) {
        for (jj = jjmin; jj <= jjmax; jj++) {
          mii2 = ii * l + m + jj;
          if (mmask[mii2] > 0) {
            mii = ii * w + jj;

            // update histogram
            for (cc = 0; cc < c; cc++) {
              v = iarray[mii * c + cc] - minval;
              if (first_value) {
                harray[cc][v] = 1;
                mine[cc] = v;
                maxe[cc] = v;
              } else {
                if (v < mine[cc]) {
                  harray[cc][v] = 1;
                  for (m2 = (v + 1); m2 < mine[cc]; m2++) harray[cc][m2] = 0;
                  mine[cc] = v;
                } else if (v > maxe[cc]) {
                  harray[cc][v] = 1;
                  for (m2 = (v - 1); m2 > maxe[cc]; m2--) harray[cc][m2] = 0;
                  maxe[cc] = v;
                } else {
                  harray[cc][v]++;
                }
              }
            }

            first_value = false;
          }
        }
      }
    }

    // Work out what rank to use for this mask
    ig1 = (int)floor(g1 * (double)(p - 1) + 0.5);
    if ((ig1 == 0) && (p > 3) && (g1 > 0.0)) ig1 = 1;
    ig2 = ig1;

    // forward operation using these updated values, and the local mask
    for (cc = 0; cc < c; cc++) {
      v1[cc] = iarray[(i * w + j) * c + cc];
      v2[cc] = v1[cc];

      if (oa1 != NULL) {
        mtotal = ig1 + 1;
        for (m = mine[cc]; m <= maxe[cc]; m++) {
          mtotal -= harray[cc][m];
          if (mtotal <= 0) break;
        }
        if (mtotal <= 0) v1[cc] = m + minval;
      }
      if (oa2 != NULL) {
        mtotal = ig2 + 1;
        for (m = maxe[cc]; m >= mine[cc]; m--) {
          mtotal -= harray[cc][m];
          if (mtotal <= 0) break;
        }
        if (mtotal <= 0) v2[cc] = m + minval;
      }
    }

    // Adjust histogram back to original mask
    if (!update_whole_mask) {
      if (p < mp) {
        m = (f - i) * l + (f - j);
        for (ii = iimin; ii <= iimax; ii++) {
          for (jj = jjmin; jj <= jjmax; jj++) {
            mii2 = ii * l + m + jj;
            if ((mask[mii2] > 0) && (mmask[mii2] <= 0)) {
              mii = ii * w + jj;
              for (cc = 0; cc < c; cc++) {
                v = iarray[mii * c + cc] - minval;
                harray[cc][v]++;
              }
            }
          }
        }
      }
    }

    // Perform inverse operation with local mask, or just write value to output array
    if (inverse) {
      m = (f - i) * l + (f - j);
      for (ii = iimin; ii <= iimax; ii++) {
        for (jj = jjmin; jj <= jjmax; jj++) {
          mii = ii * l + m + jj;
          if (mmask[mii] > 0) {
            mii2 = ii * w + jj;
            bstart = bins - pa[mii2];
            pa[mii2]++;
            if (reduce_memory) mii2 = ((ii % l) * w) + jj;
            if (bstart < 1) bstart = 1;
            for (cc = 0; cc < c; cc++) {
              mii = (mii2 * c + cc) * bins;
              if (oa1 != NULL) {
                if (v1[cc] > oa1[mii + bstart - 1]) {
                  for (b = bstart; b < bins; b++) {
                    if (v1[cc] <= oa1[mii + b]) break;
                    oa1[mii + b - 1] = oa1[mii + b];
                  }
                  oa1[mii + b - 1] = v1[cc];
                }
              }
              if (oa2 != NULL) {
                if (v2[cc] < oa2[mii + bstart - 1]) {
                  for (b = bstart; b < bins; b++) {
                    if (v2[cc] >= oa2[mii + b]) break;
                    oa2[mii + b - 1] = oa2[mii + b];
                  }
                  oa2[mii + b - 1] = v2[cc];
                }
              }
            }
          }
        }
      }
      ii = (i * w + j) * c;
      for (cc = 0; cc < c; cc++) {
        if (oa1 != NULL) forward1[ii + cc] = v1[cc];
        if (oa2 != NULL) forward2[ii + cc] = v2[cc];
      }
    } else {
      ii = (i * w + j) * c;
      for (cc = 0; cc < c; cc++) {
        if (oa1 != NULL) oa1[ii + cc] = v1[cc];
        if (oa2 != NULL) oa2[ii + cc] = v2[cc];
      }
    }

    // increment index for next time
    switch (dir) {
    case DOWN:
      if (j == (w - 1)) {
        dir = LEFT;
        j--;
      } else {
        dir = RIGHT;
        j++;
      }
      break;
    case RIGHT:
      if (j == (w - 1)) {
        dir = DOWN;
        i++;
      } else {
        j++;
      }
      break;
    case LEFT:
      if (j == 0) {
        dir = DOWN;
        i++;
      } else {
        j--;
      }
      break;
    }

    if (reduce_memory && (dir == DOWN) && inverse && (i > f)) {
      // create output for (i-f-1) row, which we are about to re-use
      // then re-initialise
      for (ii = 0; ii < w; ii++) {
        jj = (i - f - 1) * w + ii;

        // Work out what centile to use - never a greater rank than bins - 1
        ig1r = (int)floor(g1 * (double)(pa[jj] - 1) + 0.5);
        if ((ig1r == 0) && (pa[jj] > 3) && (g1 > 0.0)) ig1r = 1;
        ig1r = bins - 1 - ig1r;
        ig2r = ig1r;

        for (cc = 0; cc < c; cc++) {
          iii = jj * c + cc;
          mii2 = (((i - f - 1) % l) * w + ii) * c + cc;
          if (oa1 != NULL) {
            for (b = ig1r; b < bins; b++) {
              if (oa1[mii2 * bins + b] >= forward1[iii]) {
                oarray1[iii] = oa1[mii2 * bins + b];
                break;
              }
            }
            for (b = 0; b < bins; b++) oa1[mii2 * bins + b] = minval - 1;
          }
          if (oa2 != NULL) {
            for (b = ig2r; b < bins; b++) {
              if (oa2[mii2 * bins + b] <= forward2[iii]) {
                oarray2[iii] = oa2[mii2 * bins + b];
                break;
              }
            }
            for (b = 0; b < bins; b++) oa2[mii2 * bins + b] = maxval + 1;
          }
        }
      }
    }

  }

  // Possible read sorted value from oa1 and oa2
  if (inverse && (bins > 1)) {
    if (reduce_memory) jj = i - f; else jj = 0;
    for (i = jj; i < h; i++) {
      for (j = 0; j < w; j++) {
        ii = i * w + j;

        // Work out what centile to use - never a greater rank than bins - 1
        ig1r = (int)floor(g1 * (double)(pa[ii] - 1) + 0.5);
        if ((ig1r == 0) && (pa[ii] > 3) && (g1 > 0.0)) ig1r = 1;
        ig1r = bins - 1 - ig1r;
        ig2r = ig1r;

        for (cc = 0; cc < c; cc++) {
          iii = ii * c + cc;
          if (reduce_memory) mii2 = ((i % l) * w + j) * c + cc; else mii2 = iii;
          if (oa1 != NULL) {
            for (b = ig1r; b < bins; b++) {
              if (oa1[mii2 * bins + b] >= forward1[iii]) {
                oarray1[iii] = oa1[mii2 * bins + b];
                break;
              }
            }
          }
          if (oa2 != NULL) {
            for (b = ig2r; b < bins; b++) {
              if (oa2[mii2 * bins + b] <= forward2[iii]) {
                oarray2[iii] = oa2[mii2 * bins + b];
                break;
              }
            }
          }
        }
      }
    }
    if (oa1 != NULL) {
      delete[] oa1;
      delete[] forward1;
    }
    if (oa2 != NULL) {
      delete[] oa2;
      delete[] forward2;
    }
    delete[] pa;
  }

  // free variables
  for (cc = 0; cc < c; cc++) {
    delete[] harray[cc];
  }
  delete[] mmask;
  delete[] mask;

  return;
}

// double precision version using sorting
void xrankopen_2d(double *iarray, int w, int h, int c, int f, double t, double g1, int mlevel, double *oarray1, double *oarray2, double *imarray, double *igarray)
{
  long i, j, cc, p, mp, ii, jj, iii, l, l2;
  int m, ig1, ig1r, ig2, ig2r, mii, mii2;
  int *mmask, *mask, bins, b, bstart, *pa;
  double maxval, minval, v, vrgb, v1[4], v2[4];
  double *igrey, *oa1, *oa2, *forward1, *forward2;
  double minv, maxv, thresh, mthresh, g, ming, maxg, minrgb[4], maxrgb[4];
  double *harray[4];
  double dm, dm2;
  dir_t dir;
  bool inverse = true;
  bool reduce_memory = true;
  bool use_raw_rgb_threshold = (c > 1) && !((mlevel == 0) || (mlevel == -1));
  bool use_local_threshold; // set once we know max and min values

  // return if incorrect inputs
  if ((f < 1) || (w < 1) || (h < 1) || (c < 1) || (c > 4) || (g1 < 0.0) || (g1 > 0.5)) return;
  if ((oarray1 == NULL) && (oarray2 == NULL)) return;
  if ((imarray == NULL) || ((c > 1) && (igarray == NULL))) return;

  // Just copy image across if it is too small
  if ((f > (h / 2)) || (f > (w / 2))) {
    if (oarray1 != NULL) {
      for (i = 0; i < w * h * c; i++) {
        oarray1[i] = iarray[i];
      }
    }
    if (oarray2 != NULL) {
      for (i = 0; i < w * h * c; i++) {
        oarray2[i] = iarray[i];
      }
    }
    return;
  }


  // Calculate mask and bin sizes, and centile location
  l = 2 * f + 1;
  l2 = l * l;
  mp = l2;
  ig1 = (int)floor(g1 * (double)(mp - 1) + 0.5);
  if ((ig1 == 0) && (mp > 3) && (g1 > 0.0)) ig1 = 1;
  bins = ig1 + 1;
  if ((g1 == 0.5) && (oarray2 == NULL)) {
    // special case for median
    inverse = false;
  } else {
    // can do the whole process in one pass
    inverse = true;
  }

  // set up starting and ending points based on rectangular mask
  mmask = new int[l2];
  mask = new int[l2];
  for (i = 0; i < l2; i++) {
    mask[i] = 1;
    mmask[i] = 1;
  }

  // Check for colour input, and adjust input image and threshold
  if (c > 1) {
    igrey = igarray;
  } else {
    igrey = iarray;
  }

  // Compensation due to prior filtering, which reduces the noise
  if (mlevel == 0) {
    dm = sqrt(M_PI) * ((double)f + 0.5) * 1.2 * 0.85;
    // Further compensation to reduce the proportion of pixels being used as the mask gets larger
    dm2 = dm * erf(((double)f + 0.5) * 1.2 / 24.0);
  } else {
    // After the first level, the noise is correlated so filtering reduces it less
    // The filtering is twice as long as the previous level, so the effective filter length is about 2.0
    dm = 2.0;
    dm2 = dm * erf(((double)f + 0.5) * 1.2 / 12.0);
  }
  thresh = t;
  if (c > 1) thresh /= sqrt(3.0);
  mthresh = thresh / dm2;

  // check for input array range, use one for all components
  minval = iarray[0];
  maxval = minval;
  for (i = 1; i < w * h * c; i++) {
    if (iarray[i] > maxval) maxval = iarray[i];
    if (iarray[i] < minval) minval = iarray[i];
  }

  // use range to work out if we are in a low noise scenario
  use_local_threshold = (t > (maxval - minval) / 20.0) && (mlevel == 0);

  // initialise output if doing inverse in same loop as forward
  oa1 = oarray1;
  oa2 = oarray2;
  if (inverse) {
    if (oa1 != NULL) forward1 = new double[w * h * c];
    if (oa2 != NULL) forward2 = new double[w * h * c];
    pa = new int[w * h];
    for (i = 0; i < w * h; i++) pa[i] = 0;
    if (bins == 1) reduce_memory = false;
    if (reduce_memory) {
      // this only stores w * l rather than w * h
      if (oarray1 != NULL) oa1 = new double[w * l * c * bins];
      if (oarray2 != NULL) oa2 = new double[w * l * c * bins];
      if (oa1 != NULL) {
        for (i = 0; i < w * l * c * bins; i++) oa1[i] = minval - 1;
      }
      if (oa2 != NULL) {
        for (i = 0; i < w * l * c * bins; i++) oa2[i] = maxval + 1;
      }
    } else {
      if (bins > 1) {
        if (oarray1 != NULL) oa1 = new double[w * h * c * bins];
        if (oarray2 != NULL) oa2 = new double[w * h * c * bins];
      }
      if (oa1 != NULL) {
        for (i = 0; i < w * h * c * bins; i++) oa1[i] = minval - 1;
      }
      if (oa2 != NULL) {
        for (i = 0; i < w * h * c * bins; i++) oa2[i] = maxval + 1;
      }
    }
  }

  // create space for sorted list
  // in this case we are updating this each time
  mp = 0;
  for (cc = 0; cc < c; cc++) {
    harray[cc] = new double[l2];
  }

  // start walking through image using forwards then backwards approach to make the addition
  // of extra values to the sorted array more efficient
  i = 0;
  j = 0;
  dir = DOWN;
  while ((i < h) && (j < w)) {

    // Construct threshold mask for this location
    int iimax, iimin, jjmax, jjmin;
    iimin = i - f;
    if (iimin < 0) iimin = 0;
    iimax = i + f;
    if (iimax >= h) iimax = h - 1;
    jjmin = j - f;
    if (jjmin < 0) jjmin = 0;
    jjmax = j + f;
    if (jjmax >= w) jjmax = w - 1;
    p = mp;

    // Check distance of current value from min and max in histogram
    double tfact = 1.0;

    // This only happens in the first full size frame if doing multiresolution
    if (use_local_threshold) {
      double mrange, mminv;
      v = igrey[i * w + j];
      m = (f - i) * l + (f - j);
      maxg = v;
      ming = v;
      for (ii = iimin; ii <= iimax; ii++) {
        for (jj = jjmin; jj <= jjmax; jj++) {
          if (mask[ii * l + m + jj] > 0) { // in overall mask shape
            mii = ii * w + jj;
            if (igrey[mii] > maxg) maxg = igrey[mii];
            if (igrey[mii] < ming) ming = igrey[mii];
          }
        }
      }
      mrange = (maxg - ming);
      mminv = fabs(v - ming);
      if (fabs(maxg - v) < mminv) mminv = fabs(maxg - v);
      mrange = mrange / thresh;
      mminv = mminv / thresh;
      // These limits presume that 'thresh' has been set to four times the std of the noise.
      double low_filter, high_filter, range_limit, low_range_limit;
      low_filter = 0.7;
      high_filter = 1.6;
      range_limit = 1.57;
      low_range_limit = 0.9;
      if (mrange > range_limit) {
        // Reduce the threshold limits if the data range is very high
        // (and therefore the noise range is a small proportion of this)
        tfact = low_filter;
      } else if (mrange < low_range_limit) {
        // Increase the threshold limits if we are pretty sure that the noise occupies the whole range
        tfact = high_filter;
      } else {
        // Otherwise a blend between the two scenarios above
        tfact = high_filter - (mrange - low_range_limit) * (high_filter - low_filter) / (range_limit - low_range_limit);
      }
      if ((tfact < 1.0) && (mminv < low_filter)) {
        // As long as we are in the relatively high range case, also increase threshold if we are really close to the range edges
        // (since then we know the distribution is skewed towards other direction)
        tfact = 1.0 - mminv * (1.0 - tfact) / low_filter;
      }
    }

    // compensate for higher chance to de-select mask pixel because considering several channels
    // also, this factor isn't variable like tfact
    double tfactrgb = sqrt((double)c) * sqrt(3.0);
    if (tfact > 1.0) tfactrgb *= tfact;

    mii = i * w + j;

    // Work out what range of pixel values to include from (probably filtered) grey mask image imarray
    v = imarray[mii];
    minv = v - (mthresh * tfact);
    maxv = v + (mthresh * tfact);

    // Work out what range of pixel values to include from (un-filtered) grey image
    g = igrey[mii];
    ming = g - (thresh * tfact);
    maxg = g + (thresh * tfact);

    // Use raw threshold / average comparison to adjust mthresh
    // which allows us to ignore filtered image if it isn't from the same distribution as un-filtered data
    // the 0.5 allows for half the threshold range in the un-filtered data,
    // the 0.5/dm allows for additional half threshold range in the filtered data
    double test = thresh * (0.5 + 0.5 / dm);
    double vg = fabs(v - g);
    if (vg > test) {
      minv -= 2 * (vg - test);
      maxv += 2 * (vg - test);
    }

    if (use_raw_rgb_threshold) {
      for (cc = 0; cc < c; cc++) {
        vrgb = iarray[mii * c + cc];
        minrgb[cc] = vrgb - (thresh * tfactrgb);
        maxrgb[cc] = vrgb + (thresh * tfactrgb);
      }
    }

    // Threshold based on imarray (and possibly also imarray2 and igrey)
    for (m = 0; m < l2; m++) mmask[m] = 0;
    mmask[(l2 - 1) / 2] = 1; // always include central pixel
    p = 1;
    m = (f - i) * l + (f - j);
    for (ii = iimin; ii <= iimax; ii++) {
      for (jj = jjmin; jj <= jjmax; jj++) {
        mii2 = ii * l + m + jj;
        if (2 * mii2 == (l2 - 1)) continue;
        if (mask[mii2] > 0) { // in overall mask shape
          mii = ii * w + jj;
          if ((imarray[mii] >= minv) && (imarray[mii] <= maxv)) { // within threshold of imarray
            if ((igrey[mii] >= ming) && (igrey[mii] <= maxg)) { // within threshold of igrey
              mmask[mii2] = 1;
              p++;
            }
          }
          if (use_raw_rgb_threshold && (mmask[mii2] == 1)) {
            for (cc = 0; cc < c; cc++) {
              vrgb = iarray[mii * c + cc];
              if ((vrgb < minrgb[cc]) || (vrgb > maxrgb[cc])) {
                mmask[mii2] = 0;
                p--;
                break;
              }
            }
          }
        }
      }
    }

    // Update current sorted list from mmask
    int m2 = 0;
    m = (f - i) * l + (f - j);
    for (ii = iimin; ii <= iimax; ii++) {
      for (jj = jjmin; jj <= jjmax; jj++) {
        mii2 = ii * l + m + jj;
        if (mmask[mii2] > 0) {
          mii = ii * w + jj;
          for (cc = 0; cc < c; cc++) {
            harray[cc][m2] = iarray[mii * c + cc];
          }
          m2++;
        }
      }
    }
    for (cc = 0; cc < c; cc++) {
      sort(harray[cc], harray[cc] + m2);
    }

    // Work out what rank to use for this mask
    ig1 = (int)floor(g1 * (double)(p - 1) + 0.5);
    if ((ig1 == 0) && (p > 3) && (g1 > 0.0)) ig1 = 1;
    ig2 = p - 1 - ig1;

    // forward operation using these updated values, and the local mask
    for (cc = 0; cc < c; cc++) {
      v1[cc] = iarray[(i * w + j) * c + cc];
      v2[cc] = v1[cc];
      if (oa1 != NULL) {
        v1[cc] = harray[cc][ig1];
      }
      if (oa2 != NULL) {
        v2[cc] = harray[cc][ig2];
      }
    }

    // Perform inverse operation with local mask, or just write value to output array
    if (inverse) {
      m = (f - i) * l + (f - j);
      for (ii = iimin; ii <= iimax; ii++) {
        for (jj = jjmin; jj <= jjmax; jj++) {
          mii = ii * l + m + jj;
          if (mmask[mii] > 0) {
            mii2 = ii * w + jj;
            bstart = bins - pa[mii2];
            pa[mii2]++;
            if (reduce_memory) mii2 = ((ii % l) * w) + jj;
            if (bstart < 1) bstart = 1;
            for (cc = 0; cc < c; cc++) {
              mii = (mii2 * c + cc) * bins;
              if (oa1 != NULL) {
                if (v1[cc] > oa1[mii + bstart - 1]) {
                  for (b = bstart; b < bins; b++) {
                    if (v1[cc] <= oa1[mii + b]) break;
                    oa1[mii + b - 1] = oa1[mii + b];
                  }
                  oa1[mii + b - 1] = v1[cc];
                }
              }
              if (oa2 != NULL) {
                if (v2[cc] < oa2[mii + bstart - 1]) {
                  for (b = bstart; b < bins; b++) {
                    if (v2[cc] >= oa2[mii + b]) break;
                    oa2[mii + b - 1] = oa2[mii + b];
                  }
                  oa2[mii + b - 1] = v2[cc];
                }
              }
            }
          }
        }
      }
      ii = (i * w + j) * c;
      for (cc = 0; cc < c; cc++) {
        if (oa1 != NULL) forward1[ii + cc] = v1[cc];
        if (oa2 != NULL) forward2[ii + cc] = v2[cc];
      }
    } else {
      ii = (i * w + j) * c;
      for (cc = 0; cc < c; cc++) {
        if (oa1 != NULL) oa1[ii + cc] = v1[cc];
        if (oa2 != NULL) oa2[ii + cc] = v2[cc];
      }
    }

    // increment index for next time
    switch (dir) {
    case DOWN:
      if (j == (w - 1)) {
        dir = LEFT;
        j--;
      } else {
        dir = RIGHT;
        j++;
      }
      break;
    case RIGHT:
      if (j == (w - 1)) {
        dir = DOWN;
        i++;
      } else {
        j++;
      }
      break;
    case LEFT:
      if (j == 0) {
        dir = DOWN;
        i++;
      } else {
        j--;
      }
      break;
    }

    if (reduce_memory && (dir == DOWN) && inverse && (i > f)) {
      // create output for (i-f-1) row, which we are about to re-use
      // then re-initialise
      for (ii = 0; ii < w; ii++) {
        jj = (i - f - 1) * w + ii;

        // Work out what centile to use - never a greater rank than bins - 1
        ig1r = (int)floor(g1 * (double)(pa[jj] - 1) + 0.5);
        if ((ig1r == 0) && (pa[jj] > 3) && (g1 > 0.0)) ig1r = 1;
        ig1r = bins - 1 - ig1r;
        ig2r = ig1r;

        for (cc = 0; cc < c; cc++) {
          iii = jj * c + cc;
          mii2 = (((i - f - 1) % l) * w + ii) * c + cc;
          if (oa1 != NULL) {
            for (b = ig1r; b < bins; b++) {
              if (oa1[mii2 * bins + b] >= forward1[iii]) {
                oarray1[iii] = oa1[mii2 * bins + b];
                break;
              }
            }
            for (b = 0; b < bins; b++) oa1[mii2 * bins + b] = minval - 1;
          }
          if (oa2 != NULL) {
            for (b = ig2r; b < bins; b++) {
              if (oa2[mii2 * bins + b] <= forward2[iii]) {
                oarray2[iii] = oa2[mii2 * bins + b];
                break;
              }
            }
            for (b = 0; b < bins; b++) oa2[mii2 * bins + b] = maxval + 1;
          }
        }
      }
    }

  }

  // Possible read sorted value from oa1 and oa2
  if (inverse && (bins > 1)) {
    if (reduce_memory) jj = i - f; else jj = 0;
    for (i = jj; i < h; i++) {
      for (j = 0; j < w; j++) {
        ii = i * w + j;

        // Work out what centile to use - never a greater rank than bins - 1
        ig1r = (int)floor(g1 * (double)(pa[ii] - 1) + 0.5);
        if ((ig1r == 0) && (pa[ii] > 3) && (g1 > 0.0)) ig1r = 1;
        ig1r = bins - 1 - ig1r;
        ig2r = ig1r;

        for (cc = 0; cc < c; cc++) {
          iii = ii * c + cc;
          if (reduce_memory) mii2 = ((i % l) * w + j) * c + cc; else mii2 = iii;
          if (oa1 != NULL) {
            for (b = ig1r; b < bins; b++) {
              if (oa1[mii2 * bins + b] >= forward1[iii]) {
                oarray1[iii] = oa1[mii2 * bins + b];
                break;
              }
            }
          }
          if (oa2 != NULL) {
            for (b = ig2r; b < bins; b++) {
              if (oa2[mii2 * bins + b] <= forward2[iii]) {
                oarray2[iii] = oa2[mii2 * bins + b];
                break;
              }
            }
          }
        }
      }
    }
    if (oa1 != NULL) {
      delete[] oa1;
      delete[] forward1;
    }
    if (oa2 != NULL) {
      delete[] oa2;
      delete[] forward2;
    }
    delete[] pa;
  }

  // free variables
  for (cc = 0; cc < c; cc++) {
    delete[] harray[cc];
  }
  delete[] mmask;
  delete[] mask;

  return;
}


// Call this function using [X X2] = xrankopen2(A, t, f, Msm, M, centile, mlevel)
#define A_IN prhs[0]
#define t_IN prhs[1]
#define f_IN prhs[2]
#define Msm_IN prhs[3]
#define M_IN prhs[4]
#define centile_IN prhs[5]
#define mlevel_IN prhs[6]
#define X_OUT plhs[0]
#define X2_OUT plhs[1]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *Ad, *Adtmp, *Msmd, *Md, *Xd, *Xdtmp, *X2d, *X2dtmp, td, g;
  long w, h, c, cc, f, t, l, ndims;
  int *A, *Atmp, *X, *Xtmp, *X2, *X2tmp, *Msm, *M;
  unsigned char *uc;
  long i;
  bool is_double, is_dual;

  // check we have consistent input and output arguments
  if (nrhs < 5) {
    mexErrMsgTxt("A, t, f, Msm and M must all be defined.");
  } else if (nrhs > 7) {
    mexErrMsgTxt("No more than 7 input arguments.");
  }
  if (nlhs < 1) {
    mexErrMsgTxt("Must be at least one output argument.");
  } else if (nlhs > 2) {
    mexErrMsgTxt("Must be no more than two output arguments.");
  }
  is_dual = (nlhs == 2);

  // check input data is of the right type
  if (mxIsComplex(A_IN) || mxIsSparse(A_IN) || ((mxGetClassID(A_IN) != mxDOUBLE_CLASS) && (mxGetClassID(A_IN) != mxUINT8_CLASS) && (mxGetClassID(A_IN) != mxINT32_CLASS))) {
    mexErrMsgTxt("Input data must be a matrix of double, uint8 or int32 real values.");
  }
  c = 0;
  if (mxGetNumberOfDimensions(A_IN) == 2) {
    c = 1;
    w = (int)mxGetM(A_IN);
    h = (int)mxGetN(A_IN);
  } else if (mxGetNumberOfDimensions(A_IN) == 3) {
    w = (int)(mxGetDimensions(A_IN)[0]);
    h = (int)(mxGetDimensions(A_IN)[1]);
    c = (int)(mxGetDimensions(A_IN)[2]);
  }
  if ((c < 1) || (c > 4)) {
    mexErrMsgTxt("Input data must be a 2D matrix or a 3D matrix with colour as the third dimension.");
  }
  if (mxGetClassID(A_IN) == mxDOUBLE_CLASS) {
    is_double = true;
  } else {
    is_double = false;
  }

  // Check smoothed and raw mask images 
  if (mxIsComplex(Msm_IN) || mxIsSparse(Msm_IN) || (mxGetClassID(Msm_IN) != mxGetClassID(A_IN)) || (mxGetNumberOfDimensions(Msm_IN) != 2)) {
    mexErrMsgTxt("Msm must be an array of real values of the same type as the input data A.");
  }
  if (((int)mxGetM(Msm_IN) != w) || ((int)mxGetN(Msm_IN) != h)) {
    mexErrMsgTxt("Msm must have the same first two dimension as the input data A.");
  }
  if (mxIsComplex(M_IN) || mxIsSparse(M_IN) || (mxGetClassID(M_IN) != mxGetClassID(A_IN)) || (mxGetNumberOfDimensions(M_IN) != 2)) {
    mexErrMsgTxt("M must be an array of real values of the same type as the input data A.");
  }
  if (((int)mxGetM(M_IN) != w) || ((int)mxGetN(M_IN) != h)) {
    mexErrMsgTxt("M must have the same first two dimension as the input data A.");
  }

  // Read in threshold and filter size, these must exist
  // Check for threshold and read in immediately
  if (mxIsComplex(t_IN) || (mxGetNumberOfElements(t_IN) != 1)) {
    mexErrMsgTxt("t must be a scalar.");
  } else {
    if (is_double) {
      td = (double)mxGetScalar(t_IN);
    } else {
      td = (double)mxGetScalar(t_IN);
    }
    if (td <= 0) {
      td = 1e10;
      t = 32768;
    } else {
      t = (int)round(td);
    }
  }
  if (mxIsComplex(f_IN) || (mxGetNumberOfElements(f_IN) != 1)) {
    mexErrMsgTxt("f must be a scalar.");
  } else {
    f = (int)mxGetScalar(f_IN);
    if (f < 1) {
      mexErrMsgTxt("f must be at least one.");
    }
  }

  // Check for optional arguments, and set to default otherwise
  if (nrhs > 5) {
    if (mxIsComplex(centile_IN) || (mxGetNumberOfElements(centile_IN) != 1)) {
      mexErrMsgTxt("centile must be a scalar.");
    } else {
      g = (double)mxGetScalar(centile_IN)/100.0;
    }
    if ((g < 0.0) || (g > 0.5)) {
      mexErrMsgTxt("centile must be between 0 and 50.");
    }
  } else {
    g = 0.08;
  }
  if (nrhs > 6) {
    if (mxIsComplex(mlevel_IN) || (mxGetNumberOfElements(mlevel_IN) != 1)) {
      mexErrMsgTxt("mlevel must be a scalar.");
    } else {
      l = (int)mxGetScalar(mlevel_IN);
    }
  } else {
    l = 0;
  }

  // Get input and output data and call appropriate function
  mwSize dims[3];
  dims[0] = w;
  dims[1] = h;
  dims[2] = c;
  if (c > 1) {
    ndims = 3;
  } else {
    ndims = 2;
  }
  if (is_double) {
    Ad = (double *)mxGetData(A_IN);
    if (c > 1) {
      Adtmp = new double[w * h * c];
      // Need to swap array ordering so colours change fastest
      for (i = 0; i < w * h; i++) {
        for (cc = 0; cc < c; cc++) {
          Adtmp[i * c + cc]= Ad[cc * w * h + i];
        }
      }
    } else {
      Adtmp = Ad;
    }
    Md = (double *)mxGetData(M_IN);
    Msmd = (double *)mxGetData(Msm_IN);
    X_OUT = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
    Xd = (double *)mxGetData(X_OUT);
    if (c > 1) {
      Xdtmp = new double[w * h * c];
    } else {
      Xdtmp = Xd;
    }
    if (is_dual) {
      X2_OUT = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
      X2d = (double *)mxGetData(X2_OUT);
      if (c > 1) {
        X2dtmp = new double[w * h * c];
      } else {
        X2dtmp = X2d;
      }
      xrankopen_2d(Adtmp, w, h, c, f, td, g, l, Xdtmp, X2dtmp, Msmd, Md);
      if (c > 1) {
        // Need to swap array ordering
        for (i = 0; i < w * h; i++) {
          for (cc = 0; cc < c; cc++) {
            X2d[cc * w * h + i] = X2dtmp[i * c + cc];
          }
        }
        delete[] X2dtmp;
      }
    } else {
      xrankopen_2d(Adtmp, w, h, c, f, td, g, l, Xdtmp, NULL, Msmd, Md);
    }
    if (c > 1) {
      // Need to swap array ordering
      for (i = 0; i < w * h; i++) {
        for (cc = 0; cc < c; cc++) {
          Xd[cc * w * h + i] = Xdtmp[i * c + cc];
        }
      }
      delete[] Xdtmp;
      delete[] Adtmp;
    }
  } else if (mxGetClassID(A_IN) == mxUINT8_CLASS) {
    uc = (unsigned char *)mxGetData(A_IN);
    A = new int[w * h * c];
    // Need to swap array ordering so colours change fastest
    for (i = 0; i < w * h; i++) {
      for (cc = 0; cc < c; cc++) {
        A[i*c + cc] = (int)uc[cc * w * h + i];
      }
    }
    uc = (unsigned char *)mxGetData(M_IN);
    M = new int[w * h];
    for (i = 0; i < w * h; i++) M[i] = (int)uc[i];
    uc = (unsigned char *)mxGetData(Msm_IN);
    Msm = new int[w * h];
    for (i = 0; i < w * h; i++) Msm[i] = (int)uc[i];
    X = new int[w * h * c];
    if (is_dual) {
      X2 = new int[w * h * c];
      xrankopen_2d(A, w, h, c, f, t, g, l, X, X2, Msm, M);
      X2_OUT = mxCreateNumericArray(ndims, dims, mxUINT8_CLASS, mxREAL);
      uc = (unsigned char *)mxGetData(X2_OUT);
      // Need to swap array ordering
      for (i = 0; i < w * h; i++) {
        for (cc = 0; cc < c; cc++) {
          uc[cc * w * h + i] = (unsigned char)X2[i * c + cc];
        }
      }
      delete[] X2;
    } else {
      xrankopen_2d(A, w, h, c, f, t, g, l, X, NULL, Msm, M);
    }
    X_OUT = mxCreateNumericArray(ndims, dims, mxUINT8_CLASS, mxREAL);
    uc = (unsigned char *)mxGetData(X_OUT);
    // Need to swap array ordering
    for (i = 0; i < w * h; i++) {
      for (cc = 0; cc < c; cc++) {
        uc[cc * w * h + i] = (unsigned char)X[i * c + cc];
      }
    }
    delete[] X;
    delete[] A;
    delete[] Msm;
    delete[] M;
  } else {
    A = (int *)mxGetData(A_IN);
    if (c > 1) {
      Atmp = new int[w * h * c];
      // Need to swap array ordering so colours change fastest
      for (i = 0; i < w * h; i++) {
        for (cc = 0; cc < c; cc++) {
          Atmp[i * c + cc]= A[cc * w * h + i];
        }
      }
    } else {
      Atmp = A;
    }
    M = (int *)mxGetData(M_IN);
    Msm = (int *)mxGetData(Msm_IN);
    X_OUT = mxCreateNumericArray(ndims, dims, mxINT32_CLASS, mxREAL);
    X = (int *)mxGetData(X_OUT);
    if (c > 1) {
      Xtmp = new int[w * h * c];
    } else {
      Xtmp = X;
    }
    if (is_dual) {
      X2_OUT = mxCreateNumericArray(ndims, dims, mxINT32_CLASS, mxREAL);
      X2 = (int *)mxGetData(X2_OUT);
      if (c > 1) {
        X2tmp = new int[w * h * c];
      } else {
        X2tmp = X2;
      }
      xrankopen_2d(Atmp, w, h, c, f, t, g, l, Xtmp, X2tmp, Msm, M);
      if (c > 1) {
        // Need to swap array ordering
        for (i = 0; i < w * h; i++) {
          for (cc = 0; cc < c; cc++) {
            X2[cc * w * h + i] = X2tmp[i * c + cc];
          }
        }
        delete[] X2tmp;
      }
    } else {
      xrankopen_2d(Atmp, w, h, c, f, t, g, l, Xtmp, NULL, Msm, M);
    }
    if (c > 1) {
      // Need to swap array ordering
      for (i = 0; i < w * h; i++) {
        for (cc = 0; cc < c; cc++) {
          X[cc * w * h + i] = Xtmp[i * c + cc];
        }
      }
      delete[] Xtmp;
      delete[] Atmp;
    }
  }

}
