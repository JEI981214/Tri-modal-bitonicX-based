function [X, X2] = xrankopen2(A, t, f, Msm, M, centile, mlevel)
%XRANKOPEN2 flexible threshold-based opening and closing.
%   X = XRANKOPEN2(A, t, f, Msm, M) performs a flexible threshold-based
%   opening on A, using masks formed from the smoothed grey-scale image Msm
%   and non-smoothed version M. These are thresholded based on the level t,
%   which should be set based on the expected noise in the image A.
%
%   A can be a 2D grey-scale array or a 3D colour array, with up to four
%   colour channels.
%
%   t is a threshold which should be set to the approximate maximum range
%   of noise in the image, usually four times the standard deviation of the
%   noise.
%
%   f is the filter half length for the mask size l x l, where l = 2*f+1.
%   It should also correspond to the approximate filter size in creating
%   Msm from M.
%
%   Msm and M must be of the same data type as A, but can only be a 2D
%   array. If A is a 2D grey-scale array, M should be the same as A.
%
%   The algorithm uses a fast histogram-based technique if A is uint8 or
%   int32, or a slightly less fast sorting-based technique for doubles.
%
%   X = XRANKOPEN2(A, t, f, Msm, M, centile) also uses centile to determine
%   the centile to use in opening and closing: the default is 8. Set this
%   to 50 to output the median in X.
%
%   X = XRANKOPEN2(A, t, f, Msm, M, centile, mlevel) uses mlevel to
%   indicate the level when using multi-resolution filtering, with zero the
%   initial level, positive integers for increasingly smaller levels, and
%   negative integers for the reverse process, finishing at -1 for the last
%   level.
%
%   [X, X2] = XRANKOPEN2(A, t, f, Msm, M) outputs the closing operation in
%   X2 as well as the opening in X.
%
%   See also XBITONIC2, MXBITONIC2, RANKOPEN2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2021.

