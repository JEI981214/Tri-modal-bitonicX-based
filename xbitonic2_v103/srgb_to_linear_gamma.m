function X = srgb_to_linear_gamma(A, gamma)
%SRGB_TO_LINEAR_GAMMA remove gamma transform from SRGB image.
%   X = SRGB_TO_LINEAR_GAMMA(A) removes the standard SRGB gamma transform
%   from the SRGB image in A.
%
%   X = SRGB_TO_LINEAR_GAMMA(A, gamma) removes the gamma transform given
%   a specific value of gamma (default is 2.4).   
%  
%   A must be double precision, with an assumed range of 0.0 to 1.0
%
%   See also XBITONIC2, MXBITONIC2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2021.

% check for appropriate inputs
narginchk(1, 2);
if (nargin < 2)
  gamma = 2.4;
end

% also check for appropriate outputs
nargoutchk(1, 1);

% set up correct transform values
b = 0.055;
if (gamma == 2.4)
  % use standard values for sRGB
  c = 12.92;
  d = 0.0031308;
else
  % work out other parameters from gamma value
  c = ((1.0 + b)^gamma) * ((gamma - 1.0)^(gamma - 1.0)) / ((b^(gamma - 1.0)) * (gamma^gamma));
  d = (b / (gamma - 1.0)) / c;
end
X = zeros(size(A));
a = 1.0/(1.0 + b);
d = d * c;
c = 1.0 / c;

% low values are linear
low = (A < d);
X(low) = A(low) * c;

% high values are standard gamma
high = ~low;
X(high) = (a * (A(high) + b)).^gamma;

