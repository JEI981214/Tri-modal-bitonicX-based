function X = linear_gamma_to_srgb(A, gamma)
%LINEAR_GAMMA_TO_SRGB apply gamma transform to colour image.
%   X = LINEAR_GAMMA_TO_SRGB(A) applies the standard SRGB gamma transform
%   to the colour image in A.
%  
%   X = LINEAR_GAMMA_TO_SRGB(A, gamma) applies the gamma transform given
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
gamma = 1.0 / gamma;
a = (1.0 + b);

% apply transform
X = zeros(size(A));

% low values are linear
low = (A < d);
X(low) = c * A(low);

% high values are standard gamma
high = ~low;
X(high) = a * A(high).^gamma - b;

