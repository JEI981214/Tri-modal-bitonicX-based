function [Xn, X] = estimate_image_noise(A, skip)
%ESTIMATE_IMAGE_NOISE estimate noise STD from image.
%   Xn = ESTIMATE_IMAGE_NOISE(A) estimates the STD of the image noise in
%   A and returns the value in Xn.
%
%   Xn = ESTIMATE_IMAGE_NOISE(A, skip) processes every 'skip' pixels in the
%   image, in both directions, and repeats for all possibly offsets up to
%   skip. This removes the effect of correlated pixels, particularly for
%   SRGB data which has been created from a Bayer image.
%
%   The method is based on taking the absolute value of the discrete
%   laplacian to remove signal, then looking for a low centile of this to
%   further discourage residual signal effects.
%
%   [Xn, X] = ESTIMATE_IMAGE_NOISE(A) also returns a noise image in X.
%
%   See also XBITONIC2, MXBITONIC2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2021.

% check inputs and outputs
narginchk(1,2);
if (nargin < 2)
  skip = 1;
else
  skip = round(skip);
  if (skip < 1)
    skip = 1;
  end
end

nargoutchk(1,2);

% sort out scale for noise and ensure grey scale input
centile = 6;
cl = class(A);
if ((ndims(A)==3) && (size(A,3)>1))
  scale = 2.34 * sqrt(size(A,3));
  A = mean(double(A),3);
else
  scale = 2.34;
  A = double(A);
end
if (skip == 1)
  f = 14;
else
  f = round(14 / sqrt(skip));
end

% This process repeats for all skip pixels
C = zeros(size(A));
for x=1:skip
  for y=1:skip
    
    % create sub-image, skipping over correlated pixels
    B = A(x:skip:end, y:skip:end);
    
    % Take second derivatives in both directions to remove most of signal
    m = [1 1:size(B,1) size(B,1)];
    n = [1 1:size(B,2) size(B,2)];
    B = conv2([-1 2 -1], [-1 2 -1], B(m,n), 'valid')/6.0;
    
    % Take low rank of filtered abs of this to ignore any residual signal
    B = anisotropic2(abs(B), 2.6, 0.7);
    j = [0:(2*f)]-f;
    n = 2*f+1;
    M = double(((j'*ones(1,n)).^2 + (ones(n,1)*j).^2)<(f+0.1)^2);
    p = sum(M(:));
    c1 = round((centile/100)*(p-1)) + 1;
    if strcmp(cl, 'uint8') || strcmp(cl, 'int32')
      B = int32(B*4.0); % makes following operation much more efficient
    end
    B = rankopen2(B, M, c1);
    
    % produce full image output
    C(x:skip:end, y:skip:end) = B;
  end
end

% Possibly take max over skip rectangle
if skip > 1
  M = ones(2*round(skip/2) + 1, 2*round(skip/2) + 1);
  p = length(M(:));
  C = rankopen2(C, M, p, 0, false);
end

% output as average
if strcmp(cl, 'uint8') || strcmp(cl, 'int32')
  C = double(C)/4.0;
end
Xn = mean(C(:)) * scale;
    
% possible output as image
if (nargout > 1)
  X = C * scale;
end
