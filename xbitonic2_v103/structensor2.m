function [X, Y, Z] = structensor2(A, f, type)
%STRUCTENSOR2 returns orientation and anisotropy from the structure tensor.
%   [X, Y] = STRUCTENSOR2(A, f) calculates the local anisotropy in X, and
%   the local feature direction in Y for the image A, using the feature
%   scale given by f.
%
%   If A is a colour image, the colour tensor is calculated and combined
%   with the grey-scale version.
%
%   [X, Y] = STRUCTENSOR2(A, f, 'classic') uses the classic definition of 
%   anisotropy as (1 - l2/l1), where l1 and l2 are the eigenvalues of the
%   structure tensor. Otherwise the default behaviour is to smooth l2
%   first, which improves the response at corners.
%
%   [X, Y, Z] = STRUCTENSOR2(A, f) also returns a normalised corner
%   detection in Z, where both eigenvalues are large and hence anisotropy
%   will still be low. 
%
%   See also GRAD2, GAUSS2, ANISOTROPIC2, XRANKOPEN2, XBITONIC2, MXBITONIC2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2021.


% check for appropriate inputs
narginchk(2, 3);

% also check for appropriate outputs
nargoutchk(2, 3);

% other input checks
if ~(ndims(A) == 2) && ~(ndims(A) == 3)
  error('A must be a grey or colour image.');
end
if nargin<3
  type = 'adapted';
else
  if ~strcmp(type, 'classic')
    type = 'adapted';
  end
end
if nargout == 3
  do_corner = true;
else
  do_corner = false;
end

% calculate gradients in X and Y dimensions, and subsequent components of
% the structure tensor, from grey image
if (ndims(A)==3)
  [X, Y] = grad2(mean(double(A),3));
else
  [X, Y] = grad2(double(A));
end
Txx = X.^2;
Tyy = Y.^2;
Txy = 2*X.*Y;

% possibly also add in colour tensor
if ((ndims(A)==3) && (size(A,3)==3))
  
  % start with 1/4 of grey-scale version
  cTxx = Txx / 4.0;
  cTyy = Tyy / 4.0;
  cTxy = Txy / 4.0;
  
  % take gradient in each colour and add in a further 1/4
  for c = 1:3
    [X, Y] = grad2(double(A(:,:,c)));
    cTxx = cTxx + (X.^2) / 4.0;
    cTyy = cTyy + (Y.^2) / 4.0;
    cTxy = cTxy + (X.*Y) / 2.0;
  end
  
  % replace with colour tensor if it is large enough
  ind = (cTxx + cTyy) > 3.0*(Txx + Tyy);
  Txx(ind) = cTxx(ind) / 3.0;
  Tyy(ind) = cTyy(ind) / 3.0;
  Txy(ind) = cTxy(ind) / 3.0;
  
end

% can now filter this tensor
Txx = gauss2(Txx, f);
Tyy = gauss2(Tyy, f);
Txy = gauss2(Txy, f);

% calculate anisotropy and orientation
T1 = Txx + Tyy;
T2 = Txx - Tyy;
Y = (atan2(Txy, T2) + pi)/2;
T2 = (T2.^2 + Txy.^2).^0.5;
X = zeros(size(Y));
n = T1>0;
if strcmp(type, 'classic')
  X(n) = 1 - (T1(n) - T2(n))./(T1(n) + T2(n));
else
  X(n) = 1 - gauss2((T1(n) - T2(n)),2*f)./(T1(n) + T2(n));
  X(X<0) = 0;
end

% possibly calculate corners
if do_corner
  
  % product of eigenvalues: looking for both large, i.e. a corner
  Z = (T1 - T2) .* (T1 + T2) * 0.5;
  
  % use histogram to get below 90% average
  Zf = Z(:);
  Zfs = Zf*255.0/max(Zf);
  n = cumsum(hist(Zfs, 256))/length(Zfs);
  m = find(n>0.9, 1);
  m = mean(Zf(Zfs<m));
  
  % scale by this 
  minco = 2.0 * m;
  maxco = 40.0 * m;
  Z = (Z - minco)./(maxco - minco);
  Z(Z<0.0) = 0.0;
  
end

