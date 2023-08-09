function [X] = prolongate2(A, sz)
%PROLONGATE2 prolongation operation for use in multi-level methods.
%   X = PROLONGATE2(A) returns a double-size array X which is a
%   prolongation of A using an edge-preserving Catmull-Rom operator.
%
%   X = PROLONGATE2(A, sz) uses sz to ensure that the result of
%   prolongation is a matrix of the same size sz as that which was
%   restricted.
%
%   See also RESTRICT2, MXBITONIC2

%   Author: Graham M. Treece, University of Cambridge, UK, Sept 2021.


% check for appropriate inputs
narginchk(1, 2);

% also check for appropriate outputs
nargoutchk(1, 1);

% other input checks
if ~(ndims(A) == 2) && ~(ndims(A) == 3)
  error('A must be a grey or colour image.');
end
if (nargin < 2)
  off = [0 0];
else
  off = 1-rem(sz, 2);
end

% set up a catmull rom restriction operator and extension vectors
h = [-1 0 9 16 9 0 -1]/16;
%h = [1 8 23 32 23 8 1]/48;
hh = length(h);
m = [ones(1,1) 1:size(A,1) size(A,1)*ones(1,1)];
n = [ones(1,1) 1:size(A,2) size(A,2)*ones(1,1)];
mm = (length(m)-1)*2 + hh;
nn = (length(n)-1)*2 + hh;

% filter and remove edge samples
for k=1:size(A,3)
  
  % This version does not need the signal processing toolbox
  b = zeros(mm, length(n));
  b(1:2:mm,:) = conv2(h(1:2:hh)', double(A(m,n,k)));
  b(2:2:mm,:) = conv2(h(2:2:hh)', double(A(m,n,k)));
  c = zeros(mm, nn);
  c(:,1:2:nn) = conv2(h(1:2:hh), b);
  c(:,2:2:nn) = conv2(h(2:2:hh), b);
  X(:,:,k) = c(6:(end-5-off(1)),6:(end-5-off(2)));
  
  % this version does
  %X(:,:,k) = upfirdn(upfirdn(double(A(m,n,k)), h, 2, 1)', h, 2, 1)';

end
%X = X(6:(end-5-off(1)),6:(end-5-off(2)),:);

% ensure the right data type
if strcmp(class(A), 'uint8')
  X = uint8(X);
elseif strcmp(class(A), 'int32')
  X = int32(X);
end