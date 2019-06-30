function Y = iqdft(X, A, L)
% IQDFT Inverse discrete quaternion Fourier transform.
%
% This function computes the inverse discrete quaternion Fourier transform
% of X. See the function qdft.m for details. Because this is an inverse
% transform it divides (columns of) the result by N, the length of the
% transform.

% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(3, 3), nargoutchk(0, 1)

% We omit any check on the A and L parameters here, because the qdft
% function does this and there is no need to duplicate it here. In the
% unlikely event that -A has no meaning, an error will arise here.

if isvector(X)
    N = length(X);
else
    N = size(X, 1);
end

% We can compute the inverse transform by using the forward transform code
% with a negated axis. 

Y = qdft(X, -A, L) ./ N;

% $Id: iqdft.m 1004 2017-11-15 17:14:09Z sangwine $
