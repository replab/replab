function m = suppress(varargin)
% SUPPRESS - constructs a N-by-M vector used to suppress negative
% frequencies in constructing a discrete analytic or hyperanalytic signal.
% Parameters may be [N, M] or N, M as in Matlab functions such as zeros and
% ones to permit calling with suppress(size(f)) or suppress(1, N) etc. The
% result returned will be a row or column according to the parameters
% supplied.

% Copyright 2010-2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.
% Reference:
%
% S. Lawrence Marple, Jr.,
% "Computing the Discrete-Time 'Analytic' Signal via FFT",
% IEEE Transactions on Signal Processing, 47(9), September 1999, 2600-2603.
% DOI:10.1109/78.782222

nargoutchk(0, 1)

if length(varargin) > 2
    error('Too many parameters')
end

m = ones(varargin{:}); % Construct an initial result vector.

if ~isvector(m)
    error('Parameters must represent the size of a vector.')
end

% The result has the form [1, 2 ... 2, 1, 0 ... 0] for even length, and
% [1, 2 ... 2, 0 ... 0] for odd length, or its transpose, the number of
% twos and zeros being the same in both cases. (The omitted 1 in the
% odd-length case corresponds to the Nyquist frequency in the discrete
% Fourier transform.)

L = length(m);
K = floor((L - 1) ./2 ); % The number of zeros and twos in the result.

m(2:K+1)   = 2;
m(L-K+1:L) = 0;

end

% $Id: suppress.m 1004 2017-11-15 17:14:09Z sangwine $