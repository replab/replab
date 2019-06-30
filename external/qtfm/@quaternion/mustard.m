function r = mustard(f, g, mu, H)
% MUSTARD convolution. This is the convolution equivalent to pointwise
% multiplication of the Fourier transforms of f and g. It is dependent on
% the definition of the Fourier transform. This function assumes a
% one-sided QFT as computed by QFFT and QFFT2.

% This code handles one dimensional or two dimensional arrays f and g
% according to the parameters supplied. If they are vectors, they must be
% both row or both column vectors.

% Copyright (c) 2013-2016 Stephen J. Sangwine and Todd A. Ell.
% See the file : Copyright.m for further details.

% References:
%
% David Mustard, 'Fractional convolution',
% The Journal of the Australian Mathematical Society,
% Series B, Applied Mathematics, Vol 40, 257­265, 1998.
% doi:10.1017/ S0334270000012509
%
% De Bie, H. and De Schepper, N. and Ell, T. A. and Rubrecht, K. and
% Sangwine, S. J., 'Connecting spatial and frequency domains for the
% quaternion Fourier transform', Applied Mathematics and Computation, 271,
% 581-593, 15 November 2015, doi:10.1016/j.amc.2015.09.045. 

narginchk(3, 4), nargoutchk(0, 1)

% Parameters:
%
% f and g - quaternion vectors or matrices to be convolved.
% mu      - axis (of the QFT if the convolution is implemented pointwise in
%           the Fourier domain).
% H       - the transform axis, 'L' or 'R', corresponding to the parameters
%           of the one-sided QFFT functions.
%
% TODO Add the case of a two-sided QFFT, perhaps by allowing H to be a unit
% pure quaternion, rather than a character.

if (isvector(f) ~= isvector(g)) || ~ismatrix(f) || ~ismatrix(g)
    % Note that ismatrix tests for more than 2 dimensions, and returns true
    % for vectors. The first test checks that f and g are both vectors or
    % both not vectors, whereas the second and third tests trap cases where
    % f or g has more than two dimensions.
    error('First and second parameters must both be vectors or both matrices.')
end

if (isrow(f) && iscolumn(g)) || (iscolumn(f) && isrow(g))
    % Note that if f and g are matrices, isrow and iscolumn yield false,
    % so this error message will not be triggered.
    error('First and second parameters must be vectors of the same type (row/column)')
end

if nargin == 3
    % Supply a default value for the H parameter. This is 'L' because the
    % QFFT functions have the same default.
    H = 'L';
end

if ~ispure(mu)
    error('Third parameter must be a pure quaternion')
end
mu = unit(mu); % Ensure mu is a unit pure quaternion.

if nargin == 4
    if ~ischar(H) || length(H) > 1
        error('Fourth parameter must be a single character.')
    end
    if H ~= 'L' && H ~= 'R'
        error('Fourth parameter must be L or R.')
    end
end

% Select the appropriate convolution function and construct a handle to it.

if isvector(f)
    C = str2func('convw');
    F = str2func('fftflip');
else
    C = str2func('convw2');
    F = str2func('fftflip2');
end

switch H
    case 'L'
        [fp, fm] = opd(f, mu);
        r = C(fm, g) + C(fp, F(g));
    case 'R'
        [gp, gm] = opd(g, mu);
        r = C(f, gm) + C(F(f), gp);
    otherwise
        error('Fourth parameter has bad value, program error.')
end

end

% $Id: mustard.m 1004 2017-11-15 17:14:09Z sangwine $
