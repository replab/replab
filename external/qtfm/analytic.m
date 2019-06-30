function a = analytic(f)
% ANALYTIC   Constructs the analytic signal, for the case of a real vector.
%            The result is complex. This is not a quaternion function. It
%            is provided as part of QTFM for comparison with HYPERANALYTIC.
%
% Copyright 2007-2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Note: this function computes the same result as the Matlab function
% HILBERT available in the Signal Processing toolbox.

narginchk(1, 1), nargoutchk(0, 1)

if ~isvector(f) || isscalar(f) || ~isnumeric(f)
    error('Input parameter must be a real vector.');
end

if ~isreal(f)
    warning('Input parameter is complex. Carrying on regardless....')
end

a = ifft(fft(f) .* suppress(size(f)));

% $Id: analytic.m 1004 2017-11-15 17:14:09Z sangwine $
