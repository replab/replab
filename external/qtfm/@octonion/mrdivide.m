function r = mrdivide(a, b)
% /   Slash or right matrix divide.
% (Octonion overloading of standard Matlab function.)
%
% This function is implemented only for the case in which the
% second parameter is a scalar. In this case, the result is as
% given by the ./ function.   The reason for implementing this
% function is that some Matlab code uses / when ./ should have
% been used.    See, for example, the function cov.m in Matlab
% versions up to (at least) 7.0.4.365 (R14) Service Pack 2.

% Copyright (c) 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

if ~isscalar(b)
    error('octonion/mrdivide is implemented only for division by a scalar.');
end

r = a ./ b;

% $Id: mrdivide.m 1004 2017-11-15 17:14:09Z sangwine $
