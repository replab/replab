function d = eps(X)
% EPS Floating-point relative accuracy.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% This function will only be called for the profile eps(X) where X is a
% quaternion. It is implemented only because some Matlab functions call it,
% and if not implemented these functions will not work for quaternions.

narginchk(1, 1), nargoutchk(0, 1) 

d = eps(abs(X)); % X is a quaternion, abs(X) is real or complex, and the
                 % eps function called from this line will be the built-in.
end

% $Id: eps.m 1004 2017-11-15 17:14:09Z sangwine $
