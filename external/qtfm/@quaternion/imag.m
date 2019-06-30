function p = imag(q)
% IMAG   Imaginary part of a quaternion.
% (Quaternion overloading of standard Matlab function.)
%
% This function returns the quaternion that is the imaginary
% part of q. If q is a real quaternion, it returns zero.

% Copyright (c) 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

p = overload(mfilename, q);

end

% $Id: imag.m 1004 2017-11-15 17:14:09Z sangwine $
