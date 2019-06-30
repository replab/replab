function p = imag(o)
% IMAG   Imaginary part of an octonion.
% (Octonion overloading of standard Matlab function.)
%
% This function returns the octonion that is the imaginary
% part of o. If o is a real octonion, it returns zero.

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

p = overload(mfilename, o);

% $Id: imag.m 1004 2017-11-15 17:14:09Z sangwine $
