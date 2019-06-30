function n = length(o)
% LENGTH   Length of vector.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

n = length(o.a); % This calls the quaternion length function.

% $Id: length.m 1004 2017-11-15 17:14:09Z sangwine $
