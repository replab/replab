function b = repmat(a, m, n)
% REPMAT Replicate and tile an array.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 3), nargoutchk(0, 1) 

if nargin == 2
    b = overload(mfilename, a, m);
else
    b = overload(mfilename, a, m, n);
end

% $Id: repmat.m 1004 2017-11-15 17:14:09Z sangwine $
