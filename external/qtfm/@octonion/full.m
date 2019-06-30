function b = full(q)
% FULL  Convert sparse matrix to full matrix.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

b = overload(mfilename, q);

% $Id: full.m 1004 2017-11-15 17:14:09Z sangwine $
