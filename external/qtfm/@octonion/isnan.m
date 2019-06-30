function tf = isnan(A)
% ISNAN  True for Not-a-Number.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

tf = isnan(A.a) | isnan(A.b);

% $Id: isnan.m 1004 2017-11-15 17:14:09Z sangwine $
