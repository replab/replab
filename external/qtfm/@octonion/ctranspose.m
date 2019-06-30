function t = ctranspose(a)
% '   Octonion conjugate transpose.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

t = conj(transpose(a));

% $Id: ctranspose.m 1004 2017-11-15 17:14:09Z sangwine $
