function B = permute(A, order)
% PERMUTE Rearrange dimensions of N-D array
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

B = overload(mfilename, A, order);

% $Id: permute.m 1004 2017-11-15 17:14:09Z sangwine $

