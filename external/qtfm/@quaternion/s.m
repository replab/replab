function x = s(q)
% S(Q) Scalar part of a full quaternion.

% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% narginchk(1, 1), nargoutchk(0, 1)

x = q.w;

% $Id: s.m 1004 2017-11-15 17:14:09Z sangwine $
