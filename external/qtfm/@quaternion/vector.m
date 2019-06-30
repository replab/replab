function p = vector(q)
% VECTOR   Quaternion vector part. Synonym of V.

% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

p = q; p.w = []; % Copy and then set the scalar part to empty.

% $Id: vector.m 1004 2017-11-15 17:14:09Z sangwine $
