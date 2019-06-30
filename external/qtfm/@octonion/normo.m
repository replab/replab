function a = normo(o)
% NORMO Norm of an octonion, the sum of the squares of the components.

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

a = normq(o.a) + normq(o.b);

% $Id: normo.m 1004 2017-11-15 17:14:09Z sangwine $
