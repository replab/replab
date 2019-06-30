function tf = isfinite(A)
% ISFINITE  True for finite elements.  
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

tf = isfinite(A.a) & isfinite(A.b);

end

% $Id: isfinite.m 1004 2017-11-15 17:14:09Z sangwine $
