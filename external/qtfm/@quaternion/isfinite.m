function tf = isfinite(A)
% ISFINITE  True for finite elements.  
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isempty(A.w)
    tf = isfinite(A.x) & isfinite(A.y) & isfinite(A.z);
else
    tf = isfinite(A.w) & isfinite(vee(A));
end

% $Id: isfinite.m 1004 2017-11-15 17:14:09Z sangwine $
