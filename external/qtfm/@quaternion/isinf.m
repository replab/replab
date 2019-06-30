function tf = isinf(A)
% ISINF  True for infinite elements.  
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isempty(A.w)
    tf = isinf(A.x) | isinf(A.y) | isinf(A.z);
else
    tf = isinf(A.w) | isinf(vee(A));
end

% $Id: isinf.m 1004 2017-11-15 17:14:09Z sangwine $
