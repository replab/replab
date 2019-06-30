function tf = isnan(A)
% ISNAN  True for Not-a-Number.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isempty(A.w)
    tf = isnan(A.x) | isnan(A.y) | isnan(A.z);
else
    tf = isnan(A.w) | isnan(vee(A));
end

% $Id: isnan.m 1004 2017-11-15 17:14:09Z sangwine $
