function u = uminus(a)
% -  Unary minus.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

u = a;
if ~isempty(u.w)
    u.w = -u.w;
end
u.x = -u.x;
u.y = -u.y;
u.z = -u.z;

% $Id: uminus.m 1004 2017-11-15 17:14:09Z sangwine $
