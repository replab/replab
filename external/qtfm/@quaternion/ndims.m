function n = ndims(a)
% NDIMS   Number of array dimensions.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2008, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

n = ndims(a.x);

% An alternative, previously used is:
%
% n = length(size(q));
%
% but using the builtin function on the x-component is simpler, c.f. the
% end function, where the same approach is used.

% $Id: ndims.m 1004 2017-11-15 17:14:09Z sangwine $
