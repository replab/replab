function b = cast(q, newclass)
% CAST  Cast quaternion variable to different data type.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

if ~ischar(newclass)
    error('Second parameter must be a string.')
end

b = overload(mfilename, q, newclass);

% $Id: cast.m 1004 2017-11-15 17:14:09Z sangwine $
