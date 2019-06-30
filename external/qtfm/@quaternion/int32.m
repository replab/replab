function d = int32(q)
% INT32 Convert to signed 32-bit integer (obsolete).
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

error(['Conversion to int32 from quaternion is not possible. ',...
       'Try cast(q, ''int32'')'])

% Note: this function was replaced from version 0.9 with the convert
% function, because it is incorrect to provide a conversion function
% that returns a quaternion result.

% $Id: int32.m 1004 2017-11-15 17:14:09Z sangwine $

