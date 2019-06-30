function X = sqrtm(A)
% Matrix square root.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2008, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

X = overloadm(mfilename, A);

% TODO Implement a more accurate dedicated algorithm for this function.

% $Id: sqrtm.m 1004 2017-11-15 17:14:09Z sangwine $
