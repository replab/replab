function a = reshape(q, varargin)
% RESHAPE Change size.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

nargoutchk(0, 1)

a = overload(mfilename, q, varargin{:});

% $Id: reshape.m 1004 2017-11-15 17:14:09Z sangwine $
