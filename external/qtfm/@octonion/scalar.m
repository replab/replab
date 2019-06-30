function r = scalar(o)
% SCALAR   Octonion scalar part.
%
% This function returns zero in the case of pure octonions,
% whereas the function s gives empty if o is pure.

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

r = scalar(o.a); % The scalar part is in the 'a' quaternion.

% $Id: scalar.m 1004 2017-11-15 17:14:09Z sangwine $
