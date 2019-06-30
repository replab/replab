function Y = diff(X, n, dim)
% DIFF   Differences.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 3), nargoutchk(0, 1) 

if nargin == 1
    Y = overload(mfilename, X);
elseif nargin == 2
    Y = overload(mfilename, X, n);
else
    if ~isnumeric(dim)
        error('Dimension argument must be numeric');
    end
    
    if ~isscalar(dim) || ~ismember(dim, 1:ndims(X))
        error(['Dimension argument must be a positive'...
               ' integer scalar within indexing range.']);
    end

    Y = overload(mfilename, X, n, dim);   
end

% Note: In the absence of this file, the Matlab function of the same name is
% called for a quaternion array, but it raises an error. That is why this
% function was written, and why the Matlab code can't simply be used to
% difference quaternion arrays.

% $Id: diff.m 1004 2017-11-15 17:14:09Z sangwine $

