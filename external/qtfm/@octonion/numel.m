function n = numel(o, varargin)
% NUMEL   Number of elements in an array or subscripted array expression.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2011, 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

if nargin == 1
    % The two components of the octonion cannot have different sizes, so we
    % just pass the first of the two to the quaternion numel.
    
    n = numel(o.a);
else
    n = numel(o.a, varargin);
end
    
end

% $Id: numel.m 1004 2017-11-15 17:14:09Z sangwine $
