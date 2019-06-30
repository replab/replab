function r = cat(dim, varargin)
% CAT Concatenate arrays.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2009, 2010, 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(3, inf), nargoutchk(0, 1)

if ~isnumeric(dim)
    error('First parameter must be numeric.')
end

a = quaternion(varargin{1}); % The call to the class constructor ensures
b = quaternion(varargin{2}); % that a and b are both quaternions, even if
                             % the first or second varargin parameter was
                             % something else. Calling the constructor
                             % exploits the error checks there, that would
                             % be complex to include here for handling rare
                             % problems like catenating strings with
                             % quaternions.
if isempty(a.w)
    if isempty(b.w)
        % a and b are both pure, so the result will be too.
        r = quaternion(cat(dim, a.x, b.x), ...
                       cat(dim, a.y, b.y), ...
                       cat(dim, a.z, b.z));
    else
        % a is pure, but b isn't, so the result will need to be full and we
        % have to supply zeros for the scalar part of a with the same type
        % as the elements of a.
        t = a.x;
        r = quaternion(cat(dim, zeros(size(a), class(t)), b.w),...
                       cat(dim, t,   b.x), ...
                       cat(dim, a.y, b.y), ...
                       cat(dim, a.z, b.z));
    end
else
    if isempty(b.w)
        % b is pure, but a isn't, so the result will need to be full and we
        % have to supply zeros for the scalar part of a with the same type
        % as the elements of a.
        t = b.x;
        r = quaternion(cat(dim, a.w, zeros(size(b), class(t))),...
                       cat(dim, a.x, t  ), ...
                       cat(dim, a.y, b.y), ...
                       cat(dim, a.z, b.z));
    else
        % a and b are both full quaternion arrays.
        r = quaternion(cat(dim, a.w, b.w), ...
                       cat(dim, a.x, b.x), ...
                       cat(dim, a.y, b.y), ...
                       cat(dim, a.z, b.z));
    end
end

if nargin > 3
    r = cat(dim, r, varargin{3:end});    
end

% $Id: cat.m 1004 2017-11-15 17:14:09Z sangwine $
