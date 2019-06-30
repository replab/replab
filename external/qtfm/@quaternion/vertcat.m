function c = vertcat(varargin)
% VERTCAT Vertical concatenation.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2009 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

if length(varargin) == 1
    c = varargin{1}; % We have been passed one argument, nothing to do!
    return
end

% This is implemented recursively for simplicity, as it is unlikely to be
% used for more than a few arguments.

a = quaternion(varargin{1}); % The call to the class constructor ensures
b = quaternion(varargin{2}); % that a and b are both quaternions, even if
                             % the first or second varargin parameter was
                             % something else. Calling the constructor
                             % exploits the error checks there, that would
                             % be complex to include here for handling rare
                             % problems like catenating strings with
                             % quaternions.
if isempty(b)
    if length(varargin) == 2
        c = a;
        return
    else
        c = vertcat(a, varargin{3:end});
        return
    end
elseif isempty(a)
    if length(varargin) == 2
        c = b;
        return
    else
        c = vertcat(b, varargin{3:end});
        return
    end
end

if isempty(a.w)
    if isempty(b.w)
        % a and b are both pure, so the result will be too.
        c = a;
        c.x = [a.x; b.x];
        c.y = [a.y; b.y];
        c.z = [a.z; b.z];
    else
        % a is pure, but b isn't, so the result will need to be full and we
        % have to supply zeros for the scalar part of a with the same type
        % as the elements of a.
        c = b; % Use b to make c, as we need c to be full.
        t = a.x;
        c.w = [zeros(size(a), class(t)); b.w];
        c.x = [t;   b.x];
        c.y = [a.y; b.y];
        c.z = [a.z; b.z];
    end
else
    if isempty(b.w)
        % b is pure, but a isn't, so the result will need to be full and we
        % have to supply zeros for the scalar part of a with the same type
        % as the elements of a.
        c = a; % Use a to make c, as we need c to be full.
        t = b.x;
        c.w = [a.w; zeros(size(b), class(t))];
        c.x = [a.x; t  ];
        c.y = [a.y; b.y];
        c.z = [a.z; b.z];
    else
        % a and b are both full quaternion arrays.
        c = a;
        c.w = [a.w; b.w];
        c.x = [a.x; b.x];
        c.y = [a.y; b.y];
        c.z = [a.z; b.z];
    end
end

if length(varargin) == 2
    return
else
    c = vertcat(c, varargin{3:end});
end

% TODO Given arrays of inconsistent dimensions to concatenate, this code
% currently doesn't trap the error but passes the arrays to Matlab/vertcat
% which does. It might be better to catch the error here.

% $Id: vertcat.m 1004 2017-11-15 17:14:09Z sangwine $

