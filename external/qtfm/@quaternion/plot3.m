function h = plot3(Q, varargin)
% PLOT3  Display pure quaternion array as 3D line plot
% (Quaternion overloading of standard Matlab function.)
%
% Takes the same parameters as the Matlab function of the same name, except
% that the first three parameters (X, Y, Z) are replaced by a single
% quaternion parameter, which must be a pure quaternion vector. If varargin
% parameters are given they are passed to the Matlab plot3 function, and
% must therefore conform to the requirements of that function.

% TODO Handle more than one quaternion parameter, so that multiple lines
% can be passed to plot3, or consider whether a matrix should be plotted by
% columns?

% Copyright (c) 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

nargoutchk(0, 1)

for k = 1:length(varargin)
    if isa(varargin{k}, 'quaternion')
        error('Only the first parameter is permitted to be a quaternion.')
    end
end

if ~isempty(Q.w)
    error('Quaternion parameter must be pure.')
end

if ~isvector(Q)
    error('The first parameter must be a vector.')
end

if isempty(varargin)
    if nargout == 0
        plot3(Q.x, Q.y, Q.z);
    else
        h = plot3(Q.x, Q.y, Q.z);
    end
else
    if nargout == 0
        plot3(Q.x, Q.y, Q.z, varargin{:});
    else
        h = plot3(Q.x, Q.y, Q.z, varargin{:});
    end
end

end

% $Id: plot3.m 1004 2017-11-15 17:14:09Z sangwine $
