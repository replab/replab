function h = scatter3(Q, varargin)
% SCATTER3  Display pure quaternion array as 3D scatter plot
% (Quaternion overloading of standard Matlab function.)
%
% Takes the same parameters as the Matlab function of the same name, except
% that the first three parameters (X, Y, Z) are replaced by a single
% quaternion parameter, which must be a pure quaternion vector. If two
% parameters are given and the second is 'rgb' or 'seq', a coloured scatter
% plot is produced. 'rgb' colours the plotted points according to the
% direction of the position vector of the point from the origin using the
% convention of X = red, Y = green and Z = blue (the standard RGB colour
% space used in image processing). This is useful when the spatial position
% of points is of interest. 'seq' colours the points in a graded sequence
% of progressively darker values matching the sequence of points in the
% quaternion vector Q. This is useful when the ordering of the points is of
% interest. If more than one varargin parameter is given, or the first is
% neither 'rgb' nor 'seq', the varargin parameters are passed to the Matlab
% scatter3 function, and must therefore conform to the requirements of that
% function. For more sophisticated requirements, use the Matlab function
% and pass the X, Y and Z components of the quaternion array as the first
% three parameters.

% TODO Consider adding output parameters as per scatter4p3.

% Copyright (c) 2009, 2010, 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

nargoutchk(0, 1)

for k = 1:length(varargin)
    if isa(varargin{k}, 'quaternion')
        error('Only the first parameter is permitted to be a quaternion.')
    end
end

if ~isvector(Q)
    Q = Q(:); % This is a Matlab trick to convert an array of any number of
              % dimensions into a vector. Previous versions of this
              % function implemented an error check here and refused to
              % handle a non-vector input array.
end

if ~isempty(Q.w)
    warning(['Quaternion array is not pure, ignoring scalar part.', ...
             ' Largest scalar part is ', num2str(max(abs(scalar(Q))))]);
    Q = vector(Q);
end

if isempty(varargin)
    % Default case, no colouring other than what the standard Matlab
    % scatter3 gives.
    if nargout == 0
        scatter3(Q.x, Q.y, Q.z);
    else
        h = scatter3(Q.x, Q.y, Q.z);
    end
    return
end

v1 = varargin{1};

if length(varargin) == 1 && ischar(v1)
    
    % There is one varargin parameter and it is a character string, so we
    % need to check it for a possible parameter 'rgb' or 'seq'.
    
    if strcmp(v1, 'rgb')
        % RGB colouring.
        C = (unit(Q) + quaternion(1, 1, 1))./2;
        if size(C, 2) ~= 1, C = C.'; end % Make sure C is a column vector.
        if nargout == 0
            scatter3(Q.x, Q.y, Q.z, [], [C.x C.y C.z]); % Note 1.
        else
            h = scatter3(Q.x, Q.y, Q.z, [], [C.x C.y C.z]); % Note 1.
        end
        return
    end
    
    if strcmp(v1, 'seq')
        % Sequential colouring.
        if nargout == 0
            scatter3(Q.x, Q.y, Q.z, [], colormap(copper(length(Q)))); % Note 1.
        else
            h = scatter3(Q.x, Q.y, Q.z, [], colormap(copper(length(Q)))); % Note 1.
        end
        return
    end
end

% If we reach here, either:
%
% 1. The first varargin is a string, but it isn't one of the recognised
%    values, so we must pass it to Matlab. Of course, if it is a
%    mis-spelling of 'rgb' or 'seq' this will lead to an error in the
%    Matlab function, but to avoid that here would be too difficult.
%
% 2. varargin is either of length > 1, or the length == 1, but the first
%    varargin parameter is not a string. So we pass the varargin to the
%    Matlab function.

if nargout == 0
    scatter3(Q.x, Q.y, Q.z, varargin{:});
else
    h = scatter3(Q.x, Q.y, Q.z, varargin{:});
end

end

% Note 1. The empty fourth parameter is in place of S (see the Matlab HELP
% for scatter3. Without it, Matlab's scatter3 tries to interpret C as S and
% fails with an error.

% $Id: scatter3.m 1004 2017-11-15 17:14:09Z sangwine $
