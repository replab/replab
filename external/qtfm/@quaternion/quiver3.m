function h = quiver3(O, P, varargin)
% QUIVER3  Display pure quaternion array as 3D quiver plot
% (Quaternion overloading of standard Matlab function.)
%
% Takes the same parameters as the Matlab function of the same name, except
% that the first six parameters (X, Y, Z, U, V, W) are replaced by two
% quaternion parameters, which must be pure quaternion vectors. For more
% sophisticated requirements, use the Matlab function and pass the X, Y and
% Z components of the quaternion array as the first three parameters.

% Copyright (c) 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

nargoutchk(0, 1)

for k = 1:length(varargin)
    if isa(varargin{k}, 'quaternion')
        error('Only the first two parameters are permitted to be quaternions.')
    end
end

if ~isempty(O.w)
    % The first parameter must always be a quaternion array, so check it
    % now.
    error('Quaternion array in first parameter must be pure.')
end

if nargin == 1        
    % Since there is only one parameter, we assume zero for the origin.
    
    O = O(:); % This is a Matlab trick to convert an array of any number of
              % dimensions into a vector.
              
    Z = zeros(size(O));
    if nargout == 0
        quiver3(Z, Z, Z, O.x, O.y, O.z);
    else
        h = quiver3(Z, Z, Z, O.x, O.y, O.z);
    end
    return
end

% From here there must be at least two input parameters. These may both be
% quaternions or there may be one quaternion parameter followed by
% varargin. (The second parameter may be quaternion or not (it could be a
% string such as 'LineSpec').)

if isa(P, 'quaternion')
    
    if ~isempty(P.w)
        error('Quaternion array in second parameter must be pure.')
    end
    
    if any(size(O) ~= size(P))
        error('The two quaternion arrays must have the same size.');
    end

    if ~isvector(O)
        O = O(:); % If O isn't a vector P isn't either, because we checked
        P = P(:); % that the sizes matched just above.
    end
    
    if isempty(varargin)
        if nargout == 0
            quiver3(O.x, O.y, O.z, P.x, P.y, P.z);
        else
            h = quiver3(O.x, O.y, O.z, P.x, P.y, P.z);
        end
    else
        if nargout == 0
            quiver3(O.x, O.y, O.z, P.x, P.y, P.z, varargin{:});
        else
            h = quiver3(O.x, O.y, O.z, P.x, P.y, P.z, varargin{:});
        end
    end
else
    % The second parameter is not a quaternion, so we must pass it to the
    % Matlab quiver3 function, and assume it is valid (if it isn't quiver3
    % will raise an error). The simplest way to handle this is a recursive
    % call in which we supply a new first parameter. The code above then
    % does the rest.
              
    Z = zerosv(size(O));
    if isempty(varargin)
        if nargout == 0
            quiver3(Z, O, P);
        else
            h = quiver3(Z, O, P);
        end
    else
        if nargout == 0
            quiver3(Z, O, P, varargin{:});
        else
            h = quiver3(Z, O, P, varargin{:});
        end
    end
end

% TODO Consider how to include RGB colouring cf scatter3.m.

% $Id: quiver3.m 1004 2017-11-15 17:14:09Z sangwine $
