function q = dc(A, B)
% DC: Inverse of the Cayley-Dickson function (cd) which returns the two
% components of a quaternion or octonion. This function constructs a
% quaternion from two (complex) numbers, or an octonion from two
% quaternions, or an octonion from a numeric and a quaternion.

% Copyright (c) 2008, 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

if nargin == 1
    
    if ~isnumeric(A) && ~isa(A, 'quaternion')
        error('Parameter must be numeric or quaternion.')
    end

    if isa(A, 'quaternion')
        z = zeros(size(A), class(A.x));
        q = octonion(scalar(A), A.x, A.y, A.z, z, z, z, z);
    else % A must be numeric since we checked it above.
        q = real(A) + imag(A) .* qi;
    end
    
else
    
    % There must be two parameters.

    if any(size(A) ~= size(B))
        error('Parameters must be the same size.')
    end

    if (~isnumeric(A) && ~isa(A, 'quaternion')) || ...
       (~isnumeric(B) && ~isa(B, 'quaternion'))
        error('Parameters must be numeric or quaternion.')
    end
    
    if isnumeric(A) && isnumeric(B)
        % Both parameters are numeric, so we assume we are constructing a
        % quaternion from two complex numbers. Other possibilities must be
        % dealt with explicitly by the user.
        q = quaternion(real(A), imag(A), real(B), imag(B));
        return
    end
    
    if isa(A, 'quaternion') && isa(B, 'quaternion')
        % The result must be an octonion.
        if ~strcmp(class(A.x), class(B.x))
            error('Cannot create octonion from two quaternions of different classes')
        end
        q = octonion(scalar(A), A.x, A.y, A.z, scalar(B), B.x, B.y, B.z);
    else
        % Handle all the remaining cases (mixed numeric/quaternion) by
        % promoting both parameters to quaternion, and calling this
        % function recursively.
        q = dc(quaternion(A), quaternion(B));
    end
    
end

% $Id: dc.m 1004 2017-11-15 17:14:09Z sangwine $
