function r = eq(a, b)
% ==  Equal.
% (Quaternion overloading of standard Matlab function.)
%
% If one of the operands is not a quaternion and the other has zero vector part,
% the result is obtained by comparing the non-quaternion operand with the scalar
% part of the quaternion operand.

% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1) 

if isa(a, 'quaternion') && isa(b, 'quaternion')
    
    pa = isempty(a.w);
    pb = isempty(b.w);
    if pa && pb
        r = a.x == b.x & a.y == b.y & a.z == b.z;
    elseif pa
        % a is pure, but b isn't, so they can only be equal if the
        % scalar part of b is zero and the vector parts are equal.
        
        r = b.w == 0 & vector(a) == vector(b);
    elseif pb
        % b is pure, but a isn't, so they can only be equal if the
        % scalar part of a is zero and the vector parts are equal.
        
        r = a.w == 0   & vector(a) == vector(b);
    else
        r = a.w == b.w & vector(a) == vector(b);
    end
    
else
    % One of the arguments is not a quaternion (the other must be, or
    % Matlab would not call this function). The non-quaternion argument
    % must be a numeric (if we don't impose this restriction it could be
    % anything such as a cell array or string which makes no sense at all
    % to compare with a quaternion).
            
    if isa(a, 'quaternion') && isa(b, 'numeric')
        if isempty(a.w)
            % a has no scalar part, and b has no vector part, since it is not
            % a quaternion. The result can only be true if b is zero, and a
            % has zero vector part.
            
            r = b == 0 & a == quaternion(0, 0, 0);
        else
            % a is a full quaternion, so the result can be true only if a has
            % zero vector part, and b is equal to the scalar part of a.
            
            r = a.w == b & vector(a) == quaternion(0, 0, 0);
        end
    elseif isa(a, 'numeric') && isa(b, 'quaternion')
        r = b == a; % Swap the order and compare them using the code above.
    else
        error('Cannot compare a quaternion with a non-numeric');    
    end
end

% $Id: eq.m 1004 2017-11-15 17:14:09Z sangwine $
