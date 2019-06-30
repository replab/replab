function q = isoquaternion(z, a)
% ISOQUATERNION  Construct a quaternion from a complex number, preserving
%                the modulus and argument, and using the axis of the second
%                argument as the axis of the result. If the second argument
%                has a null vector part, the complex value z is returned as
%                a complex scalar part (note 2).
%
% TODO If the quaternion a is pure, so should the result be pure. This
% would mean that sign(qi) would return qi and not a full quaternion.
%
% Copyright (c) 2006, 2008, 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

if ~isa(a, 'quaternion')
    error('Second argument to private function isoquaternion must be a quaternion.');
end;

if isa(z,'quaternion')
    error('First argument to private function isoquaternion must not be a quaternion.');
end;

if ~isnumeric(z)
    error('First argument to private function isoquaternion must be numeric.');
end

undefined = abs(normq(v(a))) < eps; % Copied from the axis function (q.v.).
  defined = ~undefined;
  
q = zerosq(size(z), class(z)); % Preconstruct the result array.

% Now we can construct the result. Any elements for which undefined is true
% are assumed to have been computed from a quaternion with null vector
% part. We reconstruct these by putting z into the scalar part, and supply
% a null vector part. Other elements use the axis of the second argument as
% the axis of the quaternion.

% This is a private class function. We have to use subsref, subsasgn and
% substruct to handle the indexing here, because normal indexing notation
% does not work in a class function.
% See the file 'implementation notes.txt', item 8, for more details.

% It is possible for the logical arrays re or er to be all zero, in which
% case we must not try to index into the complex or quaternion array, as
% the subsasgn will select no elements.
% TODO Consider how to handle arrays of more than 2 dimensions. One way
% would be to process the array as a 1-D array and then reshape it to the
% original size.

if any(any(undefined))
    q = subsasgn(q, substruct('()', {undefined}), ...
        quaternion(z(undefined))); % A zero vector part will be supplied.
end;
if any(any(defined))
    q = subsasgn(q, substruct('()', {defined}), ...
        quaternion(real(z(defined)), ...
                   imag(z(defined)) .* axis(subsref(a, substruct('()', {defined})))));
end

% Note 1: the complex argument z may be in any of the four quadrants of the
% plane, and so may the quaternion result. This means that if the axis is
% extracted from the quaternion result, it may point in the opposite
% direction to the axis of the second argument, a.

% Note 2: if the second argument has a null vector part, it is a real or
% complex number. Therefore we return a complex number (a quaternion with
% null vector part). This can arise, for example, if we take the square
% root of a quaternion with a negative scalar part, and a null vector part.
% A classic example is sqrt(quaternion(-1,0,0,0)) which should return 0+1i.

% $Id: isoquaternion.m 1004 2017-11-15 17:14:09Z sangwine $
