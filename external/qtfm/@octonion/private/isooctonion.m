function q = isooctonion(z, a)
% ISOOCTONION  Construct an octonion from a complex number, preserving
%              the modulus and argument, and using the axis of the second
%              argument as the axis of the result. If the second argument
%              has a null vector part, the complex value z is returned as
%              a complex scalar part (note 2).
%
% TODO If the octonion a is pure, so should the result be pure. This
% would mean that sign(oi) would return oi and not a full octonion.
%
% Copyright (c) 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

if ~isa(a, 'octonion')
    error('Second argument to private function isooctonion must be an octonion.');
end;

if isa(z,'octonion')
    error('First argument to private function isooctonion must not be an octonion.');
end;

if ~isnumeric(z)
    error('First argument to private function isooctonion must be numeric.');
end

undefined = abs(normo(v(a))) < eps; % Copied from the axis function (q.v.).
  defined = ~undefined;
  
q = zeroso(size(z), class(z)); % Preconstruct the result array.

% Now we can construct the result. Any elements for which undefined is true
% are assumed to have been computed from an octonion with null vector
% part. We reconstruct these by putting z into the scalar part, and supply
% a null vector part. Other elements use the axis of the second argument as
% the axis of the octonion.

% This is a private class function. We have to use subsref, subsasgn and
% substruct to handle the indexing here, because normal indexing notation
% does not work in a class function.
% See the file 'implementation notes.txt', item 8, for more details.

% It is possible for the logical arrays re or er to be all zero, in which
% case we must not try to index into the complex or octonion array, as
% the subsasgn will select no elements.
% TODO Consider how to handle arrays of more than 2 dimensions. One way
% would be to process the array as a 1-D array and then reshape it to the
% original size.

if any(any(undefined))
    q = subsasgn(q, substruct('()', {undefined}), ...
        octonion(z(undefined))); % A zero vector part will be supplied.
end;
if any(any(defined))
    q = subsasgn(q, substruct('()', {defined}), ...
        octonion(real(z(defined)), ...
                 imag(z(defined)) .* axis(subsref(a, substruct('()', {defined})))));
end

% Note 1: the complex argument z may be in any of the four quadrants of the
% plane, and so may the octonion result. This means that if the axis is
% extracted from the octonion result, it may point in the opposite
% direction to the axis of the second argument, a. (This not is copied from
% the quaternion case, and should be checked for the octonion case.)

% Note 2: if the second argument has a null vector part, it is a real or
% complex number. Therefore we return a complex number (an octonion with
% null vector part). This can arise, for example, if we take the square
% root of an octonion with a negative scalar part, and a null vector part.
% A classic example is sqrt(octonion(-1,0,0,0)) which should return 0+1i.

% $Id: isooctonion.m 1004 2017-11-15 17:14:09Z sangwine $
