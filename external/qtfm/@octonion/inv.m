function R = inv(a)
% INV   Inverse of an octonion (matrix).
% (Octonion overloading of standard Matlab function.)

% TODO Copied from quaternion/inv with small change from normq to normo.
% EXPERIMENTAL AND GIVES INCORRECT RESULT! In any case, an octonion matrix
% has left and right inverses, so this function will probably be denoted as
% unimplemented, with a suitable warning. A better option may be to add L
% and R parameters, and give a warning if neither is specified. A third
% option would be to use [L, R] as the result, and compute both. Calling
% code can use [L, ~] or [~, R] to compute one only.

% Copyright (c) 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

[r, c] = size(a);

if r == 1 && c == 1
    R = conj(a) ./ normo(a);
else
    error(['Matrix inverse is not implemented for octonion matrices.', ...
           ' See the linv and rinv functions for left and right inverses.'])
end

% $Id: inv.m 1004 2017-11-15 17:14:09Z sangwine $
