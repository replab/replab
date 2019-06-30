function tf = isunitary(A, tol)
% ISUNITARY  True if the given matrix is unitary to within the tolerance
% given (optionally) by the second parameter.

% Copyright (c) 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

if nargin == 1
    tol = 16 .* eps; % The tolerance was not specified, supply a default.
end

[r, c] = size(A);

if r ~= c
    error('A non-square matrix cannot be unitary.');
end

% The method used is to subtract a quaternion identity matrix from A * A'.
% The result should be almost zero.  To compare it against the tolerance,
% we add the moduli of the four components.  This is guaranteed to give a
% real result, even when A is a complexified quaternion matrix.

D = A * A' - eyeq(r);

tf = all(all(abs(s(D)) + abs(x(D)) + abs(y(D)) + abs(z(D)) < tol));

% $Id: isunitary.m 1004 2017-11-15 17:14:09Z sangwine $
