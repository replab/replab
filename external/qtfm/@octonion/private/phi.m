function P = phi(A)
% PHI Computes the phi adjoint of the quaternion matrix A as in:
%
% Yongge Tian, 'Matrix representations of octonions and their applications'
% Advances in Applied Clifford Algebras, 10 (1), pp.61-90, 2000.
% [Also available in preprint form as arXiv:math/0003166, 1 April 2000.]
%
% Note: this function takes a quaternion argument, even though it is a
% private function of the octonion class. This is because of its use in the
% construction of octonion adjoint matrices. In addition, Tian's paper
% defines PHI for a single quaternion, but here we extend the idea to a
% matrix in which each quaternion is represented by a 4-by-4 block.

% Copyright (c) 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Since this is a private function (at least for the moment) we omit error
% checks.

P = [A.w, -A.x, -A.y, -A.z; ...
     A.x,  A.w, -A.z,  A.y; ...
     A.y,  A.z,  A.w, -A.x; ...
     A.z, -A.y,  A.x,  A.w];
 
if isscalar(A)
    return;
end
 
% The matrix A is not scalar, so we need to re-arrange P so that all the
% real values for each quaternion are gathered into a 4-by-4 block. We do
% this in two steps, once to re-arrange the rows, and the second to
% re-arrange the columns, making use of an intermediate array T.

[R, C] = size(A);

T = zeros(size(P));

T(1:4:end-3, :) = P(        1:1 * R, :);
T(2:4:end-2, :) = P(1 * R + 1:2 * R, :);
T(3:4:end-1, :) = P(2 * R + 1:3 * R, :);
T(4:4:end,   :) = P(3 * R + 1:4 * R, :);

P(:, 1:4:end-3) = T(:,         1:1 * C);
P(:, 2:4:end-2) = T(:, 1 * C + 1:2 * C);
P(:, 3:4:end-1) = T(:, 2 * C + 1:3 * C);
P(:, 4:4:end  ) = T(:, 3 * C + 1:4 * C);

% $Id: phi.m 1004 2017-11-15 17:14:09Z sangwine $
