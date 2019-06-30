function P = nu(A)
% NU Computes the nu (right) adjoint of the octonion matrix A as in:
%
% Yongge Tian, 'Matrix representations of octonions and their applications'
% Advances in Applied Clifford Algebras, 10 (1), pp.61-90, 2000.
% [Also available in preprint form as arXiv:math/0003166, 1 April 2000.]

% Copyright (c) 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Since this is a private function (at least for the moment) we omit error
% checks.

a = A.a; % Split the octonion matrix into its quaternion Cayley-Dickson
b = A.b; % components.

P = [tao(a), -phi(conj(b));
     phi(b), tao(conj(a))];
 
if isscalar(A)
    return;
end
 
% The matrix A is not scalar, so we need to re-arrange P so that all the
% real values for each octonion are gathered into a 8-by-8 block. We do
% this in two steps, once to re-arrange the rows, and the second to
% re-arrange the columns, making use of an intermediate array T.

[R, C] = size(A);

T = zeros(size(P));

T(1:8:end-7, :) = P(        1:1 * R, :); % This could be expressed as a for
T(2:8:end-6, :) = P(1 * R + 1:2 * R, :); % loop, but it is clearer and
T(3:8:end-5, :) = P(2 * R + 1:3 * R, :); % probably faster as inline code.
T(4:8:end-4, :) = P(3 * R + 1:4 * R, :);
T(5:8:end-3, :) = P(4 * R + 1:5 * R, :);
T(6:8:end-2, :) = P(5 * R + 1:6 * R, :);
T(7:8:end-1, :) = P(6 * R + 1:7 * R, :);
T(8:8:end,   :) = P(7 * R + 1:8 * R, :);

P(:, 1:8:end-7) = T(:,         1:1 * C);
P(:, 2:8:end-6) = T(:, 1 * C + 1:2 * C);
P(:, 3:8:end-5) = T(:, 2 * C + 1:3 * C);
P(:, 4:8:end-4) = T(:, 3 * C + 1:4 * C);
P(:, 5:8:end-3) = T(:, 4 * C + 1:5 * C);
P(:, 6:8:end-2) = T(:, 5 * C + 1:6 * C);
P(:, 7:8:end-1) = T(:, 6 * C + 1:7 * C);
P(:, 8:8:end  ) = T(:, 7 * C + 1:8 * C);

end

% $Id: nu.m 1004 2017-11-15 17:14:09Z sangwine $
