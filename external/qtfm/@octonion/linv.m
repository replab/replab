function X = linv(A, B)
% LINV   Left inverse operator of an octonion matrix. The result, if it
% exists, is the matrix that satisfies X * A = B. If B is omitted it
% defaults to an identity matrix, and X will then be the left inverse of A.

% Reference:
%
% Yongge Tian, 'Matrix representations of octonions and their applications'
% Advances in Applied Clifford Algebras, 10 (1), pp.61-90, 2000.
% [Also available in preprint form as arXiv:math/0003166, 1 April 2000.]
% [Theorem 4.13 in the arXiv version, but note the error explained in Note
% 1 below.]

% Copyright (c) 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1)

if isscalar(A) || (nargin == 2 && isscalar(B))
    error('Octonion LINV function does not accept scalar argument(s).')
end

% TODO There is currently no check that the inverse exists, because it is
% not clear how to detect this condition. Clearly a check should be added.

[R, C] = size(A);

if R ~= C
    error('Matrix must be square.');
end

if nargin == 2 && any(size(B) ~= [R, C])
    error('Second parameter must have the same size as the first.')
end

% Compute the coefficients of the characteristic polynomial of the left
% adjoint of A. MATLAB orders this with the highest power first, whereas we
% need the lowest power first, so we flip it.

P = fliplr(poly(omega(A)));

r0 = P(1); P(1) = []; % Get the zeroth coefficient and remove it from P.
r1 = P(1); P(1) = []; % Ditto for the first-order coefficient.

if nargin == 2
    PA = B; % This will be used to compute powers of A, by successive
            % multiplication on the right, at each step of the iteration.
else
    PA = eye(R); % B was not supplied, so we use an identity matrix.
end

X = r1 .* PA;

% Now work along the polynomial accumulating the powers of A scaled by the
% remaining coefficients. The final coefficient is unity, but it is simpler
% to handle it in the loop than to eliminate the multiplication by 1.

for r = P
    PA = PA * A;
    X = X + r .* PA;
end

X = (-1/r0) .* X;

end

% Note 1. Tian's formula in Theorem 4.13 is an incorrect transcription of
% his equation (4.20). Reading from the right, the powers of A have one too
% many multiplications. That is, r2 should scale A and not A^2, as is
% evident from a comparison with (4.20).

% $Id: linv.m 1004 2017-11-15 17:14:09Z sangwine $
