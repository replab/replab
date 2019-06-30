function [P, B] = tridiagonalize(A)
% Tridiagonalize Hermitian matrix A, such that P * A * P' = B and
% P' * B * P = A. B is real, P is unitary, and B has the same eigenvalues
% as A.

% Copyright (c) 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(2, 2) 

% NB: This is a reference implementation. It uses explicit Householder
% matrices. Golub and van Loan, 'Matrix Computations', 2e, 1989, section
% 5.1.4 discusses efficient calculation of Householder transformations.
% This has been tried for bidiagonalization (see bidiagonalize2.m), but
% found to be slower than the use of explicit matrices as used here. This
% requires further study.

[r, c] = size(A);

if r ~= c
    error('Cannot tridiagonalize a non-square matrix');
end

if r < 2
    error('Cannot tridiagonalize a matrix smaller than 2 by 2.');    
end

if ~ishermitian(A)
    error('Matrix to be tridiagonalized is not (accurately) Hermitian');
end

[P, B] = internal_tridiagonalizer(A);

B = check(B); % Verify the result and convert to exactly tridiagonal real
              % form.

% -------------------------------------------------------------------------

function [P, B] = internal_tridiagonalizer(A)

[r, c] = size(A);

if r ~= c
    error('Internal tridiagonalize error - non-square matrix');
end

% TODO Replace this recursive code with iterative code in order to permit
% the handling of larger matrices. Cf bidiagonalize.m.

% Compute and apply a Householder transformation to the first row and
% column of A (these are conjugates of each other because A is Hermitian).
% We omit the first element of the first column in computing the
% Householder matrix, because it is already real (A is Hermitian). We apply
% P on both sides of A so that the first row and column are nulled out
% (apart from the first two elements in each case). Note that if B is 2 by
% 2, this is all we have to do.
    
P = quaternion(eye(r));
P(2 : end, 2 : end) = householder_matrix(A(2 : end, 1), eye(r - 1, 1));
B = P * A * P';

% Apply the algorithm recursively to sub-matrices of B, except when the
% sub-matrix is 2 by 2.

if r > 2
    Q = quaternion(eye(r));
    [Q(2 : end, 2 : end), B(2 : end, 2 : end)] = ...
       internal_tridiagonalizer(B(2 : end, 2 : end));
    P = Q * P;
end

% -------------------------------------------------------------------------

function R = check(B)

% Verify results, and convert the result to exactly tridiagonal form with
% no vector part.

D = [diag(B); diag(B, +1); diag(B, -1)]; % Extract the three diagonals.
O = tril(B, -2) + triu(B, +2);           % Extract the off-diagonal part.

% Find the largest on and off tridiagonal elements. We use abs twice to
% allow for the case where B is complex (this occurs when A is a
% complexified quaternion matrix).

T1 = max(max(abs(abs(O)))); % Find the largest off tridiagonal element.
T2 =     max(abs(abs(D)));  % Find the largest     tridiagonal element.

% NB T2 and/or T1 could be exactly zero (example, if A was zero, or an
% identity matrix). Therefore we do not divide one by the other, but
% instead multiply by a tolerance.

tolerance = eps * 1.0e3; % This is empirically determined.

if T1 > T2 * tolerance
    warning('QTFM:inaccuracy', ...
            'Result of tridiagonalization was not accurately tridiagonal.');
    disp('Information: largest on- and off-tridiagonal moduli were:');
    disp(T2);
    disp(T1);
end

% Verify that the diagonal elements have neglible vector parts.

T1 = max(abs(abs(s(D)))); % The largest scalar modulus of the tridiagonal
                          % result.
T2 = max(abs(abs(v(D)))); % The largest vector modulus of the tridiagonal
                          % result.

if T2 > T1 * tolerance
    warning('QTFM:inaccuracy', ...
            'Result of tridiagonalization was not accurately scalar.')
    disp('Information: largest on-tridiagonal vector and scalar moduli were:');
    disp(T2);
    disp(T1);
end

R = s(B - O); % Subtract the off-tridiagonal part and take the scalar part
              % of the result.

% $Id: tridiagonalize.m 1004 2017-11-15 17:14:09Z sangwine $

