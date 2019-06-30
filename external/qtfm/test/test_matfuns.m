function test_matfuns
% Test code for various Matlab functions that work with quaternions and
% octonions.

% Copyright (c) 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing Matlab functions that work with quaternions ...')

% Check circshift. Shift the data and shift it back again.

A = quaternion(randi(100, 20, 10));
check(all(all(circshift(circshift(A, [3,4]), [-3,-4]) == A)), ...
             'circshift fails test');

% cov - Covariance. This test compares complex data with an isomorphic
% quaternion array. The results should be numerically identical to within
% rounding error.

T = 1e-9;

B = complex(randn(10,4), randn(10,4));
Q = real(B) + imag(B) .* qi; % Make a quaternion with the same data.
CB = cov(B);
CQ = cov(Q);

compare(real(CB), s(CQ), T, 'cov fails test in real part');
compare(imag(CB), x(CQ), T, 'cov fails test in imag part');

% dot - Vector dot product (using the quaternion conjugate)

A = randq(10, 1);
B = randq(10, 1);

compare(dot(A, B), A' * B, T, 'dot (vector dot product) fails test')

% TODO Replace flipdim with flip (new in R2013a, so we would need to make
% sure the test code insists on this release or later(

% flipdim - Flip matrix along specified dimension
% fliplr - Flip matrix in left/right direction
% flipud - Flip matrix in up/down direction

A = randi(100, 5, 4);

check(~any(any(flipdim(A, 1) - s(flipdim(quaternion(A), 1)))), ...
                                'flipdim fails test 1')
check(~any(any(flipdim(A, 2) - s(flipdim(quaternion(A), 2)))), ...
                                'flipdim fails test 2')

check(~any(any(fliplr(A) - s(fliplr(quaternion(A))))), 'fliplr fails test')
check(~any(any(flipud(A) - s(flipud(quaternion(A))))), 'flipud fails test')

% isequal  - True if arrays are numerically equal
% isscalar - True if array is a scalar
% iscolumn - True if array is a column vector
% ismatrix - True if array is a matrix
% isrow    - True if array is a row vector
% isvector - True if array is a vector

A = randq(10, 1);

check( isequal(A, A),              'isequal fails test 1')
check(~isequal(A, randq(size(A))), 'isequal fails test 2')

check( isscalar(randq), 'isscalar fails test 1')
check(~isscalar(A),     'isscalar fails test 2')

check( iscolumn(A),   'iscolumn fails test 1')
check(~iscolumn(A.'), 'iscolumn fails test 2')

check( ismatrix(randq(2,2)),   'ismatrix fails test 1')
check(~ismatrix(randq(2,2,2)), 'ismatrix fails test 2')

check( isrow(A.'), 'isrow fails test 1')
check(~isrow(A),   'isrow fails test 2')

check(isvector(A),           'isvector fails test 1')
check(~isvector(randq(2,2)), 'isvector fails test 2')

% kron - Kronecker product

A = randi(100, 4, 5);
B = randi(100, 2, 3);

check(~any(any(kron(A, B) - s(kron(quaternion(A), quaternion(B))))), ...
       'Kronecker product fails test')

% nchoosek - Combinations of a vector
% perms    - Permutations

A = randi(100, [4, 1]);

check(all(all(nchoosek(A, 2) == s(nchoosek(quaternion(A), 2)))), ...
                                                  'nchoosek fails test')
check(all(all(perms(A) == s(perms(quaternion(A))))), 'perms fails test')

% pinv - Pseudoinverse TODO pinv fails with error, needs a fix in QFTM

A = randq(4,5);
B = pinv(A);

compare(A * B, eyeq(4), T, 'pinv fails test')

% rank - Matrix rank

A = randq(4, 1);
B = randq(1, 4);
C = A * B; assert(all(size(C) == [4, 4]));

check(rank(C) == 1,          'rank fails test 1')
check(rank(randq(3,3)) == 3, 'rank fails test 2')

% rot90 - Rotate matrix 90 degrees

A = randq(4);

check(all(all(rot90(rot90(rot90(rot90(A)))) == A)), 'rot90 fails test')

% std - Standard deviation
% var - Variance

A = randq(10,1);

compare(std(A), sqrt(sum(abs(A - mean(A)).^2) ./ (9)), T, 'std fails test')
compare(var(A),      sum(abs(A - mean(A)).^2) ./ (9),  T, 'var fails test')

% trace - Sum of diagonal elements

A = randq(4, 4);

check(trace(A) == sum(diag(A)), 'trace fails test')

disp('Passed');

% -----------------------------------------------------------

disp('Testing Matlab functions that work with octonions ...')

% Check circshift. Shift the data and shift it back again.

A = octonion(randi(100, 20, 10));
check(all(all(circshift(circshift(A, [3,4]), [-3,-4]) == A)), ...
             'circshift fails test');

% cov - Covariance. This test compares complex data with an isomorphic
% octonion array. The results should be numerically identical to within
% rounding error.

T = 1e-9;

B = complex(randn(10,4), randn(10,4));
Q = real(B) + imag(B) .* oi; % Make an octonion with the same data.
CB = cov(B);
CQ = cov(Q);

compare(real(CB), part(CQ, 1), T, 'cov fails test in real part');
compare(imag(CB), part(CQ, 2), T, 'cov fails test in imag part');

% dot - Vector dot product (using the quaternion conjugate)

A = rando(10, 1);
B = rando(10, 1);

compare(dot(A, B), A' * B, T, 'dot (vector dot product) fails test')

% TODO Replace flipdim with flip (new in R2013a, so we would need to make
% sure the test code insists on this release or later(

% flipdim - Flip matrix along specified dimension
% fliplr - Flip matrix in left/right direction
% flipud - Flip matrix in up/down direction

A = randi(100, 5, 4);

check(~any(any(flipdim(A, 1) - part(flipdim(octonion(A), 1), 1))), ...
                                'flipdim fails test 1')
check(~any(any(flipdim(A, 2) - part(flipdim(octonion(A), 2), 1))), ...
                                'flipdim fails test 2')

check(~any(any(fliplr(A) - part(fliplr(octonion(A)), 1))), 'fliplr fails test')
check(~any(any(flipud(A) - part(flipud(octonion(A)), 1))), 'flipud fails test')

% isequal  - True if arrays are numerically equal
% isscalar - True if array is a scalar
% iscolumn - True if array is a column vector
% ismatrix - True if array is a matrix
% isrow    - True if array is a row vector
% isvector - True if array is a vector

A = rando(10, 1);

check( isequal(A, A),              'isequal fails test 1')
check(~isequal(A, rando(size(A))), 'isequal fails test 2')

check( isscalar(rando), 'isscalar fails test 1')
check(~isscalar(A),     'isscalar fails test 2')

check( iscolumn(A),   'iscolumn fails test 1')
check(~iscolumn(A.'), 'iscolumn fails test 2')

check( ismatrix(rando(2,2)),   'ismatrix fails test 1')
check(~ismatrix(rando(2,2,2)), 'ismatrix fails test 2')

check( isrow(A.'), 'isrow fails test 1')
check(~isrow(A),   'isrow fails test 2')

check(isvector(A),           'isvector fails test 1')
check(~isvector(rando(2,2)), 'isvector fails test 2')

% kron - Kronecker product

A = randi(100, 4, 5);
B = randi(100, 2, 3);

check(~any(any(kron(A, B) - part(kron(octonion(A), octonion(B)), 1))), ...
       'Kronecker product fails test')

% nchoosek - Combinations of a vector
% perms    - Permutations

A = randi(100, [4, 1]);

check(all(all(nchoosek(A, 2) == part(nchoosek(octonion(A), 2), 1))), ...
                                                      'nchoosek fails test')
check(all(all(perms(A) == part(perms(octonion(A)), 1))), 'perms fails test')

% pinv - Pseudoinverse - cannot work without an octonion SVD.

% A = randq(4,5);
% B = pinv(A);
% 
% compare(A * B, eyeq(4), T, 'pinv fails test')

% rank - Matrix rank - cannot work without an octonion SVD.

% A = rando(4, 1);
% B = rando(1, 4);
% C = A * B; assert(all(size(C) == [4, 4]));
% 
% check(rank(C) == 1,          'rank fails test 1')
% check(rank(rando(3,3)) == 3, 'rank fails test 2')

% rot90 - Rotate matrix 90 degrees

A = rando(4);

check(all(all(rot90(rot90(rot90(rot90(A)))) == A)), 'rot90 fails test')

% std - Standard deviation
% var - Variance

A = rando(10,1);

compare(std(A), sqrt(sum(abs(A - mean(A)).^2) ./ (9)), T, 'std fails test')
compare(var(A),      sum(abs(A - mean(A)).^2) ./ (9),  T, 'var fails test')

% trace - Sum of diagonal elements

A = rando(4, 4);

check(trace(A) == sum(diag(A)), 'trace fails test')

disp('Passed')

end

% $Id: test_matfuns.m 1004 2017-11-15 17:14:09Z sangwine $
