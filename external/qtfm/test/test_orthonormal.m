function test_orthonormal
% Test code for the orthogonal and orthonormal_basis functions.

% Copyright (c) 2010, 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing orthogonal and orthonormal_basis functions ...')

T = 1e-12;

% The tests here depend on the choices made in the code for orthogonal.
% They are not unique solutions, but since previous code may have given
% results dependent on these solutions this test code verifies them for
% backwards compatibility.

compare( qj, orthogonal( qi), T, 'orthogonal failed test 1');
compare( qk, orthogonal( qj), T, 'orthogonal failed test 2');
compare( qi, orthogonal( qk), T, 'orthogonal failed test 3');
compare(-qj, orthogonal(-qi), T, 'orthogonal failed test 4');
compare(-qk, orthogonal(-qj), T, 'orthogonal failed test 5');
compare(-qi, orthogonal(-qk), T, 'orthogonal failed test 6');

mu = unit(quaternion(1,1,1));
compare(orthogonal(mu),     unit(quaternion(1, -1, 0)), T, 'orthogonal failed test 7');
compare(orthogonal(mu, qi), unit(quaternion(0, 1, -1)), T, 'orthogonal failed test 8');
compare(orthogonal(mu, qj), unit(quaternion(-1, 0, 1)), T, 'orthogonal failed test 9');
compare(orthogonal(mu, qk), unit(quaternion(1, -1, 0)), T, 'orthogonal failed test 10');

nu = unit(quaternion(1,2,3)); % This is not parallel to any of {qi, qj, qk}.
compare(orthogonal(nu), orthogonal(nu, mu), T, 'orthogonal failed test 11');

% Now test the orthogonal function on non-scalar arguments. These tests
% were added in February 2017 when orthogonal was edited to permit
% non-scalar arguments.

A = randv(3,4);
compare(0, scalar_product(A, orthogonal(A)),     T, 'orthogonal failed test 12.');
compare(0, scalar_product(A, orthogonal(A, mu)), T, 'orthogonal failed test 13.');

% Now test with one of the parameters scalar.

compare(0, scalar_product(A, orthogonal(A, randv(size(A)))), T, 'orthogonal failed test 13.');
compare(0, scalar_product(A, orthogonal(randv(size(A)), A)), T, 'orthogonal failed test 14.');

% Now test the orthonormal basis function.

compare(eye(3), orthonormal_basis(qi), T, 'orthonormal_basis failed test 1.');

% Now we try some random values, real and complex, and test for correct
% orthogonality of the basis.

for A = 1:100
   B = orthonormal_basis(randv);
   compare(B * B.', eye(3), T, 'orthonormal_basis failed test 2.');
   C = orthonormal_basis(unit(complex(randv, randv)));
   compare(C * C.', eye(3), T, 'orthonormal_basis failed test 3.');
end

disp('Passed');

% $Id: test_orthonormal.m 1004 2017-11-15 17:14:09Z sangwine $
