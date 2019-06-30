function test_mat
% Test code for the quaternion expm, logm and sqrtm functions.

% Copyright (c) 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing expm, logm and sqrtm functions ...')

T = 1e-12;

% TODO Change the tests to use bigger matrices when the functions have been
% properly implemented rather than using adjoint matrices. (Accuracy with
% adjoint matrices is known to be low, especially for larger matrices.)

N = 6;

% Test 1. Real quaternion data.

q = randq(N);

compare(q, sqrtm(q)^2, T, 'quaternion/sqrtm failed test 1.');

compare(q, logm(expm(q)), T, 'quaternion/expm/logm failed test 1.');

% Test 2. Complex quaternion data.

b = complex(q, randq(N));

compare(b, sqrtm(b)^2, T, 'quaternion/sqrtm failed test 2.');

compare(q, logm(expm(q)), T, 'quaternion/expm/logm failed test 2.');

disp('Passed');

% $Id: test_mat.m 1004 2017-11-15 17:14:09Z sangwine $

