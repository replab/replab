function test_householder
% Test code for the Householder vector and matrix functions. These
% functions are tested by the SVD test code, but a direct test cannot do
% any harm. The householder_vector code is tested only indirectly, by calls
% from the householder_matrix function.

% Copyright (c) 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing Householder functions ...')

T = 1e-10; % Comparison threshold, arbitrary value.
N = 4;     % Size of test vectors.

% There are four sets of tests to do with real, complex, quaternion, and
% biquaternion vectors.

% TODO Add a test for octonions? Unsure whether Householder will work for
% octonions, so this will need investigating.

% Each test computes a Householder matrix, h,  for a random vector a and a
% nulling vector v, then applies h to a to give r. The checks are that h is
% orthogonal/unitary, that the norm of r matches the norm of a, and that
% all elements of r are zero apart from the first. These tests are done for
% row and column vectors a and v.

v = [1;0;0;0];

% Real case.

a = randn(N,1);

h = householder_matrix(a, v);

compare(h * h', eye(N), T, 'Householder test 1 failed.');

r = h * a; % Apply the transformation to a.

compare(r(2:end), zeros(N - 1, 1), T, 'Householder test 2 failed.');
compare(norm(r), norm(a), T, 'Householder test 3 failed.');

h = householder_matrix(a.', v.');

compare(h * h', eye(N), T, 'Householder test 4 failed.');

r = a.' * h;

compare(r(2:end), zeros(1, N - 1), T, 'Householder test 5 failed.');
compare(norm(r), norm(a), T, 'Householder test 6 failed.');

% Complex case.

a = complex(a, randn(N,1));

h = householder_matrix(a, v);

compare(h * h', eye(N), T, 'Householder test 7 failed.');

r = h * a; % Apply the transformation to a.

compare(r(2:end), zeros(N - 1, 1), T, 'Householder test 8 failed.');
compare(norm(r), norm(a), T, 'Householder test 9 failed.');

h = householder_matrix(a.', v.');

compare(h * h', eye(N), T, 'Householder test 10 failed.');

r = a.' * h;

compare(r(2:end), zeros(1, N - 1), T, 'Householder test 11 failed.');
compare(norm(r), norm(a), T, 'Householder test 12 failed.');

% Quaternion case.

a = randq(N,1);

h = householder_matrix(a, v);

compare(h * h', eyeq(N), T, 'Householder test 13 failed.');

r = h * a; % Apply the transformation to a.

compare(r(2:end), zerosq(N - 1, 1), T, 'Householder test 14 failed.');
compare(norm(r), norm(a), T, 'Householder test 15 failed.');

h = householder_matrix(a.', v.');

compare(h * h', eyeq(N), T, 'Householder test 16 failed.');

r = a.' * h;

compare(r(2:end), zerosq(1, N - 1), T, 'Householder test 17 failed.');
compare(norm(r), norm(a), T, 'Householder test 18 failed.');

% Biquaternion case.

warning('Householder tests in biquaternion case are disabled.')

% a = complex(a, randq(N,1));
% 
% h = householder_matrix(a, v);
% 
% compare(h * h', eyeq(N), T, 'Householder test 19 failed.');
% 
% r = h * a; % Apply the transformation to a.
% 
% compare(r(2:end), zerosq(N - 1, 1), T, 'Householder test 20 failed.');
% compare(norm(r), norm(a), T, 'Householder test 21 failed.');
% 
% h = householder_matrix(a.', v.');
% 
% compare(h * h', eyeq(N), T, 'Householder test 22 failed.');
% 
% r = a.' * h;
% 
% compare(r(2:end), zerosq(1, N - 1), T, 'Householder test 23 failed.');
% compare(norm(r), norm(a), T, 'Householder test 24 failed.');

disp('Passed');

end

% $Id: test_householder.m 1004 2017-11-15 17:14:09Z sangwine $
