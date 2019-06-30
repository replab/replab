function test_sum_diff
% Test code for the quaternion sum and difference and related functions.

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing sum and diff functions and similar ...');

T = 1e-10;

q = quaternion(randn(10), randn(10), randn(10), randn(10));

% Sum function.

t1 = sum(q);

compare(s(t1), sum(s(q)), T, 'quaternion/sum failed test 1.');
compare(x(t1), sum(x(q)), T, 'quaternion/sum failed test 2.');
compare(y(t1), sum(y(q)), T, 'quaternion/sum failed test 3.');
compare(z(t1), sum(z(q)), T, 'quaternion/sum failed test 4.');

clear t1

t2 = sum(q, 2);

compare(s(t2), sum(s(q), 2), T, 'quaternion/sum failed test 5.');
compare(x(t2), sum(x(q), 2), T, 'quaternion/sum failed test 6.');
compare(y(t2), sum(y(q), 2), T, 'quaternion/sum failed test 7.');
compare(z(t2), sum(z(q), 2), T, 'quaternion/sum failed test 8.');

clear t2

% Cumulative sum function.

c1 = cumsum(q);

compare(s(c1), cumsum(s(q)), T, 'quaternion/cumsum failed test 1.');
compare(x(c1), cumsum(x(q)), T, 'quaternion/cumsum failed test 2.');
compare(y(c1), cumsum(y(q)), T, 'quaternion/cumsum failed test 3.');
compare(z(c1), cumsum(z(q)), T, 'quaternion/cumsum failed test 4.');

clear c1

c2 = cumsum(q, 2);

compare(s(c2), cumsum(s(q), 2), T, 'quaternion/cumsum failed test 5.');
compare(x(c2), cumsum(x(q), 2), T, 'quaternion/cumsum failed test 6.');
compare(y(c2), cumsum(y(q), 2), T, 'quaternion/cumsum failed test 7.');
compare(z(c2), cumsum(z(q), 2), T, 'quaternion/cumsum failed test 8.');

% Diff function.

d1 = diff(q);

compare(s(d1), diff(s(q)), T, 'quaternion/diff failed test 1.');
compare(x(d1), diff(x(q)), T, 'quaternion/diff failed test 2.');
compare(y(d1), diff(y(q)), T, 'quaternion/diff failed test 3.');
compare(z(d1), diff(z(q)), T, 'quaternion/diff failed test 4.');

clear d1

d2 = diff(q, 2);

compare(s(d2), diff(s(q), 2), T, 'quaternion/diff failed test 5.');
compare(x(d2), diff(x(q), 2), T, 'quaternion/diff failed test 6.');
compare(y(d2), diff(y(q), 2), T, 'quaternion/diff failed test 7.');
compare(z(d2), diff(z(q), 2), T, 'quaternion/diff failed test 8.');

clear d2

d3 = diff(q, 2, 2);

compare(s(d3), diff(s(q), 2, 2), T, 'quaternion/diff failed test 9.');
compare(x(d3), diff(x(q), 2, 2), T, 'quaternion/diff failed test 10.');
compare(y(d3), diff(y(q), 2, 2), T, 'quaternion/diff failed test 11.');
compare(z(d3), diff(z(q), 2, 2), T, 'quaternion/diff failed test 12.');

clear d3

% Prod function. Mostly we test with empty vector parts, because we can
% then compare with the result of prod on the scalar part.

q = quaternion(randn(10,1)); % Vector.

compare(prod(s(q )), s(prod(q )), T, 'quaternion/prod failed test 1.');
compare(prod(s(q')), s(prod(q')), T, 'quaternion/prod failed test 2.');

q = quaternion(randn(10)); % Matrix.

compare(prod(s(q)   ), s(prod(q)   ), T,  'quaternion/prod failed test 3.');
compare(prod(s(q), 1), s(prod(q, 1)), T,  'quaternion/prod failed test 4.');
compare(prod(s(q), 2), s(prod(q, 2)), T,  'quaternion/prod failed test 5.');

check(prod([qi qj qk]) == qi .* qj .* qk, 'quaternion/prod failed test 6.');

% Now make a small matrix for testing on full quaternions, being careful
% about the order of multiplication (columns first).

q = quaternion(randn(2), randn(2), randn(2), randn(2));

compare(prod(prod(q)), q(1,1) .* q(2,1) .* q(1,2) .* q(2,2), T,...
                                          'quaternion/prod failed test 7.');

% Test using normed property of quaternions.
                                      
q = quaternion(randn(4), randn(4), randn(4), randn(4));

compare(abs(prod(prod(q))), prod(prod(abs(q))), T,...
                                          'quaternion/prod failed test 8.');

% Cumulative prod function. Mostly we test with empty vector parts, because
% we can then compare with the result of prod on the scalar part.

q = quaternion(randn(10,1)); % Vector.

compare(cumprod(s(q )), s(cumprod(q )), T, 'quaternion/cumprod failed test 1.');
compare(cumprod(s(q')), s(cumprod(q')), T, 'quaternion/cumprod failed test 2.');

q = quaternion(randn(10)); % Matrix.

compare(cumprod(s(q)   ), s(cumprod(q)   ), T, 'quaternion/cumprod failed test 3.');
compare(cumprod(s(q), 1), s(cumprod(q, 1)), T, 'quaternion/cumprod failed test 4.');
compare(cumprod(s(q), 2), s(cumprod(q, 2)), T, 'quaternion/cumprod failed test 5.');

check(cumprod([qi qj qk]) == [0 + qi, 0 + qi .* qj, qi .* qj .* qk],...
                                               'quaternion/cumprod failed test 6.');

% Test using normed property of quaternions.
                                      
q = quaternion(randn(4), randn(4), randn(4), randn(4));

compare(abs(cumprod(cumprod(q))), cumprod(cumprod(abs(q))), T,...
                                          'quaternion/cumprod failed test 7.');
disp('Passed');

% $Id: test_sum_diff.m 1004 2017-11-15 17:14:09Z sangwine $

