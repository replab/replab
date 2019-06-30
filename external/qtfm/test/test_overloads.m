function test_overloads
% Test code for various quaternion overloadings of Matlab functions.

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing quaternion overloadings of various Matlab functions ...')

T = 1e-13;

% Test data, real and complex.

q = quaternion(randn(100), randn(100), randn(100), randn(100));
b = quaternion(complex(randn(100)),complex(randn(100)),...
               complex(randn(100)),complex(randn(100)));

test_f (@diag,    q, T, 'quaternion/diag failed test 1.');
test_f2(@diag, 1, q, T, 'quaternion/diag failed test 2.');
test_f (@diag,    b, T, 'quaternion/diag failed test 3.');
test_f2(@diag, 1, b, T, 'quaternion/diag failed test 4.');

test_f (@tril,    q, T, 'quaternion/tril failed test 1.');
test_f2(@tril, 1, q, T, 'quaternion/tril failed test 2.');
test_f (@tril,    b, T, 'quaternion/tril failed test 3.');
test_f2(@tril, 1, b, T, 'quaternion/tril failed test 4.');

test_f (@triu,    q, T, 'quaternion/triu failed test 1.');
test_f2(@triu, 1, q, T, 'quaternion/triu failed test 2.');
test_f (@triu,    b, T, 'quaternion/triu failed test 3.');
test_f2(@triu, 1, b, T, 'quaternion/triu failed test 4.');

test_f (@mean,    q, T, 'quaternion/mean failed test 1.');
test_f2(@mean, 2, q, T, 'quaternion/mean failed test 2.');
test_f (@mean,    b, T, 'quaternion/mean failed test 3.');
test_f2(@mean, 2, b, T, 'quaternion/mean failed test 4.');

test_f (@sum,     q, T, 'quaternion/sum failed test 1.');
test_f2(@sum,  2, q, T, 'quaternion/sum failed test 2.');
test_f (@sum,     b, T, 'quaternion/sum failed test 3.');
test_f2(@sum,  2, b, T, 'quaternion/sum failed test 4.');

test_f (@ceil,    q, T, 'quaternion/ceil failed test 1.');
test_f (@ceil,    b, T, 'quaternion/ceil failed test 2.');

test_f (@floor,   q, T, 'quaternion/floor failed test 1.');
test_f (@floor,   b, T, 'quaternion/floor failed test 2.');

test_f (@fix,     q, T, 'quaternion/fix failed test 1.');
test_f (@fix,     b, T, 'quaternion/fix failed test 2.');

test_f (@round,   q, T, 'quaternion/round failed test 1.');
test_f (@round,   b, T, 'quaternion/round failed test 2.');

disp('Passed');

function test_f(f, q, T, M)
% Tests function with handle f using quaternion data q, tolerance T and
% error message M.

compare(f(q),    quaternion(f(s(q)), f(x(q)), f(y(q)), f(z(q))), T, M);
compare(f(v(q)), quaternion(         f(x(q)), f(y(q)), f(z(q))), T, M);

function test_f2(f, k, q, T, M)
% Tests function with handle f using quaternion data q, second parameter k,
% tolerance T and error message M.

compare(f(q, k),    ...
    quaternion(f(s(q), k), f(x(q), k), f(y(q), k), f(z(q), k)), T, M);
compare(f(v(q), k), ...
    quaternion(            f(x(q), k), f(y(q), k), f(z(q), k)), T, M);

% $Id: test_overloads.m 1004 2017-11-15 17:14:09Z sangwine $
