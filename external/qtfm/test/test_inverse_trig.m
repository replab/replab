function test_inverse_trig
% Test code for the inverse quaternion trigonometric functions.

% Copyright (c) 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing inverse trigonometric functions ...');

T = 1e-14;

% Test 1. Real quaternion data.

q = quaternion(randn(100,100), randn(100,100), randn(100,100), randn(100,100));

compare(sin(asin(q)), q, T,...
    'quaternion/asin failed test 1.');
compare(cos(acos(q)), q, T,...
    'quaternion/acos failed test 1.');
compare(tan(atan(q)), q, T,...
    'quaternion/atan failed test 1.');

% Test 2. Complex quaternion data.

% The functions do not support this at present so there is no test.

disp('Passed');

% $Id: test_inverse_trig.m 1004 2017-11-15 17:14:09Z sangwine $

