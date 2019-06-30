function test_sqrt
% Test code for the quaternion sqrt function.

% Copyright (c) 2006, 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing square root function ...')

T = 1e-12;

% Test 1. Real quaternion data.

q = randq(100);

compare(q, sqrt(q) .^2, T, 'quaternion/sqrt failed test 1.');

% Test 2. Complex quaternion data.

b = complex(q, randq(100));

compare(b, sqrt(b) .^2, T, 'quaternion/sqrt failed test 2.');

% Test 3. Complex quaternion data with undefined axis.

b = quaternion(scalar(b));

compare(b, sqrt(b) .^2, T, 'quaternion/sqrt failed test 3.');

% Test 4. Complex quaternion data with mixture of defined and undefined
% axes.

b = complex(q, randq(100));
n = rand(100) > 0.5; % Random logical array, 50:50 true:false.
b(n) = quaternion(scalar(b(n))); % Null out the vector parts.

compare(b, sqrt(b) .^2, T, 'quaternion/sqrt failed test 4.');

% Test 5. Quaternions with real scalar parts and null vector parts.
% Some of the values are negative and should have imaginary roots.

b = quaternion(randn(100));
compare(b, sqrt(b).^2, T, 'quaternion/sqrt failed test 5.');

% Test 6. Real octonion data.

q = rando(100) .* randn(100);

compare(q, sqrt(q) .^2, T, 'octonion/sqrt failed test 6.');

% % Test 7. Complex octonion data.
% 
% b = complex(q, rando(100) .* randn(100));
% 
% compare(b, sqrt(b) .^2, T, 'octonion/sqrt failed test 7.');
% 
% % Test 8. Complex octonion data with undefined axis.
% 
% b = octonion(scalar(b));
% 
% compare(b, sqrt(b) .^2, T, 'octonion/sqrt failed test 8.');

% % Test 9. Complex octonion data with mixture of defined and undefined
% % axes.
% 
% b = complex(q, rando(100));
% n = rand(100) > 0.5; % Random logical array, 50:50 true:false.
% b(n) = octonion(scalar(b(n))); % Null out the vector parts.
% 
% compare(b, sqrt(b) .^2, T, 'octonion/sqrt failed test 9.');
warning('Some complex octonion square root tests are disabled.')

% Test 10. Octonions with real scalar parts and null vector parts.
% Some of the values are negative and should have imaginary roots.

b = octonion(randn(100));
compare(b, sqrt(b).^2, T, 'octonion/sqrt failed test 10.');
          
disp('Passed');

% $Id: test_sqrt.m 1004 2017-11-15 17:14:09Z sangwine $
