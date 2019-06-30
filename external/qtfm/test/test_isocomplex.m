function test_isocomplex
% Test code for various quaternion functions that use the isocomplex
% method (see @quaternion/private/iscomplex.m).

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing quaternion functions based on isocomplex numbers ...')

T = 1e-10;

% Test data, real and complex. Some of the complex tests are commented out
% because the functions concerned are not yet implemented for complex
% quaternions and will raise an error.

c = complex(randn(10), randn(10)); r = real(c);

test_f (@sin, r, T, 'quaternion/sin failed test 1.');
test_f (@sin, c, T, 'quaternion/sin failed test 2.');

test_f (@cos, r, T, 'quaternion/cos failed test 1.');
test_f (@cos, c, T, 'quaternion/cos failed test 2.');

test_f (@tan, r, T, 'quaternion/tan failed test 1.');
test_f (@tan, c, T, 'quaternion/tan failed test 2.');

test_f (@asin, r, T, 'quaternion/asin failed test 1.');
%test_f (@asin, c, T, 'quaternion/asin failed test 2.');

test_f (@acos, r, T, 'quaternion/acos failed test 1.');
%test_f (@acos, c, T, 'quaternion/acos failed test 2.');

test_f (@atan, r, T, 'quaternion/atan failed test 1.');
%test_f (@atan, c, T, 'quaternion/atan failed test 2.');

test_f (@sinh, r, T, 'quaternion/sinh failed test 1.');
test_f (@sinh, c, T, 'quaternion/sinh failed test 2.');

test_f (@cosh, r, T, 'quaternion/cosh failed test 1.');
test_f (@cosh, c, T, 'quaternion/cosh failed test 2.');

test_f (@tanh, r, T, 'quaternion/tanh failed test 1.');
test_f (@tanh, c, T, 'quaternion/tanh failed test 2.');

test_f (@asinh, r, T, 'quaternion/asinh failed test 1.');
%test_f (@asinh, c, T, 'quaternion/asinh failed test 2.');

test_f (@acosh, r, T, 'quaternion/acosh failed test 1.');
%test_f (@acosh, c, T, 'quaternion/acosh failed test 2.');

test_f (@atanh, r, T, 'quaternion/atanh failed test 1.');
%test_f (@atanh, c, T, 'quaternion/atanh failed test 2.');

test_f (@log, r, T, 'quaternion/log failed test 1.'); % TODO this test fails because of pi ambiguity in the imaginary part.
%test_f (@log, c, T, 'quaternion/log failed test 2.');

%test_f (@sign, r, T, 'quaternion/sign failed test 1.'); % Testing the sign function like this is incorrect TODO.
%test_f (@sign, c, T, 'quaternion/sign failed test 2.');

disp('Passed');

function test_f(f, q, T, M)
% Tests function with handle f using real/complex data q, tolerance T and
% error message M.

compare(quaternion(f(q)), f(quaternion(q)), T, M);
end

end

% $Id: test_isocomplex.m 1004 2017-11-15 17:14:09Z sangwine $
