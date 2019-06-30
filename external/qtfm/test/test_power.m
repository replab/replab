function test_power
% Test code for the quaternion and octonion power functions (.^).

% This also tests the exponential and log functions, plus sqrt/conj etc.

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing power function (.^) ...');

T = 1e-10;

q = quaternion(randn(100), randn(100), randn(100), randn(100));

% Scalar powers handled as special cases.

compare(q,       q.^1,   T, 'quaternion/power failed test 1.');
compare(q.*q,    q.^2,   T, 'quaternion/power failed test 2.');
compare(sqrt(q), q.^0.5, T, 'quaternion/power failed test 3.');
compare(q,  (q.^-1).^-1, T, 'quaternion/power failed test 4.');
compare(q,(q.^-0.5).^-2, T, 'quaternion/power failed test 5.');

% General scalar power. 3 is not handled as a special case, so we can
% compare it with cubing explicitly.

compare(q.*q.*q, q.^3, T,   'quaternion/power failed test 6.');

% Scalar raised to a vector power.

compare(qi.^[0 1 2 3], [quaternion( 1,  0, 0, 0), ...
                        quaternion( 0,  1, 0, 0), ...
                        quaternion(-1,  0, 0, 0), ...
                        quaternion( 0, -1, 0, 0)],...
                        T, 'quaternion/power failed test 7.');
                    
% Test the octonion power function.
                    
o = rando(100);

% Scalar powers handled as special cases.

compare(o,       o.^1,   T, 'octonion/power failed test 1.');
compare(o.*o,    o.^2,   T, 'octonion/power failed test 2.');
compare(sqrt(o), o.^0.5, T, 'octonion/power failed test 3.');
compare(o,  (o.^-1).^-1, T, 'octonion/power failed test 4.');
compare(o,(o.^-0.5).^-2, T, 'octonion/power failed test 5.');

p = o .* o;
for j = 3:9
    p = p .* o;
    compare(o .^ j, p, T, ['octonion/power failed test 6 with j = ' num2str(j)]);
    compare(o .^ -j, p.^-1, T, ['octonion/power failed test 7 with j = ' num2str(j)]);    
end

disp('Passed');

% $Id: test_power.m 1004 2017-11-15 17:14:09Z sangwine $
