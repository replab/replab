function test_factorq
% Test code for the quaternion factorq function.

% Copyright (c) 2016 Stephen J. Sangwine and Todd A. Ell
% (Contributed code.)
% See the file : Copyright.m for further details.

disp('Testing factorq function ...');

T = 1e-13;

mu = randv;

q = randq(10);

[f, g] = factorq(q, mu, 'R');

compare(q, f .* g, T, 'quaternion/factorq failed test 1.');

[f, g] = factorq(q, mu, 'L');

compare(q, g .* f, T, 'quaternion/factorq failed test 2.');

[f, ~, theta] = factorq(q, mu, 'R');

compare(q, f .* exp(mu .* theta), T, 'quaternion/factorq failed test 3.');

[f, ~, theta] = factorq(q, mu, 'L');

compare(q, exp(mu .* theta) .* f, T, 'quaternion/factorq failed test 4.');

disp('Passed');

end

% $Id: test_factorq.m 1004 2017-11-15 17:14:09Z sangwine $
