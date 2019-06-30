function test_slerp
% Test code for the slerp and spherical mean functions.

% This also tests the exponential and log functions, plus sqrt/conj etc
% because it uses the power function.

% Copyright (c) 2008, 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing slerp, spherical mean and randf/randvmf functions ...');

T = 1e-12;

compare(slerp(qi, qj, 0.5), ...
        unit(qi+qj),           T, 'quaternion/slerp failed test 1.');
compare(slerp(qi, qj, 0:0.5:1), ...
        [qi, unit(qi+qj), qj], T, 'quaternion/slerp failed test 2.');
    
% Test with random data, using slerp to interpolate half-way between the
% element values.
    
p = unit(quaternion(randn(2), randn(2), randn(2), randn(2)));
q = unit(quaternion(randn(2), randn(2), randn(2), randn(2)));

compare(slerp(p, q, 0.5), unit(p+q), T, 'quaternion/slerp failed test 3.');

% TODO Add a more complex test that does not depend on interpolating half
% way between given values. The problem is how to compute the results
% without using slerp, in order to check the results?

% Now test spherical mean and randf/randvmf functions ...

% The first two tests compare the spherical mean against slerp (the former
% is a generalisation of the latter). We can only do this with a pair of
% quaternions, but we test both full and pure cases.

p = randq;
q = randq;

compare(slerp(p, q, 0.5), spherical_mean([p, q]), T, ...
                               'quaternion/spherical_mean failed test 1.');
p = randv;
q = randv;
compare(slerp(p, q, 0.5), spherical_mean([p, q]), T, ...
                               'quaternion/spherical_mean failed test 2.');

% To do a more sophisticated test is more tricky. It makes sense to use the
% randf and randvmf functions, to see that spherical_mean gives a result
% close to the mean direction parameter. As often, the choice of the
% tolerance parameter here is rough and ready.

mu = randq;
p = randvmf(mu, 10, 1, 1000);
nu = spherical_mean(p);
compare(mu, nu, 0.05, 'quaternion/spherical_mean failed test 3.');

mu = randv;
p = randf(mu, 10, 1, 1000);
nu = spherical_mean(p);
compare(mu, nu, 0.05, 'quaternion/spherical_mean failed test 4.');

disp('Passed');

% TODO Add tests for the octonion slerp function when it works (requires
% log to be implemented).

% $Id: test_slerp.m 1004 2017-11-15 17:14:09Z sangwine $
