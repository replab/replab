function test_quaternion_identities
% Test code for quaternion identities. Modelled on the octonion case, but
% there are fewer identities to test.

% Copyright (c) 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing quaternion identities ...');

x = randq(10); y = randq(10); z = randq(10);

% Check the commutator function.

compare(commutator(x, y), ...
       -commutator(y, x), 1e-12, 'Quaternion commutator fails test 1');
compare(commutator(x, y, 'diff'), ...
       -commutator(y, x, 'diff'), ...
                          1e-12, 'Quaternion commutator fails test 2');
compare(x .* y .* commutator(x, y, 'prod'), y .* x, ...
                          1e-12, 'Quaternion commutator fails test 3.');

                      
% Hall identity.

c2 = commutator(x, y).^2;

compare(c2 .* z, z .* c2, 1e-12, 'Quaternion Hall identity fails.')

disp('Passed');

end

% $Id: test_quaternion_identities.m 1004 2017-11-15 17:14:09Z sangwine $
