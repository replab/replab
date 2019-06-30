function test_octonion_identities
% Test code for the octonion identities. This code tests the Moufang and
% alternative identities, providing another check on the multiplication
% tables, but also verifying powers and similar.

% Copyright (c) 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing octonion identities ...');

% References:
%
% Richard D. Schafer, Introduction to Non-Associative Algebras,
% Academic Press, 1966
%
% Marshall Hall, Projective Planes, Transactions American Mathematical
% Society, 54(2), pp. 229-277, September 1943.

x = rando(10); y = rando(10); z = rando(10);

% Start by checking the commutator and associator functions.

compare(commutator(x, y), ...
       -commutator(y, x), 1e-12, 'Octonion commutator fails test 1');
compare(commutator(x, y, 'diff'), ...
       -commutator(y, x, 'diff'), ...
                          1e-12, 'Octonion commutator fails test 2');
compare(x .* y .* commutator(x, y, 'prod'), y .* x, ...
                          1e-12, 'Octonion commutator fails test 3.');
        
compare(associator(x, y, z), (x .* y) .* z - x .* (y .* z), 1e-12, ...
                                 'Octonion associator fails test 1.');
compare(associator(x, y, z, 'diff'), ...
                             (x .* y) .* z - x .* (y .* z), 1e-12, ...
                                 'Octonion associator fails test 2.');
compare(((x .* y) .* z) .* associator(x, y, z, 'prod'), x .* (y .* z), ...
                          1e-12, 'Octonion associator fails test 3.');                          
   
% Alternative laws (Schafer pp. 23, 27, 28.)

compare(x .* (x .* y), x.^2 .* y, 1e-12, 'Octonion alternative law 1 fails')
compare((y .* x) .* x, y .* x.^2, 1e-12, 'Octonion alternative law 2 fails')

compare(x .* (y .* x), ...
       (x .* y) .* x, 1e-12, 'Octonion alternative law 3 fails')

% Moufang identities

compare((x .* y) .* (z .* x), ...
         x .* (y .* z) .* x, 1e-12, 'Octonion Moufang identity 1 fails')
compare((x .* y) .* z, (x .* z.^-1) .* (z .* y .* z), ...
                             1e-12, 'Octonion Moufang identity 2 fails')
compare(x .* (y .* z), (x .* y .* x) .* (x.^-1 .* z), ...
                             1e-12, 'Octonion Moufang identity 3 fails')
compare((x .* z .* x) .* y, x .* (z .* (x .* y)), ...
                             1e-12, 'Octonion Moufang identity 4 fails')
compare(y .* (x .* z .* x), ((y .* x) .* z) .* x, ...
                             1e-12, 'Octonion Moufang identity 5 fails')

% Associator identity (Schafer section 2.4, page 13)

a = rando(10);

compare(a .* associator(x, y, z) + associator(a, x, y) .* z, ...
        associator(a .* x, y, z) - associator(a, x .* y, z) + ...
        associator(a, x, y .* z), 1e-12, ...
        'Octonion associator identity fails')
    
% Hall identity.

c2 = commutator(x, y).^2;

compare(c2 .* z, z .* c2, 1e-12, 'Octonion Hall identity fails.')
   
disp('Passed');

end

% $Id: test_octonion_identities.m 1004 2017-11-15 17:14:09Z sangwine $
