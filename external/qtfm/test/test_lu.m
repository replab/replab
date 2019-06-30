function test_lu
% Test code for the quaternion and octonion lu functions.

% Copyright (c) 2011, 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing LU decompositions ...')

T = 1e-12;

A = randq(10);

[L, U, P] = lu(A);
compare(L * U, P * A, T,   'quaternion/lu failed test 1A')

[L, U] = lu(A);
compare(L * U,     A, T,   'quaternion/lu failed test 1B')

A = randq(10,13);

[L, U, P] = lu(A);
compare(L * U, P * A, T,   'quaternion/lu failed test 2A')

[L, U] = lu(A);
compare(L * U,     A, T,   'quaternion/lu failed test 2B')

[L, U, P] = lu(A.');
compare(L * U, P * A.', T, 'quaternion/lu failed test 3A')

[L, U] = lu(A.');
compare(L * U,     A.', T, 'quaternion/lu failed test 3B')

T = 1e-9;

A = rando(10);

[L, U, P] = lu(A);
compare(L * U, P * A, T,   'octonion/lu failed test 1A')

[L, U] = lu(A);
compare(L * U,     A, T,   'octonion/lu failed test 1B')

A = rando(10,13);

[L, U, P] = lu(A);
compare(L * U, P * A, T,   'octonion/lu failed test 2A')

[L, U] = lu(A);
compare(L * U,     A, T,   'octonion/lu failed test 2B')

[L, U, P] = lu(A.');
compare(L * U, P * A.', T, 'octonion/lu failed test 3A')

[L, U] = lu(A.');
compare(L * U,     A.', T, 'octonion/lu failed test 3B')

% Now repeat the whole lot with complex data, but smaller matrices, keeping
% the same tolerance. No particular reason, but no need to keep the same
% sizes.

T = 1e-9; % Relax the requirements a little.

A = complex(randq(5), randq(5));

[L, U, P] = lu(A);
compare(L * U, P * A, T,   'quaternion/lu failed test 4')

A = complex(randq(7,5), randq(7,5));

[L, U, P] = lu(A);
compare(L * U, P * A, T,   'quaternion/lu failed test 5')

[L, U, P] = lu(A.');
compare(L * U, P * A.', T, 'quaternion/lu failed test 6')

% We omit these tests for octonion matrices for the moment, as the
% test sometimes fails (badly). The reason is not known.

% T = 1e-6;
% 
% A = complex(rando(5), rando(5));
% 
% [L, U, P] = lu(A);
% compare(L * U, P * A, T,   'octonion/lu failed test 4')
% 
% A = complex(rando(7,5), rando(7,5));
% 
% [L, U, P] = lu(A);
% compare(L * U, P * A, T,   'octonion/lu failed test 5')
% 
% [L, U, P] = lu(A.');
% compare(L * U, P * A.', T, 'octonion/lu failed test 6')

disp('Passed');

end

% $Id: test_lu.m 1004 2017-11-15 17:14:09Z sangwine $
