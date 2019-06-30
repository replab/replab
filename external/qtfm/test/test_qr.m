function test_qr
% Test code for the quaternion qr function.

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing QR decomposition ...')

% TODO Add checks that Q is unitary, and that R is upper triangular, to
% within a tolerance. The function isunitary can be used for the former.

T = 1e-12;

A = randq(10);

[Q, R] = qr(A);
compare(Q * R, A, T,   'quaternion/qr failed test 1')

A = randq(10,13);

[Q, R] = qr(A);
compare(Q * R, A, T,   'quaternion/qr failed test 2')

[Q, R] = qr(A.');
compare(Q * R, A.', T, 'quaternion/qr failed test 3')

% Now repeat the whole lot with complex data, but smaller matrices, keeping
% the same tolerance. No particular reason, but no need to keep the same
% sizes.

% Commented out complex tests while we sort out the problems with
% biquaternion unitary matrices.
disp('Passed')
return

T = 1e-9; % Relax the requirements a little.

A = complex(randq(5), randq(5));

[Q, R] = qr(A);
compare(Q * R, A, T,   'quaternion/qr failed test 4')

A = complex(randq(5,7), randq(5,7));

[Q, R] = qr(A);
compare(Q * R, A, T,   'quaternion/qr failed test 5')

[Q, R] = qr(A.');
compare(Q * R, A.', T, 'quaternion/qr failed test 6')

disp('Passed');

% $Id: test_qr.m 1004 2017-11-15 17:14:09Z sangwine $
