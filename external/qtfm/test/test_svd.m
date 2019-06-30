function test_svd
% Test code for the quaternion svd function.

% Copyright (c) 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing singular value decompositions ...')

T = 1e-10;

% This is quite a long function because of the need to test the
% function with square, wide and tall arrays, and with 1, 2 or
% 3 output parameters, and with normal and economy decompositions.

% Square array.

disp('  svd on real quaternions ...')
disp('  svd on square array ...')

A = quaternion(randn(5,5),randn(5,5),randn(5,5),randn(5,5));

% Normal decomposition.

[U, B, V] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 1A')

[U, B] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 1B')

S = svd(A);
compare(S, diag(B), T, 'quaternion/svd failed test 1C')

% Economy (0) decomposition.

[U, B, V] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 2A')

[U, B] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 2B')

S = svd(A,0);
compare(S, diag(B), T, 'quaternion/svd failed test 2C')

% Economy (''econ'') decomposition.

[U, B, V] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 3A')

[U, B] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 3B')

S = svd(A,'econ');
compare(S, diag(B), T, 'quaternion/svd failed test 3C')

% Wide array.

disp('  svd on wide   array ...')

A = quaternion(randn(5,7),randn(5,7),randn(5,7),randn(5,7));

% Normal decomposition.

[U, B, V] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 4A')

[U, B] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 4B')

S = svd(A);
compare(S, diag(B), T, 'quaternion/svd failed test 4C')

% Economy (0) decomposition.

[U, B, V] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 5A')

[U, B] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 5B')

S = svd(A,0);
compare(S, diag(B), T, 'quaternion/svd failed test 5C')

% Economy (''econ'') decomposition.

[U, B, V] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 6A')

[U, B] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 6B')

S = svd(A,'econ');
compare(S, diag(B), T, 'quaternion/svd failed test 6C')

% Tall array.

disp('  svd on tall   array ...')

A = A.';

% Normal decomposition.

[U, B, V] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 7A')

[U, B] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 7B')

S = svd(A);
compare(S, diag(B), T, 'quaternion/svd failed test 7C')

% Economy (0) decomposition.

[U, B, V] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 8A')

[U, B] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 8B')

S = svd(A,0);
compare(S, diag(B), T, 'quaternion/svd failed test 8C')

% Economy (''econ'') decomposition.

[U, B, V] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 9A')

[U, B] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 9B')

S = svd(A,'econ');
compare(S, diag(B), T, 'quaternion/svd failed test 9C')

% Now repeat the whole lot with complex data ........

% Temporarily removed due to incorrect results from SVD on complex
% quaternions (see the svd.m function for a TODO note).

% TODO
warning('SVD is not tested on complex quaternion matrices.')

disp('Passed')
return

T = 1e-8; % The error margin for complex data is much more
          % relaxed than for real data, because the results
          % are not so accurate.

% Square array.

disp('  svd on complex quaternions ...')
disp('  svd on square array ...')

A = quaternion(complex(randn(5,5),randn(5,5)),...
               complex(randn(5,5),randn(5,5)),...
               complex(randn(5,5),randn(5,5)),...
               complex(randn(5,5),randn(5,5)));

% Normal decomposition.

[U, B, V] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 10A')

[U, B] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 10B')

S = svd(A);
compare(S, diag(B), T, 'quaternion/svd failed test 10C')

% Economy (0) decomposition.

[U, B, V] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 11A')

[U, B] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 11B')

S = svd(A,0);
compare(S, diag(B), T, 'quaternion/svd failed test 11C')

% Economy (''econ'') decomposition.

[U, B, V] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 12A')

[U, B] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 12B')

S = svd(A,'econ');
compare(S, diag(B), T, 'quaternion/svd failed test 12C')

% Wide array.

disp('  svd on wide   array ...')

A = quaternion(complex(randn(5,7),randn(5,7)),...
               complex(randn(5,7),randn(5,7)),...
               complex(randn(5,7),randn(5,7)),...
               complex(randn(5,7),randn(5,7)));

% Normal decomposition.

[U, B, V] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 13A')

[U, B] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 13B')

S = svd(A);
compare(S, diag(B), T, 'quaternion/svd failed test 13C')

% Economy (0) decomposition.

[U, B, V] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 14A')

[U, B] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 14B')

S = svd(A,0);
compare(S, diag(B), T, 'quaternion/svd failed test 14C')

% Economy (''econ'') decomposition.

[U, B, V] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 15A')

[U, B] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 15B')

S = svd(A,'econ');
compare(S, diag(B), T, 'quaternion/svd failed test 15C')

% Tall array.

disp('  svd on tall   array ...')

A = A.';

% Normal decomposition.

[U, B, V] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 16A')

[U, B] = svd(A);
compare(U*B*V', A, T,  'quaternion/svd failed test 16B')

S = svd(A);
compare(S, diag(B), T, 'quaternion/svd failed test 16C')

% Economy (0) decomposition.

[U, B, V] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 17A')

[U, B] = svd(A,0);
compare(U*B*V', A, T,  'quaternion/svd failed test 17B')

S = svd(A,0);
compare(S, diag(B), T, 'quaternion/svd failed test 17C')

% Economy (''econ'') decomposition.

[U, B, V] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 18A')

[U, B] = svd(A,'econ');
compare(U*B*V', A, T,  'quaternion/svd failed test 18B')

S = svd(A,'econ');
compare(S, diag(B), T, 'quaternion/svd failed test 18C')

disp('Passed');

% $Id: test_svd.m 1004 2017-11-15 17:14:09Z sangwine $
