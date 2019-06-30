function test_conv
% Test code for the quaternion convolution functions.

% Copyright (c) 2006, 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing convolutions 1D and 2D ...')

T = 1e-10;

% Test the 1D conv function.
% The method is to construct real vectors and compare the Matlab conv
% function on them with the quaternion conv function operating on
% quaternions with zero vector parts.

A = randn(1,13);
B = randn(1,5);
C = randn(13,1);

compare(conv(A, B), s(conv(quaternion(A), quaternion(B))), T, ...
        'quaternion/conv failed test 1')
    
compare(conv(B, A), s(conv(quaternion(B), quaternion(A))), T, ...
        'quaternion/conv failed test 2')

compare(conv(B, C), s(conv(quaternion(B), quaternion(C))), T, ...
        'quaternion/conv failed test 3')

compare(conv(C, B), s(conv(quaternion(C), quaternion(B))), T, ...
        'quaternion/conv failed test 4')

compare(conv(A, C), s(conv(quaternion(A), quaternion(C))), T, ...
        'quaternion/conv failed test 5')
    
% To test that the right coefficients of quaternion/conv works we supply a
% vector of ones for the first parameter.

compare(conv(A, B), s(conv({quaternion(ones(size(B))), quaternion(B)}, ...
        quaternion(A))), T, 'quaternion/conv failed test 6')
    
compare(conv(B, A), s(conv({quaternion(ones(size(A))), quaternion(A)}, ...
        quaternion(B))), T, 'quaternion/conv failed test 7')

compare(conv(B, C), s(conv({quaternion(ones(size(C))), quaternion(C)}, ...
        quaternion(B))), T, 'quaternion/conv failed test 8')

% Check that the special cases of scalar arguments work OK.

check(all(conv({qi, qk}, [qj qj]) == [-1 -1]), ...
                            'quaternion/conv failed test 9');

compare(conv(C, B), s(conv({quaternion(1), B}, ...
        quaternion(C))), T, 'quaternion/conv failed test 10')
    
compare(conv(B, C), s(conv({B, quaternion(1)}, ...
        quaternion(C))), T, 'quaternion/conv failed test 11')

% Check a row vector with a row vector.

compare(conv(B, B), s(conv(quaternion(B), quaternion(B))), ...
                   T, 'quaternion/conv failed test 12');

%% Test the 2D conv function.
% The method is to construct real matrices and compare the Matlab conv
% function on them with the quaternion conv function operating on
% quaternions with zero vector parts.

A = randn(4,3);
B = randn(6,7);
C = randn(3,4);

compare(conv2(A, B), s(conv2(quaternion(A), quaternion(B))), T, ...
        'quaternion/conv failed test 13')
    
compare(conv2(B, A), s(conv2(quaternion(B), quaternion(A))), T, ...
        'quaternion/conv failed test 14')

compare(conv2(B, C), s(conv2(quaternion(B), quaternion(C))), T, ...
        'quaternion/conv failed test 15')

compare(conv2(C, B), s(conv2(quaternion(C), quaternion(B))), T, ...
        'quaternion/conv failed test 16')

compare(conv2(A, C), s(conv2(quaternion(A), quaternion(C))), T, ...
        'quaternion/conv failed test 17')
    
% Test the function with left and right coefficient arrays.

compare(conv2(A .* C', B), s(conv2({quaternion(A), quaternion(C')}, ...
                                    quaternion(B))), T,...
        'quaternion/conv failed test 18')

% Check that the special case of scalar arguments works OK.

compare(B, s(conv2({quaternion(1), quaternion(1)}, ...
                    quaternion(B))), T, 'quaternion/conv failed test 19')

%% TODO Add test code for the convw, convw2 and mustard functions.

% Ideas - these should ideally be tested against the qfft and qfft2.
% However, when these tests are done, the FFT functions may not have been
% tested.

disp('Passed');


% $Id: test_conv.m 1004 2017-11-15 17:14:09Z sangwine $

