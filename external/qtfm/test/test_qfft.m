function test_qfft
% Test code for the fast quaternion Fourier transform.

% This code tests the following functions:
%
%  qfft  qfft2
% iqfft iqfft2

% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% We have to test all of the above functions separately, since they do not
% call each other, as is the case with the corresponding dqft functions. In
% addition, we need to verify the results of the FFTs against the DFT code
% in the corresponding qdft functions, e.g. for qfft2, this is qdft2. For
% each transform we need to verify that it inverts correctly for all four
% combinations of real/complex data and real/complex axis, and for left and
% right exponentials, and we need to verify each against the corresponding
% DFT code. We do this with the default axes defined in the private
% function dft_axis, and with random axes, to check that the decompositions
% used in the qfft code are correct.

disp('Testing fast quaternion Fourier transforms ...')

% Define one real and one complex quaternion array.

q = randq(10) .* randn(10);
b = complex(randq(10), randq(10)) .* randn(10);

T = 1e-10;

% Construct a set of test axes, some fixed, some random.

RA = [unit(quaternion(1,1,1)), qi, qj, qk, randv, randv, randv, randv];
CA = [complex(quaternion(1,1,1), quaternion(0,1,-1)), ...
      unit(complex(quaternion(0,1,1), quaternion(1,1,1))), ...
      unit(qi + qj + qk .* 1i), unit(qi + randv .* 1i), ...
      unit(complex(randv, randv)), unit(complex(randv, randv)), ...
      unit(complex(randv, randv)), unit(complex(randv, randv))];
  
assert(length(RA) == length(CA));

if  isreal(CA); error('Complex axis is not complex.'); end
if ~isreal(RA); error('Real axis is complex.'); end

for k = 1:length(RA)
    
    R = RA(k);
    C = CA(k);
    
    % 1D FFT code against its own inverse ....
    
    % Verify correct transform and inverse for a real quaternion array with a
    % real quaternion axis.
    
    compare(q, iqfft(qfft(q, R, 'L'), R, 'L'), T, ...
        ['qfft failed test 1L with axis: ', char(R)]);
    compare(q, iqfft(qfft(q, R, 'R'), R, 'R'), T, ...
        ['qfft failed test 1R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a real quaternion array with a
    % complex axis.
    
    compare(q, iqfft(qfft(q, C, 'L'), C, 'L'), T, ...
        ['qfft failed test 2L with axis: ', char(R)]);
    compare(q, iqfft(qfft(q, C, 'R'), C, 'R'), T, ...
        ['qfft failed test 2R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a complex axis.
    
    compare(b, iqfft(qfft(b, C, 'L'), C, 'L'), T, ...
        ['qfft failed test 3L with axis: ', char(R)]);
    compare(b, iqfft(qfft(b, C, 'R'), C, 'R'), T, ...
        ['qfft failed test 3R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a real axis.
    
    compare(b, iqfft(qfft(b, R, 'L'), R, 'L'), T, ...
        ['qfft failed test 4L with axis: ', char(R)]);
    compare(b, iqfft(qfft(b, R, 'R'), R, 'R'), T, ...
        ['qfft failed test 4R with axis: ', char(R)]);
    
    % 2D FFT code against its own inverse ....
    
    % Verify correct transform and inverse for a real quaternion array with a
    % real quaternion axis.
    
    compare(q, iqfft2(qfft2(q, R, 'L'), R, 'L'), T, ...
        ['qfft2 failed test 1L with axis: ', char(R)]);
    compare(q, iqfft2(qfft2(q, R, 'R'), R, 'R'), T, ...
        ['qfft2 failed test 1R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a real quaternion array with a
    % complex axis.
    
    compare(q, iqfft2(qfft2(q, C, 'L'), C, 'L'), T, ...
        ['qfft2 failed test 2L with axis: ', char(R)]);
    compare(q, iqfft2(qfft2(q, C, 'R'), C, 'R'), T, ...
        ['qfft2 failed test 2R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a complex axis.
    
    compare(b, iqfft2(qfft2(b, C, 'L'), C, 'L'), T, ...
        ['qfft2 failed test 3L with axis: ', char(R)]);
    compare(b, iqfft2(qfft2(b, C, 'R'), C, 'R'), T, ...
        ['qfft2 failed test 3R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a real axis.
    
    compare(b, iqfft2(qfft2(b, R, 'L'), R, 'L'), T, ...
        ['qfft2 failed test 4L with axis: ', char(R)]);
    compare(b, iqfft2(qfft2(b, R, 'R'), R, 'R'), T, ...
        ['qfft2 failed test 4R with axis: ', char(R)]);
    
    % 1D FFT code against DFT inverse ....
    
    % Verify correct transform and inverse for a real quaternion array with a
    % real quaternion axis.
    
    compare(q, iqdft(qfft(q, R, 'L'), R, 'L'), T, ...
        ['qfft failed test 5L with axis: ', char(R)]);
    compare(q, iqdft(qfft(q, R, 'R'), R, 'R'), T, ...
        ['qfft failed test 5R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a real quaternion array with a
    % complex axis.
    
    compare(q, iqdft(qfft(q, C, 'L'), C, 'L'), T, ...
        ['qfft failed test 6L with axis: ', char(R)]);
    compare(q, iqdft(qfft(q, C, 'R'), C, 'R'), T, ...
        ['qfft failed test 6R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a complex axis.
    
    compare(b, iqdft(qfft(b, C, 'L'), C, 'L'), T, ...
        ['qfft failed test 7L with axis: ', char(R)]);
    compare(b, iqdft(qfft(b, C, 'R'), C, 'R'), T, ...
        ['qfft failed test 7R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a real axis.
    
    compare(b, iqdft(qfft(b, R, 'L'), R, 'L'), T, ...
        ['qfft failed test 8L with axis: ', char(R)]);
    compare(b, iqdft(qfft(b, R, 'R'), R, 'R'), T, ...
        ['qfft failed test 8R with axis: ', char(R)]);
    
    % 2D FFT code against DFT inverse ....
    
    % Verify correct transform and inverse for a real quaternion array with a
    % real quaternion axis.
    
    compare(q, iqdft2(qfft2(q, R, 'L'), R, 'L'), T, ...
        ['qfft2 failed test 5L with axis: ', char(R)]);
    compare(q, iqdft2(qfft2(q, R, 'R'), R, 'R'), T, ...
        ['qfft2 failed test 5R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a real quaternion array with a
    % complex axis.
    
    compare(q, iqdft2(qfft2(q, C, 'L'), C, 'L'), T, ...
        ['qfft2 failed test 6L with axis: ', char(R)]);
    compare(q, iqdft2(qfft2(q, C, 'R'), C, 'R'), T, ...
        ['qfft2 failed test 6R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a complex axis.
    
    compare(b, iqdft2(qfft2(b, C, 'L'), C, 'L'), T, ...
        ['qfft2 failed test 7L with axis: ', char(R)]);
    compare(b, iqdft2(qfft2(b, C, 'R'), C, 'R'), T, ...
        ['qfft2 failed test 7R with axis: ', char(R)]);
    
    % Verify correct transform and inverse for a complex quaternion array with
    % a real axis.
    
    compare(b, iqdft2(qfft2(b, R, 'L'), R, 'L'), T, ...
        ['qfft2 failed test 8L with axis: ', char(R)]);
    compare(b, iqdft2(qfft2(b, R, 'R'), R, 'R'), T, ...
        ['qfft2 failed test 8R with axis: ', char(R)]);

end

disp('Passed');

% $Id: test_qfft.m 1004 2017-11-15 17:14:09Z sangwine $
