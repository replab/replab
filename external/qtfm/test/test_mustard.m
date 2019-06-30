function test_mustard
% Test code for the quaternion mustard function. This also tests convw and
% convw2, since these are called by the mustard code. These tests should be
% done after the QFFT code has been tested, otherwise we are testing
% convolutions against untested QFFTs.

% Copyright (c) 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing mustard convolution (and wrap around convolutions) ...')

T = 1e-10;

% 1D tests.

f = randq(10, 1) .* randn(10, 1);
g = randq(10, 1) .* randn(10, 1);

mu = randv;

M = mustard(f, g, mu, 'L'); 
N = iqfft(qfft(f, mu, 'L') .* qfft(g, mu, 'L'), mu, 'L');

compare(M, N, T, '1D mustard convolution fails test 1')

M = mustard(f, g, mu, 'R'); 
N = iqfft(qfft(f, mu, 'R') .* qfft(g, mu, 'R'), mu, 'R');

compare(M, N, T, '1D mustard convolution fails test 2')

% 2D tests.

f = randq(5, 5) .* randn(5, 5);
g = randq(5, 5) .* randn(5, 5);

mu = randv;

M = mustard(f, g, mu, 'L'); 
N = iqfft2(qfft2(f, mu, 'L') .* qfft2(g, mu, 'L'), mu, 'L');

compare(M, N, T, '2D mustard convolution fails test 1')

M = mustard(f, g, mu, 'R'); 
N = iqfft2(qfft2(f, mu, 'R') .* qfft2(g, mu, 'R'), mu, 'R');

compare(M, N, T, '2D mustard convolution fails test 2')

disp('Passed');

end

% $Id: test_mustard.m 1004 2017-11-15 17:14:09Z sangwine $
