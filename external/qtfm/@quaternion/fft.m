function Y = fft(X)
% FFT Discrete Quaternion Fourier transform.
% (Quaternion overloading of standard Matlab function, but only one parameter.)
% (The parameters N and dim of the standard function are not yet implemented.)
% 
% This function implements a default quaternion FFT. See the related function
% QFFT, which implements transforms with left or right exponentials and a
% user-specified axis.

% Copyright (c) 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Implementation note: keeping QFFT separate means that the quaternion-specific
% parameters (axis and left/right) are kept separate from the Matlab standard
% parameters (N and dim) which might be added here at a later date.

narginchk(1, 1), nargoutchk(0, 1) 

Y = qfft(X, dft_axis(isreal(X)), 'L');

% $Id: fft.m 1004 2017-11-15 17:14:09Z sangwine $

