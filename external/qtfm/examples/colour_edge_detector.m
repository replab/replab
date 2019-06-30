% Example script to calculate the colour image edge detector described in:
%
% Sangwine, S. J.,
% Colour image edge detector based on quaternion convolution,
% Electronics Letters, 34(10), May 14 1998, 969-971.
% http://dx.doi.org/10.1049/el:19980697
%
% This script is not intended to demonstrate best use of Matlab - the
% intention is to demonstrate some of the features of the QTFM toolbox
% as an example of how the toolbox can be used to program a quaternion
% algorithm at a high level using vectorized coding, and also to show how
% to convert Matlab image arrays to quaternion matrices and vice versa.
%
% The script reads the lena image from an external website. The URL is
% given below. If access to this is not possible, replace the URL in the
% imread statement below with a suitable local filename and path, assuming
% you have a local copy of the lena image file.
%
% Copyright (c) : Steve Sangwine, May 2007

% Read the lena image into an array. Because this is a colour image the
% array will have three dimensions, the last of which has three index
% values.

url = 'http://sipi.usc.edu/database/misc/4.2.04.tiff';

disp(['Reading the lena image from ', url]);

qlena = cast(imreadq(url), 'double') ./ 255;

% Save the lena image after conversion to a quaternion array. It can then
% be reloaded later without accessing the external URL.

disp('Saving the quaternion lena image to local file qlena.mat')

save qlena.mat qlena

disp('Processing image with colour edge detector ...');

% Convert the array into a pure quaternion matrix, with the RGB components
% in the x, y and z components of the quaternion matrix. At the same time,
% convert the quaternion components into normalised double form, rather
% than the uint8 form read in from the file.

% qlena = convert(quaternion(lena(:,:,1), ...
%                            lena(:,:,2), ...
%                            lena(:,:,3)), 'double') ./ 256;
%                        
% Define a 3x3 quaternion mask pair as defined in the paper cited above.

mu = unit(quaternion(1,1,1)); % Equivalent to quaternion(1,1,1) ./ sqrt(3);

R = exp(mu .* pi ./ 4) ./ sqrt(6); % This is a direct coding of equation 2
                                   % in the paper.

% Create the left mask. Notice that we can supply a 3x1 matrix of zeros to
% the quaternion constructor to construct the (quaternion) zeros in the
% middle row, and we can apply the quaternion conjugate function to a 1x3
% matrix [R R R] to construct the bottom row, and join all three rows
% together using standard Matlab operators (which are overloaded to handle
% quaternions).

left = [R R R; quaternion(zeros(1, 3)); conj([R R R])]; % The left mask.

right = conj(left); % Notice how we apply the conjugate operation to the
                    % whole left mask in order to create the right mask.
                    
% Now convolve the lena image with the mask. This uses the quaternion
% overloading of the standard conv2 function. Warning: this is very slow -
% an alternative is to use a quaternion Fourier transform to compute the
% convolution much faster.

tic;
fqlena = conv2({left, right}, qlena); % Notice the argument form {} to
toc;                                  % indicate a left/right convolution
                                      % (special to the quaternion version
                                      % of conv2).
                                      
result = v(fqlena); % Ignore the scalar part of the result.

image(result) % Display the filtered result in an image window.
axis off
axis image

disp('Writing edge detector result to file ...');
imwrite(result, 'result.tif');

disp('Finished.');

% $Id: colour_edge_detector.m 1004 2017-11-15 17:14:09Z sangwine $

