% Example script to calculate the vector phase correlation described in:
%
% Sangwine, S. J., Ell, T. A. and Moxey C. E., Vector Phase Correlation,
% Electronics Letters, 37(25), 6 December 2001, 1513-5.
% http://dx.doi.org/10.1049/el:20011035
%
% This script is not intended to explain how images are represented or read
% into quaternion arrays. For this, see the accompanying script
% colour_edge_detector.m.
%
% This script reads the two files noisy_image1.png and noisy_image2.png
% from the author's website. This is to avoid distributing the files with
% the toolbox itself as this would make the distribution zip file much
% larger and most users will not need the images.
%
% Copyright (c) : Steve Sangwine, May 2007

% Read the two input images, display them on screen and save them to disk.

disp('Reading two images from http://privatewww.essex.ac.uk/ ...');

n1 = imread('http://privatewww.essex.ac.uk/~sjs/images/noisy_image1.png');
n2 = imread('http://privatewww.essex.ac.uk/~sjs/images/noisy_image2.png');

imwrite(n1, 'noisy_image1.png');
imwrite(n2, 'noisy_image2.png');

figure(1)
image(n1)
axis off
axis image

figure(2)
image(n2)
axis off
axis image

f = convert(quaternion(n1(:,:,1), ...
                       n1(:,:,2), ...
                       n1(:,:,3)), 'double') ./ 256;
g = convert(quaternion(n2(:,:,1), ...
                       n2(:,:,2), ...
                       n2(:,:,3)), 'double') ./ 256;

disp('Calculating vector phase correlation ...');

mu = unit(quaternion(1,1,1)); % Transform axis, a unit pure quaternion.

% Calculate the quantities needed to evaluate equation 4 in the paper. The
% Matlab implementations of the FFTs differ in scale factor from the FFTs
% in equations 2 and 3 of the paper (because they follow the Matlab
% convention of applying the scale factor only to the inverse transform,
% rather than splitting it equally between the forward and inverse) so we
% have to make an adjustment for this:

MN  = prod(size(f)); % This gives the number of pixels in each image, which
                     % is the product MN in equations 2 and 3.

FL  =  qfft2(f, mu, 'L') ./ sqrt(MN); % This is F subscript L in the paper.
FIL = iqfft2(f, mu, 'L') .* sqrt(MN); % This is F subscript -L.
GR  =  qfft2(g, mu, 'R') ./ sqrt(MN); % This is G subscript R.

% Now decompose GR into components parallel and perpendicular to the
% transform axis, mu.

GRpar  = (GR - mu .* GR .* mu) ./ 2; % Parallel/perpendicular
GRperp =  GR - GRpar;                % decomposition.

RR = conj(FL) .* GRpar + conj(FIL) .* GRperp;

% Now we have the hypercomplex cross-power spectrum RR. We compute the
% hypercomplex phase correlation directly according to equation 6 in the
% paper. Note that the unit function from the QTFM toolbox divides each
% element of RR by its modulus. We don't have to write this out explicitly.

p = iqfft2(unit(RR), mu, 'R') .* sqrt(MN);

% Now display p as in Figure 3 in the paper. (The figure in the paper was
% created with Maple, rather than Matlab, so the appearance will be
% somewhat different.)

pcs = abs(p); % Phase correlation surface.

figure(3)
mesh(pcs)

disp('Mean value of phase correlation surface is:')
disp(mean(mean(pcs)))

[C, I] = max(pcs); % Find the maximum element in the phase correlation
[D, J] = max(C);   % surface, and its index value.

% The displayed value for the peak in the phase correlation surface differs
% slightly from that given in the paper. This is almost certainly due to
% the fact that the two PNG image files have a gamma value of 1, and Matlab
% probably assumes 2.2. Therefore the numerical values of the pixels used
% in the computation here are not the same as those used in the software
% that was used to compute the values in the paper (this software used the
% gamma values found in the file, rather than assuming a standard value of
% 2.2).

disp(['Peak value of phase correlation surface is:', D])
disp(D)

disp('Maximum value in phase correlation surface occurs at:')
disp([I(J), J])

% $Id: vector_phase_correlation.m 1004 2017-11-15 17:14:09Z sangwine $
