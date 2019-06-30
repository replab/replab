function test_qdft
% Test code for the discrete quaternion Fourier transform.

% This code tests the following functions:
%
%  qdft  qdft2
% iqdft iqdft2
%
% It also verifies indirectly many of the basic quaternion operations
% since the qdft depends on them.

% Copyright (c) 2005, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing discrete quaternion Fourier transforms ...')

T = 1e-12;

RA =        unit(quaternion(1,1,1));  % Real axis.
CA = complex(RA, quaternion(1,0,-1)); % Complex axis.

if  isreal(CA); error('Complex axis is not complex.'); end
if ~isreal(RA); error('Real axis is complex.'); end

% Test the one-dimensional qdft code for the case of row and column
% vectors. This is a simple test to check that the code handles the
% cases identically, but not a comprehensive test, since the tests on
% the two-dimensional case below also exercise the one-dimensional code.

% Define one real and one complex quaternion vector.

q = randq(1,10) .* randn(1,10);
b = complex(randq(1,10), randq(1,10)) .* randn(1,10);

% Tests 1 and 2. Verify that row and column vectors transform identically
% by taking the forward transform of a row vector and the inverse transform
% of the transposed spectrum. The result should be the transpose of the
% original vector.

compare(q.', iqdft(qdft(q, RA, 'L').', RA, 'L'), T, 'qdft failed test 1L.');
compare(q.', iqdft(qdft(q, RA, 'R').', RA, 'R'), T, 'qdft failed test 1R.');

compare(b.', iqdft(qdft(b, CA, 'L').', CA, 'L'), T, 'qdft failed test 2L.');
compare(b.', iqdft(qdft(b, CA, 'R').', CA, 'R'), T, 'qdft failed test 2R.');

% From here on we are testing the one-dimensional code using the
% two-dimensional code.

% Redefine q and b as a real and a complex quaternion matrix.

q = randq(10) .* randn(10);
b = complex(randq(10), randq(10)) .* randn(10);

% Test 3. Verify correct transform and inverse for a real quaternion
% array with a real quaternion axis.

compare(q, iqdft2(qdft2(q, RA, 'L'), RA, 'L'), T, 'qdft failed test 3L.');
compare(q, iqdft2(qdft2(q, RA, 'R'), RA, 'R'), T, 'qdft failed test 3R.');

% Test 4. Verify correct transform and inverse for a real quaternion
% array with a complex axis.

compare(q, iqdft2(qdft2(q, CA, 'L'), CA, 'L'), T, 'qdft failed test 4L.');
compare(q, iqdft2(qdft2(q, CA, 'R'), CA, 'R'), T, 'qdft failed test 4R.');

% Test 5. Verify correct transform and inverse for a complex quaternion
% array with a complex axis.

compare(b, iqdft2(qdft2(b, CA, 'L'), CA, 'L'), T, 'qdft failed test 5L.');
compare(b, iqdft2(qdft2(b, CA, 'R'), CA, 'R'), T, 'qdft failed test 5R.');

% Test 6. Verify correct transform and inverse for a complex quaternion
% array with a real axis.

compare(b, iqdft2(qdft2(b, RA, 'L'), RA, 'L'), T, 'qdft failed test 6L.');
compare(b, iqdft2(qdft2(b, RA, 'R'), RA, 'R'), T, 'qdft failed test 6R.');

% Verify the 1D qdft code against a matrix exponential DFT. This is a
% strong test since the latter is computed without quaternion arithmetic.
% It does depend on the adjoint function, so it also verifies that the
% adjoint function gives a correct result.

q = randq(1,10) .* randn(1,10);
b = complex(randq(1,10), randq(1,10)) .* randn(1,10);

RA = randv;
CA = unit(complex(randv, randv));

f = qdft(q, RA, 'L');
Q = matdft([f.w; f.x; f.y; f.z], adjoint(-RA, 'real'));
compare(q, quaternion(Q(1,:), Q(2,:), Q(3,:), Q(4,:)) ./ length(q), ...
                                                 T, 'qdft failed test 7.');

f = qdft(q, CA, 'L');
Q = matdft([f.w; f.x; f.y; f.z], adjoint(-CA, 'real'));
compare(q, quaternion(Q(1,:), Q(2,:), Q(3,:), Q(4,:)) ./ length(q), ...
                                                 T, 'qdft failed test 8.');
                                             
% TODO Consider whether the matrix exponential comparison can be extended
% to the two-dimensional case.

disp('Passed');

end

% Quaternion DFT based on matrix exponential, used in the last few tests
% above to provide an assurance that the QDFT code is correctly
% implementing the quaternion DFT (this is important since the QFFT
% functions are tested against the QDFT).
%
% Reference:
%
% Sangwine, S. J. and Ell, T. A., "Complex and Hypercomplex Discrete
% Fourier Transforms Based on Matrix Exponential Form of Euler's Formula",
% e-print arXiv:1001.4379, 25 January 2010, available at
% http://arxiv.org/abs/arxiv:1001.4379.

function F = matdft(f, J)
% f is an array of four rows, these corresponding to the w, x, y, z parts
% of a quaternion. J must be an adjoint matrix representation of the axis.
M = size(f, 2);
F = zeros(size(f));
for m = 1:M
    for u = 1:M
        F(:,u) = F(:,u) + ...
            expm(-J .* 2 .* pi .* (m-1) .* (u-1)./ M) * f(:,m);
    end
end
end

% $Id: test_qdft.m 1004 2017-11-15 17:14:09Z sangwine $
