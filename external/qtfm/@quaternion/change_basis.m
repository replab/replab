function P = change_basis(Q, B)
% CHANGE_BASIS changes the basis of the quaternion Q, to the basis B.
% Q may be a vector or matrix of quaternions, or a scalar quaternion.

% Copyright (c) 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% TODO This could be generalised to the 4D full quaternion case too, either
% using a 4x4 matrix for B, or by permitting B to be a 4-vector of
% quaternions. The latter would be a neater parameter profile. The choice
% between 3D and 4D would then be made by whether the B parameter was pure
% or full. However changing B would necessitate changes elsewhere. The
% matrix form for B permits a simple check on orthogonality, and also
% simple inversion of the basis.

narginchk(2, 2), nargoutchk(0, 1)

if ~isa(Q, 'quaternion')
    error('Q must be a quaternion, or a vector or matrix of quaternions.')
end

% Verify that B is an orthonormal basis.

if any(size(B) ~= [3, 3])
    error('The basis, B, must be a 3 by 3 matrix.');
end

if max(max(B * B.' - eye(3))) > 10 * eps
    warning('QTFM:inaccuracy', ...
            'The basis matrix is not accurately orthogonal.');
end

% Construct three pure quaternions from B.

V1 = quaternion(B(1,1), B(1,2), B(1,3));
V2 = quaternion(B(2,1), B(2,2), B(2,3));
V3 = quaternion(B(3,1), B(3,2), B(3,3));

% Change the basis of Q. This is done by resolving the vector part of Q
% into the directions of the three basis vectors V1, V2 and V3.

P = Q; % Create P as a copy of Q. This means if Q is full, so will P be,
       % and the scalar part will carry over from Q. If Q is pure, the
       % scalar part of P will be empty, as it should be.

P.x = scalar_product(Q, V1);
P.y = scalar_product(Q, V2);
P.z = scalar_product(Q, V3);

% $Id: change_basis.m 1004 2017-11-15 17:14:09Z sangwine $
