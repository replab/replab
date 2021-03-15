classdef DivisionAlgebra < replab.domain.VectorSpace
% Describes an encoding of a division algebra using real or complex matrices
%
% The three division algebras supported are the real numbers (R), complex numbers (C) and quaternions (H).
% They can be encoded using real matrix blocks (R) or complex matrix blocks (C).
%
% A complex scalar is written ``s = a + i b`` where ``a`` and ``b` are reals and ``i`` is the imaginary unit.
% A quaternion scalar is written ``s = a + i b + j c + k d``, where ``i, j , k`` are the fundamental quaternion units,
% obeying the relations ``i^2 = j^2 = k^2 = i j k = -1``.
%
% It supports the following encodings:
%
% * ``'R->C'``: A real scalar is encoded as one complex scalar.
% * ``'C->R'``: A complex scalar is encoded as a 2x2 real matrix block.
% * ``'H->C'``: A quaternion scalar is encoded as a 2x2 complex matrix block.
% * ``'H->R:rep'`` : A quaternion scalar is encoded as a 4x4 real matrix block, using the ``rep`` encoding.
% * ``'H->R:equivariant'`` : A quaternion scalar is encoded as a 4x4 real matrix block, using the ``equivariant`` encoding.
%
% The ``'C->R'`` encoding represents a complex number $a + i b$ as::
%
%   [ a -b
%     b  a ]
%
% It has dimension 2 and matrixSize 2 as well. Note that another encoding could be ``[a b; -b a]``;
% the encoding we choose seems more widespread but is otherwise arbitrary.
%
% For the quaternions, we need three matrices $I, J, K$ obeying those relations. If we further require the matrices
% to be sparse, there are 48 different possible encodings, see
% R. W. Farebrother, J. Groß, and S.-O. Troschke, "Matrix representation of quaternions", Linear Algebra and its Applications, vol. 362, pp. 251–255, Mar. 2003
% https://doi.org/10.1016/S0024-3795(02)00535-9
%
% Before introducing our choice, let us remark that the quaternion $q$ can be encoded used the complex matrix::
%
%   Quaternion a+i*b+j*c+k*d, encoding using a 2x2 complex matrix (H->C)
%
%     [ a+i*b c+i*d
%      -c+i*d a-i*b ]
%
% Now, expanding each complex coefficient using the complex algebra real encoding described above, we obtain::
%
%   Quaternion equivariant encoding (H->R:equivariant)
%
%     [  a -b  c -d
%        b  a  d  c
%       -c -d  a  b
%        d -c -b  a  ]
%
% This is the encoding that will be looked for when sampling from quaternionic-type commutant matrices, thus we named it that way.
%
% Now, quaternionic-type representations will have another encoding of the quaternion algebra that commutes with the equivariant encoding. It is::
%
%   Quaternion representation encoding (H->R:rep)
%
%     [  a -b -c -d
%        b  a -d  c
%        c  d  a -b
%        d -c  b  a  ]

    properties (SetAccess = protected)
        name % ('complex', 'quaternion.rep', 'quaternion.equivariant'): Name of the division algebra
        basis % (double(matrixSize,matrixSize,dimension)): Division algebra basis
        dualBasis % (double(matrixSize,matrixSize,dimension)): Dual abasis under the product ``(D,B) = trace(D'*B)``
        indexMatrix % (double(matrixSize,matrix)): Index matrix
    end

    methods

        function self = DivisionAlgebra(name)
            switch name
              case 'complex'
                name = 'C->R';
              case 'quaternion.equivariant'
                name = 'H->R:equivariant';
              case 'quaternion.rep'
                name = 'H->R:rep';
            end
            switch name
              case 'R->C'
                indexMatrix = [1];
              case 'C->R'
                indexMatrix = [1 -2; 2 1];
              case 'H->C'
                indexMatrix = [1+2i 3+4i; -3+4i 1-2i];
              case 'H->R:equivariant'
                indexMatrix = [1 -2 3 -4; 2 1 4 3; -3 -4 1 2; 4 -3 -2 1];
              case 'H->R:rep'
                indexMatrix = [1 -2 -3 -4; 2 1 -4 3; 3 4 1 -2; 4 -3 2 1];
              otherwise
                error('Unknown division algebra %s', name);
            end
            dim = max([real(indexMatrix(:)); imag(indexMatrix(:))]);
            ms = size(indexMatrix, 1);
            basis = zeros(ms, ms, dim);
            dualBasis = zeros(ms, ms, dim);
            for j = 1:dim
                realPart = (real(indexMatrix) == j) - (real(indexMatrix) == -j);
                imagPart = (imag(indexMatrix) == j) - (imag(indexMatrix) == -j);
                basis(:,:,j) = realPart + 1i*imagPart;
                dualBasis(:,:,j) = basis(:,:,j);
                f = trace(dualBasis(:,:,j)'*basis(:,:,j));
                dualBasis(:,:,j) = dualBasis(:,:,j)/f;
            end
            self.field = 'R';
            self.name = name;
            self.basis = basis;
            self.dualBasis = dualBasis;
            self.indexMatrix = indexMatrix;
        end

        function m = matrixSize(self)
            m = size(self.basis, 1);
        end

        function d = dimension(self)
            d = size(self.basis, 3);
        end

        function X1 = decode(self, X)
        % Decodes the division algebra elements present in a matrix encoding
        %
        % This method accepts a matrix of size ``(ms*d1) x (ms*d2)`` as argument, where
        % ``ms`` is `.matrixSize` and ``d1``, ``d2`` are integers.
        %
        % Args:
        %   X (double(ms*d1,ms*d2)): (Block) matrix to decode
        %
        % Returns:
        %   double(d1,d2,d): Decoded elements where ``d`` is `.dimension`
            ms = self.matrixSize;
            d = self.dimension;
            D = self.dualBasis;
            d1 = size(X, 1)/ms;
            d2 = size(X, 2)/ms;
            X1 = zeros(d1, d2, d);
            X = reshape(permute(reshape(X, [ms d1 ms d2]), [1 3 2 4]), [ms*ms d1*d2]);
            for i = 1:d
                X1(:,:,i) = real(reshape(reshape(conj(D(:,:,i)), [1 ms*ms])*X, [d1 d2]));
            end
        end

        function X1 = encode(self, X)
        % Takes a matrix of division algebra elements and encodes them in a matrix
        %
        % This method accepts a tensor of size ``d1 x d2 x d`` as argument, where
        % ``d`` is `.dimension` and ``d1``, ``d2`` are integers.
        %
        % Args:
        %   X (double(d1,d2,d)): Tensor encoding the matrix
        %
        % Returns:
        %   double(ms*d1, ms*d2): Encoded matrix
            ms = self.matrixSize;
            d = self.dimension;
            D = self.dualBasis;
            B = self.basis;
            d1 = size(X, 1);
            d2 = size(X, 2);
            X1 = zeros(ms*ms, d1*d2);
            for i = 1:d
                X1 = X1 + reshape(B(:,:,i), [ms*ms 1])*reshape(X(:,:,i), [1 d1*d2]);
            end
            X1 = reshape(permute(reshape(X1, [ms ms d1 d2]), [1 3 2 4]), [ms*d1 ms*d2]);
        end

        function X1 = project(self, X)
        % Projects the given matrix in this division algebra
        %
        % This method accepts a matrix of size ``(ms*d1) x (ms*d2)`` as argument, where
        % ``ms`` is `.matrixSize` and ``d1``, ``d2`` are integers. The projection is then done
        % on each ``ms x ms`` block.
        %
        % Args:
        %   X (double(ms*d1,ms*d2)): (Block) matrix to project
        %
        % Returns:
        %   double(ms*d1, ms*d2): Projected matrix
            ms = self.matrixSize;
            d = self.dimension;
            D = self.dualBasis;
            B = self.basis;
            d1 = size(X, 1)/ms;
            d2 = size(X, 2)/ms;
            X1 = zeros(ms*ms, d1*d2);
            X = reshape(permute(reshape(X, [ms d1 ms d2]), [1 3 2 4]), [ms*ms d1*d2]);
            for i = 1:d
                X1 = X1 + reshape(B(:,:,i), [ms*ms 1])*real(reshape(conj(D(:,:,i)), [1 ms*ms])*X);
            end
            X1 = reshape(permute(reshape(X1, [ms ms d1 d2]), [1 3 2 4]), [ms*d1 ms*d2]);
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, lhs, rhs)
            b = all(all(lhs == rhs));
        end

        function l = laws(self)
            l = replab.laws.DivisionAlgebraLaws(self);
        end

        function X = sample(self)
            y = randn(self.dimension, 1);
            X = reshape(self.basis, [self.matrixSize*self.matrixSize self.dimension]) * y;
            X = reshape(X, [self.matrixSize self.matrixSize]);
        end

    end

end
