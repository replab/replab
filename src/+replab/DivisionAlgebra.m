classdef DivisionAlgebra < replab.domain.VectorSpace
% Represents a division algebra over the reals, encoded as matrices
%
% It supports four division algebras/encodings:
%
% * ``'complex'`` : Complex numbers
% * ``'quaternion.rep'`` : The quaternion division algebra, with the encoding used for representations
% * ``'quaternion.equivariant'`` : The quaternion division algebra, with the encoding used for equivariant spaces
%
% The ``'complex'`` algebra represents a complex number $a + i b$ encoded in a ``2 x 2`` matrix algebra::
%
%   [ a -b
%     b  a ]
%
% It has dimension 2 and matrixSize 2 as well. Note that another encoding could be ``[a b; -b a]``; the encoding we choose seems more widespread but is
% otherwise arbitrary.
%
% The quaternion algebra/encodings represent quaternions of the form $q = a + i b + j c + k d$, where $i, j, k$ are the fundamental quaternion units, obeying the relations
% $i^2 = j^2 = k^2 = i j k = -1$. Thus, we need three matrices $I, J, K$ obeying those relations. If we further require the matrices to be sparse, there are
% 48 different possible encodings, see
% R. W. Farebrother, J. Groß, and S.-O. Troschke, "Matrix representation of quaternions", Linear Algebra and its Applications, vol. 362, pp. 251–255, Mar. 2003
% https://doi.org/10.1016/S0024-3795(02)00535-9
%
% Before introducing our choice, let us remark that the quaternion $q$ can be encoded used the complex matrix::
%
%   Quaternion a+i*b+j*c+k*d, encoding using a 2x2 complex matrix
%
%     [ a+i*b c+i*d
%      -c+i*d a-i*b ]
%
% Now, expanding each complex coefficient using the complex algebra real encoding described above, we obtain::
%
%   Quaternion equivariant encoding (quaternion.equivariant)
%
%     [  a -b -c -d
%        b  a  d  c
%       -c -d  a  b
%        d -c -b  a  ]
%
% This is the encoding that will be looked for when sampling from quaternionic-type commutant matrices, thus we named it that way.
%
% Now, quaternionic-type representations will have another encoding of the quaternion algebra that commutes with the equivariant encoding. It is::
%
%   Quaternion representation encoding (quaternion.rep)
%
%     [  a -b -c -d
%        b  a -d  c
%        c  d  a -b
%        d -c  b  a  ]

    properties (SetAccess = protected)
        name % ('complex', 'quaternion.rep', 'quaternion.equivariant'): Name of the division algebra
        indexMatrix % (integer(matrixSize,matrixSize)): Matrix of signed indices
        basis % (double(matrixSize,matrixSize,basis)): Division algebra basis
        dualBasis % (double(matrixSize,matrixSize,basis)): Dual basis under the product ``(D,B) = trace(D'*B)``
        dimension % (integer): Algebra dimension
        matrixSize % (integer): Row/column size of the matrix representation
    end

    methods

        function self = DivisionAlgebra(name)
            switch name
              case 'complex'
                indexMatrix = [1 -2; 2 1];
              case 'quaternion.equivariant'
                indexMatrix = [1 -2 3 -4; 2 1 4 3; -3 -4 1 2; 4 -3 -2 1];
              case 'quaternion.rep'
                indexMatrix = [1 -2 -3 -4; 2 1 -4 3; 3 4 1 -2; 4 -3 2 1];
              otherwise
                error('Unknown division algebra %s', name);
            end
            dimension = max(abs(indexMatrix(:)));
            matrixSize = size(indexMatrix, 1);
            basis = zeros(matrixSize, matrixSize, dimension);
            dualBasis = zeros(matrixSize, matrixSize, dimension);
            for i = 1:dimension
                basis(:,:,i) = (indexMatrix == i) - (indexMatrix == -i);
                dualBasis(:,:,i) = basis(:,:,i)/sum(sum(abs(indexMatrix) == i));
            end
            self.field = 'R';
            self.name = name;
            self.indexMatrix = indexMatrix;
            self.basis = basis;
            self.dualBasis = dualBasis;
            self.dimension = dimension;
            self.matrixSize = matrixSize;
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
                X1 = X1 + reshape(B(:,:,i), [ms*ms 1])*(reshape(D(:,:,i), [1 ms*ms])*X);
            end
            X1 = reshape(permute(reshape(X1, [ms ms d1 d2]), [1 3 2 4]), [ms*d1 ms*d2]);
        end

    end

    methods % Implementations

        function b = eqv(self, lhs, rhs)
            b = all(all(lhs == rhs));
        end

        function X = sample(self)
            y = randn(self.dimension, 1);
            X = reshape(self.basis, [self.matrixSize*self.matrixSize self.dimension]) * y;
            X = reshape(X, [self.matrixSize self.matrixSize]);
        end

    end

end
