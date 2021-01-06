classdef ComplexTypeMatrices < replab.domain.VectorSpace
% Represents a vector space of real matrices that encode complex entries
%
% A complex number is written ``c = a + b * 1i``, and can be represented
% by a matrix::
%
%   [a -b
%    b  a]
%
% There is an alternative complex encoding obtained by changing the sign of ``b``,
% corresponding to the complex conjugation.

    properties (SetAccess = protected)
        nR % (integer): Row size, must be a multiple of 2
        nC % (integer): Column size, must be a multiple of 2
    end

    methods

        function self = ComplexTypeMatrices(nR, nC)
            assert(mod(nR, 2) == 0);
            assert(mod(nC, 2) == 0);
            self.nR = nR;
            self.nC = nC;
            self.field = 'R';
        end

        %% Str methods

        function s = headerStr(self)
            s = sprintf('%d x %d real matrices encoding complex coefficient blocks', self.nR, self.nC);
        end

        %% Domain methods

        function b = eqv(self, X, Y)
            b = all(all(X == Y));
        end

        function X = sample(self)
            A = randn(self.nR/2, self.nC/2);
            B = randn(self.nR/2, self.nC/2);
            X = replab.domain.ComplexTypeMatrices.toMatrix(A, B)
        end

    end

    methods (Static)

        function [basisA, basisB] = basis
            basisA = [ 1  0
                       0  1];
            basisB = [ 0 -1
                       1  0];
        end

        function M = project(M)
        % Projects a generic matrix
            [A B] = replab.domain.ComplexTypeMatrices.fromMatrix(M);
            M = replab.domain.ComplexTypeMatrices.toMatrix(A, B);
        end

        function [A, B] = fromMatrix(M)
        % Given a matrix encoding complex blocks, returns the complex elements
        %
        % Args:
        %   M (double(\*,\*)): Matrix part of a `+replab.+domain.ComplexTypeMatrices` domain
        %
        % Returns
        % -------
        %   A: double(\*,\*)
        %     Real part
        %   B: double(\*,\*)
        %     Imaginary part
            nR = size(M, 1);
            nC = size(M, 2);
            assert(mod(nR, 2) == 0);
            assert(mod(nC, 2) == 0);
            A = (M(1:2:nR, 1:2:nC) + M(2:2:nR, 2:2:nC))/2;
            B = (M(2:2:nR, 1:2:nC) - M(1:2:nR, 2:2:nC))/2;
        end

        function M = toMatrix(A, B)
        % Returns the complex type matrix encoding the given complex coefficients
        %
        % The complex number is implicitly given as C = A + B i
        %
        % Args:
        %   A (double(\*,\*)): Real part
        %   B (double(\*,\*)): Imaginary part
        %
        % Returns:
        %   double(\*,\*): The matrix encoding the complex coefficient blocks
            [basisA, basisB] = replab.domain.ComplexTypeMatrices.basis;
            if issparse(A) && issparse(B)
                basisA = sparse(basisA);
                basisB = sparse(basisB);
            end
            M = kron(A, basisA) + kron(B, basisB);
        end

    end

end
