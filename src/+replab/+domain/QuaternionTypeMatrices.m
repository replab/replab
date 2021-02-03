classdef QuaternionTypeMatrices < replab.domain.VectorSpace
% Represents a vector space of real matrices that encode quaternion entries
%
% A quaternion is written q = a + b*i + c*j + d*k, with
%     i^2 = j^2 = k^2 = ijk = -1.
%
% Such a quaternion q can be encoded in a matrix algebra as such::
%
%   Group algebra
%   [a -b -c -d
%    b  a -d  c
%    c  d  a -b
%    d -c  b  a]
%
% There are 48 such encodings, see https://en.wikipedia.org/wiki/Quaternion#Matrix_representations
%
% The "group algebra" encoding is the one in John Vince, Quaternions for Computer Graphics, Springer 2011, page 67
% also given in Maehara and Murota, doi:10.1007/s13160-010-0007-8
%
% Commutant spaces, however, obey a different algebra, which commutes with the group algebra::
%
%   Commutant algebra
%   [ a -b  c -d
%     b  a  d  c
%    -c -d  a  b
%     d -c -b  a]
%
%
% Example:
%   >>> A = randn(2,2); B = randn(2,2); C = randn(2,2); D = randn(2,2);
%   >>> X1 = replab.domain.QuaternionTypeMatrices.toMatrix(A, B, C, D, 'group');
%   >>> X2 = replab.domain.QuaternionTypeMatrices.toMatrix(A, B, C, D, 'commutant');
%   >>> [A1, B1, C1, D1] = replab.domain.QuaternionTypeMatrices.fromMatrix(X1, 'group');
%   >>> [A2, B2, C2, D2] = replab.domain.QuaternionTypeMatrices.fromMatrix(X2, 'commutant');
%   >>> all(all(A == A1)) && all(all(A == A2))
%       1
%   >>> all(all(B == B1)) && all(all(B == B2))
%       1
%   >>> all(all(C == C1)) && all(all(C == C2))
%       1
%   >>> all(all(D == D1)) && all(all(D == D2))
%       1

    properties (SetAccess = protected)
        nR % (integer): Row size, must be a multiple of 4
        nC % (integer): Column size, must be a multiple of 4
        type % ('group' or 'commutant'): Encoding type
    end

    methods

        function self = QuaternionTypeMatrices(nR, nC, type)
        % Constructs a vector space of real matrices encoding quaternion-valued matrices
        %
        % Args:
        %   nR (integer): Row size, must be a multiple of 4
        %   nC (integer): Column size, must be a multiple of 4
        %   type ('group' or 'commutant'): Encoding type
            assert(mod(nR, 4) == 0);
            assert(mod(nC, 4) == 0);
            assert(strcmp(type, 'group') || strcmp(type, 'commutant'));
            self.field = 'R';
            self.nR = nR;
            self.nC = nC;
            self.type = type;
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = sprintf('%d x %d real matrices encoding quaternion coefficient blocks', self.nR, self.nC);
        end

        % Domain

        function b = eqv(self, X, Y)
            b = all(all(X == Y));
        end

        function X = sample(self)
            A = randn(self.nR/4, self.nC/4);
            B = randn(self.nR/4, self.nC/4);
            C = randn(self.nR/4, self.nC/4);
            D = randn(self.nR/4, self.nC/4);
            X = replab.domain.QuaternionTypeMatrices.toMatrix(A, B, C, D, type)
        end

    end

    methods

        function M = project(self, M)
        % Projects a generic matrix
        %
        % Args:
        %   M (double(\*,\*)): Matrix to project
        %
        % Returns:
        %   double(\*,\*): Matrix projected on the matrix subspace
            [A, B, C, D] = replab.domain.QuaternionTypeMatrices.fromMatrix(M, self.type);
            M = replab.domain.QuaternionTypeMatrices.toMatrix(A, B, C, D, self.type);
        end

    end

    methods (Static)

        function M = projectType(M, type)
        % Projects a generic matrix
        %
        % Args:
        %   M (double(\*,\*)): Matrix to project
        %   type ('group' or 'commutant'): Encoding type
        %
        % Returns:
        %   double(\*,\*): Matrix projected on the matrix subspace
            [A, B, C, D] = replab.domain.QuaternionTypeMatrices.fromMatrix(M, type);
            M = replab.domain.QuaternionTypeMatrices.toMatrix(A, B, C, D, type);
        end

        function [A, B, C, D] = fromMatrix(M, type)
        % Given a matrix encoding quaternion coefficient blocks, returns the quaternion elements
        %
        % Args:
        %   M (double(\*,\*)): Matrix part of a `.QuaternionTypeMatrices` domain
        %   type ('group' or 'commutant'): Encoding type
        %
        % Returns
        % -------
        %   A: double(\*,\*)
        %     Real part
        %   B: double(\*,\*)
        %     Pure quaternion part 'i'
        %   C: double(\*,\*)
        %     Pure quaternion part 'j'
        %   D: double(\*,\*)
        %     Pure quaternion part 'k'
            nR = size(M, 1);
            nC = size(M, 2);
            assert(mod(nR, 4) == 0);
            assert(mod(nC, 4) == 0);
            switch type
              case 'group'
                A = (M(1:4:nR, 1:4:nC) + M(2:4:nR, 2:4:nC) + M(3:4:nR, 3:4:nC) + M(4:4:nR, 4:4:nC))/4;
                B = (M(2:4:nR, 1:4:nC) - M(1:4:nR, 2:4:nC) + M(4:4:nR, 3:4:nC) - M(3:4:nR, 4:4:nC))/4;
                C = (M(3:4:nR, 1:4:nC) - M(1:4:nR, 3:4:nC) + M(2:4:nR, 4:4:nC) - M(4:4:nR, 2:4:nC))/4;
                D = (M(4:4:nR, 1:4:nC) + M(3:4:nR, 2:4:nC) - M(2:4:nR, 3:4:nC) - M(1:4:nR, 4:4:nC))/4;
              case 'commutant'
                A = (M(1:4:nR, 1:4:nC) + M(2:4:nR, 2:4:nC) + M(3:4:nR, 3:4:nC) + M(4:4:nR, 4:4:nC))/4;
                B = (M(2:4:nR, 1:4:nC) - M(1:4:nR, 2:4:nC) - M(4:4:nR, 3:4:nC) + M(3:4:nR, 4:4:nC))/4;
                C = -(M(3:4:nR, 1:4:nC) + M(1:4:nR, 3:4:nC) + M(2:4:nR, 4:4:nC) - M(4:4:nR, 2:4:nC))/4;
                D = (M(4:4:nR, 1:4:nC) - M(3:4:nR, 2:4:nC) + M(2:4:nR, 3:4:nC) - M(1:4:nR, 4:4:nC))/4;
              otherwise
                error('Unknown type');
            end
        end

        function [basisA, basisB, basisC, basisD] = basis(type)
        % Returns a basis of a quaternion encoding using real matricesx
        %
        % Args:
        %   type ('group' or 'commutant'): Encoding type
        %
        % Returns
        % -------
        %   basisA: double(\*,\*)
        %     Basis of the real part
        %   B: double(\*,\*)
        %     Basis of the pure quaternion part 'i'
        %   C: double(\*,\*)
        %     Basis of the pure quaternion part 'j'
        %   D: double(\*,\*)
        %     Basis of the pure quaternion part 'k'
            switch type
              case 'group'
                basisA = [ 1  0  0  0
                           0  1  0  0
                           0  0  1  0
                           0  0  0  1];
                basisB = [ 0 -1  0  0
                           1  0  0  0
                           0  0  0 -1
                           0  0  1  0];
                basisC = [ 0  0 -1  0
                           0  0  0  1
                           1  0  0  0
                           0 -1  0  0];
                basisD = [ 0  0  0 -1
                           0  0 -1  0
                           0  1  0  0
                           1  0  0  0];
              case 'commutant'
                basisA = [ 1  0  0  0
                           0  1  0  0
                           0  0  1  0
                           0  0  0  1];
                basisB = [ 0 -1  0  0
                           1  0  0  0
                           0  0  0  1
                           0  0 -1  0];
                basisC = [ 0  0  1  0
                           0  0  0  1
                          -1  0  0  0
                           0 -1  0  0];
                basisD = [ 0  0  0 -1
                           0  0  1  0
                           0 -1  0  0
                           1  0  0  0];
              otherwise
                error('Unknown type');
            end

        end

        function M = toMatrix(A, B, C, D, type)
        % Returns the quaternion matrix encoding the given quaternion coefficients
        %
        % The quaternion is implicitly given as Q = A + B i + C j + D k
        %
        % Args:
        %   A (double(\*,\*)): Real part
        %   B (double(\*,\*)): Pure quaternion 'i' part
        %   C (double(\*,\*)): Pure quaternion 'j' part
        %   D (double(\*,\*)): Pure quaternion 'k' part
        %
        % Returns:
        %   double(\*,\*): The matrix encoding the quaternion coefficient blocks
            [basisA, basisB, basisC, basisD] = replab.domain.QuaternionTypeMatrices.basis(type);
            if issparse(A) && issparse(B) && issparse(C) && issparse(D)
                basisA = sparse(basisA);
                basisB = sparse(basisB);
                basisC = sparse(basisC);
                basisD = sparse(basisD);
            end
            M = kron(A, basisA) + kron(B, basisB) + kron(C, basisC) + kron(D, basisD);
        end

    end

end
