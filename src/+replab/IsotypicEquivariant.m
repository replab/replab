classdef IsotypicEquivariant < replab.SubEquivariant
% Equivariant space between two harmonized isotypic components containing equivalent representations
%
% If the two isotypic components contains equivalent irreducible representations, the matrices in this equivariant
% space have the following form:
%
% $ X = \sum_i M_i \otimes R_i \otimes A_i $
%
% where
%
% - $M_i$ represents the multiplicity space,
% - $R_i$ is a constant matrix representing the representation space,
% - $A_i$ is a constant matrix encoding the division algebra.
%
% and `.isZero` is false.
%
% Otherwise, when the irreducible representations in both components are inequivalent,
% ``M`` is always empty, and the equivariant space contains a single element, the zero matrix; then `.isZero` is true.
%
% Example:
%   >>> S3 = replab.S(3);
%   >>> natDec = S3.naturalRep.decomposition;
%   >>> std1 = natDec.component(2);
%   >>> nat2Dec = S3.naturalRep.tensorPower(2).decomposition;
%   >>> std2 = nat2Dec.component(3);
%   >>> E = std1.isotypicEquivariantFrom(std2);
%   >>> X = E.sample;
%   >>> M = E.projectAndFactor(X);
%   >>> length(M)
%       1
%   >>> R = E.R('double');
%   >>> A = E.A('double');
%   >>> tol = 1e-14;
%   >>> norm(kron(M{1}, kron(R{1}, A{1})) - X, 'fro') <= tol
%       1
%   >>> norm(E.reconstruct(M) - X, 'fro') <= tol
%       1

    properties (SetAccess = protected)
        R_internal % (cell(1,\*) of double(\*,\*) or `.cyclotomic`(\*,\*)): Representation space basis
        A_internal % (cell(1,\*) of double(\*,\*) or `.cyclotomic`(\*,\*)): Division algebra basis
    end

    methods

        function self = IsotypicEquivariant(parent, repR, repC, special, R_internal, A_internal)
            self@replab.SubEquivariant(parent, repR, repC, special);
            self.R_internal = R_internal;
            self.A_internal = A_internal;
        end

    end

    methods % Implementations

        % Domain

        function l = laws(self)
            l = replab.laws.IsotypicEquivariantLaws(self);
        end

    end

    methods

        function d = divisionAlgebraDimension(self)
        % Returns the dimension of the division algebra encoded in this block
        %
        % As a special case, if the block is zero, the returned size is zero.
        %
        % Returns:
        %   integer: Dimension of the division algebra
            d = length(self.A_internal);
        end

        function b = isZero(self)
        % Returns whether this equivariant space contains only the zero matrix
        %
        % This happens when the isotypic components `.repR` and `.repC` correspond to inequivalent irreducible representations
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> rep = S3.naturalRep;
        %   >>> triv = rep.decomposition.component(1);
        %   >>> std = rep.decomposition.component(2);
        %   >>> E = triv.isotypicEquivariantFrom(std);
        %   >>> E.isZero
        %       1
        %
        % Returns:
        %   logical: True if the equivariant space is trivial
            error('Abstract');
        end

    end

    methods % Factorization

        function res = R(self, type)
        % Returns the representation space basis
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   cell(1,\*) of double(\*,\*) or `.cyclotomic`(\*,\*): The representation space basis
            if nargin < 2 || isempty(type)
                type = 'double';
            end
            res = cellfun(@(X) replab.numerical.convert(X, type), self.R_internal, 'uniform', 0);
        end

        function res = A(self, type)
        % Returns the division algebra basis
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   cell(1,\*) of double(\*,\*) or `.cyclotomic`(\*,\*): The division algebra basis
            if nargin < 2 || isempty(type)
                type = 'double';
            end
            res = cellfun(@(X) replab.numerical.convert(X, type), self.A_internal, 'uniform', 0);
        end

        function X = reconstruct(self, M, type)
        % Reconstructs an equivariant space element returned by a factorization method
        %
        % Args:
        %   M (cell(1,\*) of double(\*,\*) or cyclotomic(\*,\*)): Factorized element
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or cyclotomic(\*,\*): Equivariant matrix
            error('Abstract');
        end

    end

    methods % Projection

        function [M, err] = projectAndFactor(self, X, type)
        % Projects the given matrix in the commutant algebra and factors it
        %
        % It returns the decomposition of the projection.
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix in the isotypic component space to project
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   M: cell(1,\*) of double(\*,\*) or `.cyclotomic`(\*,\*)
        %     The part containing the degrees of freedom of the commutant algebra
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned projection to the invariant subspace in Frobenius norm
            error('Abstract');
        end

    end

    methods % Projection from parent space

        function [M, err] = projectAndFactorFromParent(self, X, type)
        % Projects the given matrix in the parent representation space into the commutant algebra and factors it
        %
        % It returns the decomposition of the projection.
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix in the parent representation space to project
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   M: cell(1,\*) of double(\*,\*) or `.cyclotomic`(\*,\*)
        %     The part containing the degrees of freedom of the commutant algebra
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned projection to the invariant subspace in Frobenius norm
            error('Abstract');
        end

    end

    methods (Static, Access = protected)

        function E = make_exact(parent, repR, repC, special)
            error('Not implemented'); % TODO
        end

        function E = make_double(parent, repR, repC, special, parentSample)
            if isempty(parentSample)
                sub = replab.SubEquivariant(parent, repR, repC, special);
                [X, err] = sub.sample; % TODO: proper error estimation here
            else
                X = repR.projection('double/sparse') * parentSample * repC.injection('double/sparse');
                assert(replab.globals.yolo);
                eR = repR.errorBound;
                eC = repC.errorBound;
                cR = repR.conditionNumberEstimate; % condition number of repR
                cC = repC.conditionNumberEstimate; % condition number of repC
                sX = replab.numerical.norm2UpperBound(X);
                err = sX*(eR*cC + cR*eC);
            end
            err = max(err, 10*sqrt(repR.dimension*repC.dimension)*eps);
            assert(err < 1e-10, 'Error to big');
            if norm(X, 'fro') < err
                E = [];
                return
            end
            if repR.overR || repR.irrep(1).frobeniusSchurIndicator == 1
                ird = repR.irrepDimension;
                R = replab.IsotypicEquivariant.kronSecondFactor(X, ird, ird);
                R_internal = {R};
                A_internal = {1};
            elseif repR.irrep(1).frobeniusSchurIndicator == -2 % quaternion-type
                ird = repR.irrepDimension;
                [A1, A2, A3, A4] = replab.domain.QuaternionTypeMatrices.basis('commutant');
                [X1, X2, X3, X4] = replab.domain.QuaternionTypeMatrices.fromMatrix(X, 'commutant');
                R1 = replab.IsotypicEquivariant.kronSecondFactor(X1, ird, ird);
                R2 = replab.IsotypicEquivariant.kronSecondFactor(X2, ird, ird);
                R3 = replab.IsotypicEquivariant.kronSecondFactor(X3, ird, ird);
                R4 = replab.IsotypicEquivariant.kronSecondFactor(X4, ird, ird);
                R1 = R1/norm(R1,'fro')*sqrt(ird);
                R2 = R2/norm(R2,'fro')*sqrt(ird);
                R3 = R3/norm(R3,'fro')*sqrt(ird);
                R4 = R4/norm(R4,'fro')*sqrt(ird);
                A_internal = {A1 A2 A3 A4};
                R_internal = {R1 R2 R3 R4};
            elseif repR.irrep(1).frobeniusSchurIndicator == 0 % complex-type
                ird = repR.irrepDimension;
                [A1, A2] = replab.domain.ComplexTypeMatrices.basis;
                [X1, X2] = replab.domain.ComplexTypeMatrices.fromMatrix(X);
                R1 = replab.IsotypicEquivariant.kronSecondFactor(X1, ird, ird);
                R2 = replab.IsotypicEquivariant.kronSecondFactor(X2, ird, ird);
                R1 = R1/norm(R1,'fro')*sqrt(ird);
                R2 = R2/norm(R2,'fro')*sqrt(ird);
                A_internal = {A1 A2};
                R_internal = {R1 R2};
            else
                error('Unknown type');
            end
            E = replab.equi.IsotypicEquivariant_nontrivial(parent, repR, repC, special, R_internal, A_internal);
            E
        end

    end

    methods (Static)

        function X1 = kronFirstFactor(X, X2)
        % Factors a Kronecker product, second step
        %
        % Assuming that X is approximately ``X = kron(X1, X2)``, with ``X2`` known,
        % this estimates the first factor ``X1`` and returns it.
        %
        % Args:
        %   X (double(\*,\*)): Matrix to factor
        %   X2 (double(\*,\*)): Second factor
        %
        % Returns:
        %   double(\*,\*): An estimation of the first factor in the Kronecker product
            F = sum(sum(conj(X2).*X2));
            nRows2 = size(X2, 1);
            nCols2 = size(X2, 2);
            nRows1 = size(X, 1) / nRows2;
            nCols1 = size(X, 2) / nCols2;
            X1 = zeros(nRows1, nCols1);
            for r = 1:nRows1
                Ir = (r-1)*nRows2+(1:nRows2);
                for c = 1:nCols1
                    Ic = (c-1)*nCols2+(1:nCols2);
                    M = X(Ir, Ic);
                    X1(r, c) = sum(sum(M.*conj(X2)))/F;
                end
            end
        end

        function X2 = kronSecondFactor(X, nRows2, nCols2)
        % Factors a Kronecker product, first step
        %
        % Assuming that X is approximately ``X = kron(X1, X2)``, with ``X2`` a matrix of size ``nRows2`` x ``nCols2``,
        % this estimates the factor ``X2`` and returns it.
        %
        % Args:
        %   X (double(\*,\*)): Matrix to factor
        %   nRows2 (integer): Number of rows of the second factor
        %   nCols2 (integer): Number of columns of the second factor
        %
        % Returns:
        %   double(\*,\*): An estimation of the second factor in the Kronecker product
            nRows1 = size(X, 1) / nRows2;
            nCols1 = size(X, 2) / nCols2;
            assert(round(nRows1) == nRows1 && round(nCols1) == nCols1);
            X2 = zeros(nRows2, nCols2);
            for r = 1:nRows1
                Ir = (r-1)*nRows2+(1:nRows2);
                for c = 1:nCols1
                    Ic = (c-1)*nCols2+(1:nCols2);
                    M = X(Ir, Ic);
                    X2 = X2 + M/sign(M(1,1));
                end
            end
        end

        function E = make(repR, repC, varargin)
        % Returns the space of equivariant linear maps between two isotypic components
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``repC.image(g) * X = X * repR.image(g)``
        %
        % Args:
        %   repR (`+replab.Isotypic`): Isotypic component on the target/row space
        %   repC (`+replab.Isotypic`): Isotypic component on the source/column space
        %
        % Keyword Args:
        %   special (charstring): Special structure, see help on `.Equivariant`
        %   type ('exact', 'double' or 'double/sparse'): Whether to obtain an exact equivariant space ('double' and 'double/sparse' are equivalent)
        %   parent (`.Equivariant`, optional): Equivariant space from ``repC.parent`` to ``repR.parent``
        %   parentSample (double(\*,\*), optional): A generic sample from the parent space
        %
        % Returns:
        %   `+replab.IsotypicEquivariant`: The equivariant vector space
            assert(isa(repR, 'replab.Isotypic'));
            assert(isa(repC, 'replab.Isotypic'));
            assert(repR.field == repC.field);
            args = struct('special', '', 'parent', [], 'type', 'double', 'parentSample', []);
            args = replab.util.populateStruct(args, varargin);
            parent = args.parent;
            if isempty(parent)
                switch args.special
                  case 'commutant'
                    parent = repR.parent.commutant(args.type);
                  case 'hermitian'
                    parent = repR.parent.hermitianInvariant(args.type);
                  case 'trivialRows'
                    parent = repR.parent.equivariantFrom(repC.parent, 'type', args.type);
                  case 'trivialCols'
                    parent = repR.parent.equivariantFrom(repC.parent, 'type', args.type);
                  case ''
                    parent = repR.parent.equivariantFrom(repC.parent, 'type', args.type);
                  otherwise
                    error('Invalid special structure');
                end
            end
            if repR.dimension == 0 || repC.dimension == 0
                E = [];
            elseif repR.irrepDimension ~= repC.irrepDimension
                E = [];
            elseif repR.irrep(1).frobeniusSchurIndicator ~= repC.irrep(1).frobeniusSchurIndicator
                E = [];
            else
                if strcmp(args.type, 'exact')
                    assert(repR.isExact);
                    assert(repC.isExact);
                    assert(parent.isExact);
                    E = replab.IsotypicEquivariant.make_exact(parent, repR, repC, args.special);
                elseif strcmp(args.type, 'double') || strcmp(args.type, 'double/sparse')
                    E = replab.IsotypicEquivariant.make_double(parent, repR, repC, args.special, args.parentSample);
                else
                    error('Wrong type value');
                end
            end
            if isempty(E)
                E = replab.equi.IsotypicEquivariant_trivial(parent, repR, repC, args.special);
            end
        end

    end

end
