classdef IsotypicCommutant < replab.SubEquivariant
% Commutant of a harmonized isotypic component
%
% Matrices in this commutant space have the following form:
% $ X = \sum_i M_i \otimes D_i \otimes A_i $
% where $M_i$ represents the multiplicity space, $D_i$ is a constant representing the representation space,
% and $A_i$ encodes the division algebra (in the case of representation over C, it is trivial).
%
% Subclasses of this class implement specialized projection methods that take advantage of Schur's lemma. In the real case, specializations
% are provided for the three real division algebras, corresponding to real-type, complex-type and quaternion-type representations.

    properties (SetAccess = protected)
        divisionAlgebraDimension % (integer): Size of a block in the division algebra encoding
    end

    methods (Access = protected)

        function [M, D, A] = projectAndFactorFromParent_exact(self, X)
            P = self.repR.projection(type, 'exact');
            I = self.repR.injection(type, 'exact');
            [M, D, A] = self.projectAndFactor_exact(P * X * I);
        end

        function [M, D, A, err] = projectAndFactorFromParent_double_sparse(self, X)
            P = self.repR.projection(type, 'double/sparse');
            I = self.repR.injection(type, 'double/sparse');
            [M, D, A, err] = self.projectAndFactor_double_sparse(P * X * I);
        end

        function [M, D, A] = projectAndFactor_exact(self, X)
            error('Exact projection not implemented');
        end

        function [M, D, A, err] = projectAndFactor_double_sparse(self, X)
            error('Abstract');
        end

    end

    methods

        function self = IsotypicCommutant(isotypic, divisionAlgebraDimension)
            self@replab.SubEquivariant(isotypic.parent.commutant, isotypic, isotypic, 'commutant');
            self.divisionAlgebraDimension = divisionAlgebraDimension;
        end

    end

    methods % Implementations

        % Equivariant

        function [X1, err] = project(self, X, type)
            if nargin < 3
                type = 'double';
            end
            [M, D, A, err] = self.projectAndFactor(X, type);
            X1 = kron(M{1}, kron(D{1}, A{1}));
            for i = 2:length(M)
                X1 = X1 + kron(M{i}, kron(D{i}, A{i}));
            end
        end

    end


    methods % Projection

        function s = reducedBlockSize(self)
        % Returns the size of a commutant algebra element block, with the redundancy eliminated
        %
        % Returns:
        %   integer: Block size
            s = self.repR.multiplicity * self.divisionAlgebraDimension;
        end

        function [X1, err] = projectFromParent(self, X, type)
        % Projects the given matrix in the parent representation space into the commutant algebra
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix in the parent representation space to project
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   X1: double(\*,\*) or `.cyclotomic`(\*,\*)
        %     Projected matrix
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned projection to the invariant subspace in Frobenius norm
            if nargin < 3
                type = 'double';
            end
            [M, D, A, err] = self.projectAndFactorFromParent(self, X, type);
            X1 = kron(M{1}, kron(D{1}, A{1}));
            for i = 2:length(M)
                X1 = X1 + kron(M{i}, kron(D{i}, A{i}));
            end
        end

        function [M, D, A, err] = projectAndFactor(self, X, type)
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
        %   D: cell(1,\*) of double(\*,\*) or `.cyclotomic`(\*,\*)
        %     The constant part that only depends on the isotypic component structure
        %   A: cell(1,\*) of double(\*,\*) or `.cyclotomic`(\*,\*)
        %     The part that encodes the division algebra
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned projection to the invariant subspace in Frobenius norm
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
                case 'exact'
                  [M, D, A] = self.projectAndFactor_exact(X);
                  err = 0;
              case 'double'
                [M, D, A, err] = self.projectAndFactor_double_sparse(X);
                M = cellfun(@full, M, 'uniform', 0);
                D = cellfun(@full, D, 'uniform', 0);
                A = cellfun(@full, A, 'uniform', 0);
              case 'double/sparse'
                [M, D, A, err] = self.projectAndFactor_double_sparse(X);
              otherwise
                error('Unknown type %s', type);
            end
        end

        function [M, D, A, err] = projectAndFactorFromParent(self, X, type)
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
        %   D: cell(1,\*) of double(\*,\*) or `.cyclotomic`(\*,\*)
        %     The constant part that only depends on the isotypic component structure
        %   A: cell(1,\*) of double(\*,\*) or `.cyclotomic`(\*,\*)
        %     The part that encodes the division algebra
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned projection to the invariant subspace in Frobenius norm
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
                case 'exact'
                  [M, D, A] = self.projectAndFactorFromParent_exact(X);
                  err = 0;
              case 'double'
                [M, D, A, err] = self.projectAndFactorFromParent_double_sparse(X);
                M = cellfun(@full, M, 'uniform', 0);
                D = cellfun(@full, D, 'uniform', 0);
                A = cellfun(@full, A, 'uniform', 0);
              case 'double/sparse'
                [M, D, A, err] = self.projectAndFactorFromParent_double_sparse(X);
              otherwise
                error('Unknown type %s', type);
            end
        end

    end

end
