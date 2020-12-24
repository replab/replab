classdef IsotypicCommutant < replab.Equivariant
% Commutant of a harmonized isotypic component

    properties (SetAccess = protected)
        divisionAlgebraDimension % (integer): Size of a block in the division algebra encoding
    end

    methods (Access = protected)

        function part1 = projectAndFactor_exact(self, X)
            error('Exact projection not implemented');
        end

        function [part1, err] = projectAndFactor_double_sparse(self, X)
            error('Abstract');
        end

    end

    methods

        function self = IsotypicCommutant(isotypic, divisionAlgebraDimension)
            self = self@replab.Equivariant(isotypic, isotypic, 'commutant');
            self.divisionAlgebraDimension = divisionAlgebraDimension;
        end

        function s = reducedBlockSize(self)
        % Returns the size of a commutant algebra element block, without repetition
        %
        % Returns:
        %   integer: Block size
            s = self.repR.multiplicity * self.divisionAlgebraDimension;
        end

        function [part1, part2, err] = projectAndFactor(self, X, type)
        % Projects the given matrix in the commutant algebra and factors it
        %
        % It returns two arguments ``part1`` and ``part2`` such that the result of the projection
        % is ``kron(part1, part2)``; and ``part2`` depends only on the isotypic component structure.
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix in the isotypic component space to project
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   part1: double(\*,\*) or `.cyclotomic`(\*,\*)
        %     The block of size `reducedBlockSize`, removing the redundancy due to the irrep dimension
        %   part2: double(\*,\*) or `.cyclotomic`(\*,\*)
        %     The constant part that only depends on the isotypic component structure
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned projection to the invariant subspace in Frobenius norm
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
                case 'exact'
                  part1 = self.projectAndFactor_exact(X);
                  if nargout >= 2
                      part2 = replab.cyclotomic.eye(self.repR.irrepDimension);
                  end
                  if nargout >= 3
                      err = 0;
                  end
              case 'double'
                if nargout >= 3
                    [part1, err] = self.projectAndFactor_double_sparse(X);
                else
                    part1 = self.projectAndFactor_double_sparse(X);
                end
                part1 = full(part1);
                if nargout >= 2
                    part2 = eye(self.repR.irrepDimension)
                end
              case 'double/sparse'
                if nargout >= 3
                    [part1, err] = self.projectAndFactor_double_sparse(X);
                else
                    part1 = self.projectAndFactor_double_sparse(X);
                end
                if nargout >= 2
                    part2 = speye(self.repR.irrepDimension)
                end
              otherwise
                error('Unknown type %s', type);
            end
        end

        function [part1, part2, err] = projectAndFactorFromParent(self, X, type)
        % Projects the given matrix in the parent representation space into the commutant algebra and factors it
        %
        % It returns two arguments ``part1`` and ``part2`` such that the result of the projection
        % is ``kron(part1, part2)``; and ``part2`` depends only on the isotypic component structure.
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix in the parent representation space to project
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   part1: double(\*,\*) or `.cyclotomic`(\*,\*)
        %     The block of size `reducedBlockSize`, removing the redundancy due to the irrep dimension
        %   part2: double(\*,\*) or `.cyclotomic`(\*,\*)
        %     The constant part that only depends on the isotypic component structure
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned projection to the invariant subspace in Frobenius norm
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            if strcmp(type, 'exact')
                P = self.repR.projection(type, 'exact');
                I = self.repR.injection(type, 'exact');
            else
                P = self.repR.projection(type, 'double/sparse');
                I = self.repR.injection(type, 'double/sparse');
            end
            X1 = P*X*I;
            if strcmp(type, 'double')
                X1 = full(X1);
            end
            switch nargout
              case 2
                [part1, part2] = self.projectAndReduce(X1, type);
              case 3
                [part1, part2, err] = self.projectAndReduce(X1, type);
              otherwise
                part1 = self.projectAndReduce(X1, type);
            end
        end

    end

end
