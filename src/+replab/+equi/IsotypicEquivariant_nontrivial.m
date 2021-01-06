classdef IsotypicEquivariant_nontrivial < replab.IsotypicEquivariant

    methods

        function self = IsotypicEquivariant_nontrivial(parent, repR, repC, special, R_internal, A_internal)
            self@replab.IsotypicEquivariant(parent, repR, repC, special, R_internal, A_internal);
        end

    end

    methods % Implementations

        % IsotypicEquivariant

        function b = isZero(self)
        % Returns whether this equivariant space contains only the zero matrix
        %
        % This happens when the isotypic components `.repR` and `.repC` correspond to inequivalent irreducible representations
        %
        % Returns:
        %   logical: True if the equivariant space is trivial
            b = false;
        end

        function X = reconstruct(self, M, type)
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            R = self.R(type);
            A = self.A(type);
            X = kron(M{1}, kron(R{1}, A{1}));
            for i = 2:length(M)
                X = X + kron(M{i}, kron(R{i}, A{i}));
            end
        end

        function [M, err] = projectAndFactor(self, X, type)
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
              case 'exact'
                M = self.projectAndFactor_exact(X);
                err = 0;
              case 'double'
                M = self.projectAndFactor_double_sparse(X);
                M = cellfun(@full, M, 'uniform', 0);
              case 'double/sparse'
                M = self.projectAndFactor_double_sparse(X);
              otherwise
                error('Unknown type %s', type);
            end
            if nargout > 1
                assert(replab.globals.yolo);
                eR = self.repR.errorBound;
                eC = self.repC.errorBound;
                cR = self.repR.conditionNumberEstimate; % condition number of repR
                cC = self.repC.conditionNumberEstimate; % condition number of repC
                sX = replab.numerical.norm2UpperBound(X);
                err = sX*(eR*cC + cR*eC);
            end
        end

        function [M, err] = projectAndFactorFromParent(self, X, type)
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
                case 'exact'
                  M = self.projectAndFactorFromParent_exact(X);
                  err = 0;
              case 'double'
                M = self.projectAndFactorFromParent_double_sparse(X);
                M = cellfun(@full, M, 'uniform', 0);
              case 'double/sparse'
                M = self.projectAndFactorFromParent_double_sparse(X);
              otherwise
                error('Unknown type %s', type);
            end
            if nargout > 1
                assert(replab.globals.yolo);
                eR = self.repR.errorBound;
                eC = self.repC.errorBound;
                cR = self.repR.conditionNumberEstimate; % condition number of repR
                cC = self.repC.conditionNumberEstimate; % condition number of repC
                sX = 0;
                R = self.R;
                A = self.A;
                X1 = kron(M{1}, kron(R{1}, A{1}));
                for i = 1:length(M)
                    sX = sX + replab.numerical.norm2UpperBound(M{i}) * replab.numerical.norm2UpperBound(R{i}) * replab.numerical.norm2UpperBound(A{i});
                end
                err = sX*(eR*cC + cR*eC);
            end
        end

    end

    methods (Access = protected)

        function M = projectAndFactor_exact(self, X)
            parentX = self.repR.injection('exact') * X * self.repC.projection('exact');
            M = self.projectAndFactorFromParent(parentX, 'exact');
        end

        function M = projectAndFactor_double_sparse(self, X)
            parentX = self.repR.injection('double/sparse') * X * self.repC.projection('double/sparse');
            M = self.projectAndFactorFromParent(parentX, 'double/sparse');
        end

        function M = projectAndFactorFromParent_exact(self, parentX)
            X = self.projectFromParent(parentX, 'exact');
            R = self.R('exact');
            if self.repR.overR || self.repR.irrep(1).frobeniusSchurIndicator == 1
                M = cell(1, 1);
                M{1} = replab.IsotypicEquivariant.kronFirstFactor(X, R{1});
            elseif repR.irrep(1).frobeniusSchurIndicator == -2 % quaternion-type
                [X1, X2, X3, X4] = replab.domain.QuaternionTypeMatrices.fromMatrix(X, 'commutant');
                M = cell(1, 4);
                M{1} = replab.IsotypicEquivariant.kronFirstFactor(X1, R{1});
                M{2} = replab.IsotypicEquivariant.kronFirstFactor(X2, R{2});
                M{3} = replab.IsotypicEquivariant.kronFirstFactor(X3, R{3});
                M{4} = replab.IsotypicEquivariant.kronFirstFactor(X4, R{4});
            elseif repR.irrep(1).frobeniusSchurIndicator == 0 % complex-type
                [X1, X2] = replab.domain.ComplexTypeMatrices.fromMatrix(X);
                M = cell(1, 2);
                M{1} = replab.IsotypicEquivariant.kronFirstFactor(X1, R{1});
                M{2} = replab.IsotypicEquivariant.kronFirstFactor(X2, R{2});
            end
        end

        function M = projectAndFactorFromParent_double_sparse(self, parentX)
            X = self.projectFromParent(parentX, 'double/sparse');
            R = self.R('double/sparse');
            if self.repR.overR || self.repR.irrep(1).frobeniusSchurIndicator == 1
                M = cell(1, 1);
                M{1} = replab.IsotypicEquivariant.kronFirstFactor(X, R{1});
            elseif repR.irrep(1).frobeniusSchurIndicator == -2 % quaternion-type
                [X1, X2, X3, X4] = replab.domain.QuaternionTypeMatrices.fromMatrix(X, 'commutant');
                M = cell(1, 4);
                M{1} = replab.IsotypicEquivariant.kronFirstFactor(X1, R{1});
                M{2} = replab.IsotypicEquivariant.kronFirstFactor(X2, R{2});
                M{3} = replab.IsotypicEquivariant.kronFirstFactor(X3, R{3});
                M{4} = replab.IsotypicEquivariant.kronFirstFactor(X4, R{4});
            elseif repR.irrep(1).frobeniusSchurIndicator == 0 % complex-type
                [X1, X2] = replab.domain.ComplexTypeMatrices.fromMatrix(X);
                M = cell(1, 2);
                M{1} = replab.IsotypicEquivariant.kronFirstFactor(X1, R{1});
                M{2} = replab.IsotypicEquivariant.kronFirstFactor(X2, R{2});
            end
        end

    end

end
