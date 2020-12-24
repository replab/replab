classdef IsotypicSimpleCommutant < replab.IsotypicCommutant

    methods (Access = protected)

        function X1 = blockFromParent(self, X, type)
        % Changes the basis and projects a block on this isotypic component
        %
        % Args:
        %   X (double(\*,\*)): Matrix to project on this commutant algebra in the basis of the original representation
        %
        % Returns:
        %   double(\*,\*): Block corresponding to the isotypic component
            m = self.repR.multiplicity;
            id = self.repR.irrepDimension;
            P = self.repR.projection(type);
            I = self.repR.injection(type);
            range = 1:id:m*id;
            X1 = P(range,:)*X*I(:,range);
            for i = 2:id
                range = i:id:m*id;
                X1 = X1 + P(range,:)*X*I(:,range);
            end
            X1 = X1/id;
        end

        function X1 = block(self, X)
        % Returns the block of a matrix projected in the commutant algebra
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*)): Matrix to project on this commutant algebra
        %
        % Returns:
        %   double(\*,\*) or `.cyclotomic`(\*,\*): The projected block
            m = self.repR.multiplicity;
            id = self.repR.irrepDimension;
            X1 = X(1:id:m*id, 1:id:m*id);
            for i = 2:id
                X1 = X1 + X(i:id:m*id, i:id:m*id);
            end
            X1 = X1/id;
        end

        function part1 = projectAndFactorFromParent_exact(self, X)
            part1 = self.blockFromParent(X, 'exact');
        end

        function [part1, err] = projectAndFactorFromParent_double_sparse(self, X)
            part1 = self.blockFromParent(X, 'double/sparse');
            err = NaN;
        end

        function part1 = projectAndFactor_exact(self, X)
            part1 = self.block(X);
        end

        function [part1, err] = projectAndFactor_double_sparse(self, X)
            part1 = self.block(X);
            err = NaN;
        end

        function X1 = project_exact(self, X)
            X1 = kron(self.block(X), replab.cyclotomic.eye(self.repR.irrepDimension));
        end

        function [X1, err] = project_double_sparse(self, X)
            X1 = kron(self.block(X), speye(self.repR.irrepDimension));
        end

    end

    methods

        function self = IsotypicSimpleCommutant(isotypic)
            self = self@replab.IsotypicCommutant(isotypic, 1);
        end

    end

end
