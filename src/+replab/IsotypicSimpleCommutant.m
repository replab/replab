classdef IsotypicSimpleCommutant < replab.IsotypicCommutant

    methods

        function self = IsotypicSimpleCommutant(isotypic)
            self = self@replab.IsotypicCommutant(isotypic, 1);
        end

    end

    methods (Access = protected)

        function X1 = blockFromParent(self, X, type)
        % Changes the basis and projects a block on this isotypic component
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*)): Matrix to project on this commutant algebra in the basis of the original representation
        %   type ('double', 'double/sparse' or 'exact'): Type
        %
        % Returns:
        %   double(\*,\*) or `.cyclotomic`(\*,\*): Block corresponding to the isotypic component
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

        function [M, D, A] = projectAndFactorFromParent_exact(self, X)
            M = {self.blockFromParent(X, 'exact')};
            D = {replab.cyclotomic.eye(self.repR.irrepDimension)};
            A = {replab.cyclotomic.eye(1)};
        end

        function [M, D, A, err] = projectAndFactorFromParent_double_sparse(self, X)
            M = {self.blockFromParent(X, 'double/sparse')};
            D = {speye(self.repR.irrepDimension)};
            A = {1};
            err = inf;
        end

        function [M, D, A] = projectAndFactor_exact(self, X)
            M = {self.block(X)};
            D = {replab.cyclotomic.eye(self.repR.irrepDimension)};
            A = {replab.cyclotomic.eye(1)};
        end

        function [M, D, A, err] = projectAndFactor_double_sparse(self, X)
            M = {self.block(X)};
            D = {speye(self.repR.irrepDimension)};
            A = {1};
            err = inf;
        end

    end

end
