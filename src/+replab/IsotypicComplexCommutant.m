classdef IsotypicComplexCommutant < replab.IsotypicCommutant

    methods

        function self = IsotypicComplexCommutant(isotypic)
            self = self@replab.IsotypicCommutant(isotypic, 2);
        end

    end

    methods (Static)

        function X = basisA
            X = eye(2);
        end

        function X = basisB
            X = [0 -1; 1 0];
        end

    end

    methods (Access = protected)

        function [A, B] = blockFromParent(self, X, type)
        % Changes the basis and projects a block on this isotypic component
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*)): Matrix to project on this commutant algebra in the basis of the original representation
        %   type ('double', 'double/sparse' or 'exact'): Type
        %
        % Returns
        % -------
        %   A:
        %    double(\*,\*) or `.cyclotomic`(\*,\*): The real part of the projected block
        %   B:
        %    double(\*,\*) or `.cyclotomic`(\*,\*): The imaginary part of the projected block
            m = self.repR.multiplicity;
            id = self.repR.irrepDimension;
            P = self.repR.projection(type);
            I = self.repR.injection(type);
            range = 1:id:m*id;
            A = P(range,:)*X*I(:,range) + P(range+1,:)*X*I(:,range+1);
            B = P(range+1,:)*X*I(:,range) - P(range,:)*X*I(:,range+1);
            for i = 3:2:id
                range = i:id:m*id;
                A = A + P(range,:)*X*I(:,range) + P(range+1,:)*X*I(:,range+1);
                B = B + P(range+1,:)*X*I(:,range) - P(range,:)*X*I(:,range+1);
            end
            A = A/id;
            B = B/id;
        end

        function [A, B] = block(self, X)
        % Returns the block of a matrix projected in the commutant algebra
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*)): Matrix to project on this commutant algebra
        %
        % Returns
        % -------
        %   A:
        %    double(\*,\*) or `.cyclotomic`(\*,\*): The real part of the projected block
        %   B:
        %    double(\*,\*) or `.cyclotomic`(\*,\*): The imaginary part of the projected block
            m = self.repR.multiplicity;
            id = self.repR.irrepDimension;
            range = 1:id:m*id;
            A = X(range, range) + X(range+1, range+1);
            B = X(range+1, range) - X(range, range+1);
            for i = 3:2:id
                r = i:id:m*id;
                A = A + X(range, range) + X(range+1, range+1);
                B = B + X(range+1, range) - X(range, range+1);
            end
            A = A/id;
            B = B/id;
        end

        function [M, D, A] = projectAndFactorFromParent_exact(self, X)
            [A, B] = self.blockFromParent(X);
            M = {A B};
            I = replab.cyclotomic.eye(self.repR.irrepDimension/2);
            D = {I I};
            A = {replab.cyclotomic.eye(2) replab.fromDoubles([0 -1; 1 0])};
        end

        function [M, D, A] = projectAndFactorFromParent_double_sparse(self, X)
            [A, B] = self.blockFromParent(X);
            M = {A B};
            I = speye(self.repR.irrepDimension/2);
            D = {I I};
            A = {speye(2) [0 -1; 1 0]};
        end

        function [M, D, A] = projectAndFactor_exact(self, X)
            [A, B] = self.block(X);
            M = {A B};
            I = replab.cyclotomic.eye(self.repR.irrepDimension/2);
            D = {I I};
            A = {replab.cyclotomic.eye(2) replab.fromDoubles([0 -1; 1 0])};
        end

        function [M, D, A] = projectAndFactor_double_sparse(self, X)
            [A, B] = self.block(X);
            M = {A B};
            I = speye(self.repR.irrepDimension/2);
            D = {I I};
            A = {speye(2) [0 -1; 1 0]};
        end

    end

end
