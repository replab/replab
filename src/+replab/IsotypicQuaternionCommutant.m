classdef IsotypicQuaternionCommutant < replab.IsotypicCommutant

    methods

        function self = IsotypicQuaternionCommutant(isotypic)
            self = self@replab.IsotypicCommutant(isotypic, 4);
        end

    end

    methods (Access = protected)

        function [A, B, C, D] = block(self, X)
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
        %    double(\*,\*) or `.cyclotomic`(\*,\*): The 'i' part of the projected block
        %   C:
        %    double(\*,\*) or `.cyclotomic`(\*,\*): The 'j' part of the projected block
        %   D:
        %    double(\*,\*) or `.cyclotomic`(\*,\*): The 'k' part of the projected block
            m = self.repR.multiplicity;
            id = self.repR.irrepDimension;
            r1 = 1:id:m*id;
            r2 = r1+1;
            r3 = r1+2;
            r4 = r1+3;
            A = X(r1,r1) + X(r2,r2) + X(r3,r3) + X(r4,r4);
            B = X(r2,r1) - X(r1,r2) - X(r4,r3) + X(r3,r4);
            C = X(r3,r1) - X(r1,r3) - X(r2,r4) + X(r4,r2);
            D = X(r4,r1) - X(r3,r2) + X(r2,r3) - X(r1,r4);
            for i = 5:4:id
                r1 = i:id:m*id;
                r2 = (i:id:m*id)+1;
                r3 = (i:id:m*id)+2;
                r4 = (i:id:m*id)+3;
                A = A + X(r1,r1) + X(r2,r2) + X(r3,r3) + X(r4,r4);
                B = B + X(r2,r1) - X(r1,r2) - X(r4,r3) + X(r3,r4);
                C = C + X(r3,r1) - X(r1,r3) - X(r2,r4) + X(r4,r2);
                D = D + X(r4,r1) - X(r3,r2) + X(r2,r3) - X(r1,r4);
            end
            A = A/id;
            B = B/id;
            C = C/id;
            D = D/id;
        end

        function [A, B, C, D] = blockFromParent(self, X, type)
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
        %    double(\*,\*) or `.cyclotomic`(\*,\*): The 'i' part of the projected block
        %   C:
        %    double(\*,\*) or `.cyclotomic`(\*,\*): The 'j' part of the projected block
        %   D:
        %    double(\*,\*) or `.cyclotomic`(\*,\*): The 'k' part of the projected block
            m = self.repR.multiplicity;
            id = self.repR.irrepDimension;
            P = self.repR.projection(type);
            I = self.repR.injection(type)
            R = 1:id:m*id;
            A = P(R  ,:)*X*I(:, R  ) + P(R+1,:)*X*I(:, R+1) + P(R+2,:)*X*I(:, R+2) + P(R+3,:)*X*I(:, R+3);
            B = P(R+1,:)*X*I(:, R  ) - P(R  ,:)*X*I(:, R+1) - P(R+3,:)*X*I(:, R+2) + P(R+2,:)*X*I(:, R+3);
            C = P(R+2,:)*X*I(:, R  ) - P(R  ,:)*X*I(:, R+2) - P(R+1,:)*X*I(:, R+3) + P(R+3,:)*X*I(:, R+1);
            D = P(R+3,:)*X*I(:, R  ) - P(R+2,:)*X*I(:, R+1) + P(R+1,:)*X*I(:, R+2) - P(R  ,:)*X*I(:, R+3);
            % shape of things that commute with our representation quaternion encoding
            % [ a -b -c -d
            %   b  a  d -c
            %   c -d  a  b
            %   d  c -b  a]
            for i = 5:4:id
                R = i:id:m*id;
                A = A + P(R  ,:)*X*I(:, R  ) + P(R+1,:)*X*I(:, R+1) + P(R+2,:)*X*I(:, R+2) + P(R+3,:)*X*I(:, R+3);
                B = B + P(R+1,:)*X*I(:, R  ) - P(R  ,:)*X*I(:, R+1) - P(R+3,:)*X*I(:, R+2) + P(R+2,:)*X*I(:, R+3);
                C = C + P(R+2,:)*X*I(:, R  ) - P(R  ,:)*X*I(:, R+2) - P(R+1,:)*X*I(:, R+3) + P(R+3,:)*X*I(:, R+1);
                D = D + P(R+3,:)*X*I(:, R  ) - P(R+2,:)*X*I(:, R+1) + P(R+1,:)*X*I(:, R+2) - P(R  ,:)*X*I(:, R+3);
            end
            A = A/id;
            B = B/id;
            C = C/id;
            D = D/id;
        end

        function [M, D, A] = projectAndFactorFromParent_exact(self, X)
            [A, B, C, D] = self.blockFromParent(X);
            M = {A B C D};
            I = replab.cyclotomic.eye(self.repR.irrepDimension/4);
            D = {I I I I};
            A = cellfun(@(M) replab.cyclotomic.fromDoubles(M), self.basis, 'uniform', 0);
        end

        function [M, D, A, err] = projectAndFactorFromParent_double_sparse(self, X)
            [A, B, C, D] = self.blockFromParent(X);
            M = {A B C D};
            I = speye(self.repR.irrepDimension/4);
            D = {I I I I};
            A = self.basis;
            err = inf;
        end

        function [M, D, A] = projectAndFactor_exact(self, X)
            [A, B, C, D] = self.block(X);
            M = {A B C D};
            I = replab.cyclotomic.eye(self.repR.irrepDimension/4);
            D = {I I I I};
            A = cellfun(@(M) replab.cyclotomic.fromDoubles(M), self.basis, 'uniform', 0);
        end

        function [M, D, A, err] = projectAndFactor_double_sparse(self, X)
            [A, B, C, D] = self.block(X);
            M = {A B C D};
            I = speye(self.repR.irrepDimension/4);
            D = {I I I I};
            A = self.basis;
            err = inf;
        end

    end

    methods (Static)

        function B = basis
            B = {eye(4) ...
                 [ 0 -1  0  0
                   1  0  0  0
                   0  0  0  1
                   0  0 -1  0] ...
                 [ 0  0 -1  0
                   0  0  0 -1
                   1  0  0  0
                   0  1  0  0] ...
                 [ 0  0  0 -1
                   0  0  1  0
                   0 -1  0  0
                   1  0  0  0]};
        end

    end

end
