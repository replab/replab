classdef IsotypicComplexCommutant < replab.IsotypicCommutant

    methods

        function self = IsotypicComplexCommutant(isotypic)
            self = self@replab.IsotypicCommutant(isotypic, 2);
        end

        function [A B] = block(self, X)
        % Returns the block of a matrix projected in the commutant algebra
        %
        % Args:
        %   X (double(\*,\*)): Matrix to project on this commutant algebra
        %
        % Returns
        % -------
        %   A:
        %    double(\*,\*): The real part of the projected block
        %   B:
        %    double(\*,\*): The imaginary part of the projected block
            m = self.repR.multiplicity;
            id = self.repR.irrepDimension;
            A = zeros(m, m);
            B = zeros(m, m);
            for i = 1:2:id
                r = i:id:m*id;
                A = A + X(r, r) + X(r+1, r+1);
                B = B + X(r+1, r) - X(r, r+1);
            end
            A = A/id;
            B = B/id;
        end

        function [A B] = blockFromParent(self, X)
        % Changes the basis and projects a block on this isotypic component
        %
        % Args:
        %   X (double(\*,\*)): Matrix to project on this commutant algebra in the basis of the original representation
        %
        % Returns
        % -------
        %   A:
        %    double(\*,\*): The real part of the projected block
        %   B:
        %    double(\*,\*): The imaginary part of the projected block
            m = self.repR.multiplicity;
            id = self.repR.irrepDimension;
            Er = self.repR.E_internal;
            Br = self.repR.B_internal;
            A = zeros(m, m);
            B = zeros(m, m);
            for i = 1:2:id
                range = shift+(i:id:m*id);
                A = A + Er(range,:)*X*Br(:,range) + Er(range+1,:)*X*Br(:,range+1);
                B = B + Er(range+1,:)*X*Br(:,range) - Er(range,:)*X*Br(:,range+1);
            end
            A = A/id;
            B = B/id;
        end

        function X1 = projectAndReduceFromParent(self, X)
            [A B] = self.blockFromParent(X);
            X1 = kron(A, eye(2)) + kron(B, [0 -1; 1 0]);
        end

        function X1 = projectAndReduce(self, X)
            [A B] = self.block(X);
            X1 = kron(A, eye(2)) + kron(B, [0 -1; 1 0]);
        end

        function [X1 err] = project(self, X)
            id = self.repR.irrepDimension;
            [A B] = self.block(X);
            X1 = kron(A, eye(id)) + kron(B, kron(eye(id/2), [0 -1; 1 0]));
            err = NaN;
        end

    end

end
