classdef IsotypicComplexCommutant < replab.IsotypicCommutant
    
    methods
        
        function self = IsotypicSimpleCommutant(isotypic)
            self = self@replab.IsotypicCommutant(isotypic);
            self.divisionAlgebraDimension = 1;
        end

        function [A B] = block(self, X)
        % Returns the block of a matrix projected in the commutant algebra
        %
        % Args:
        %   X (double): Matrix to project on this commutant algebra
        %
        % Returns
        % -------
        %   A:
        %    double: The real part of the projected block
        %   B:
        %    double: The imaginary part of the projected block
            m = self.rep.multiplicity;
            id = self.rep.irrepDimension;
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
        
        function X = projectAndReduce(self, X)
            [A B] = self.block(X);
            X = kron(A, eye(2)) + kron(B, [0 -1; 1 0]);
        end
        
        function X = project(self, X)
            id = self.rep.irrepDimension;
            [A B] = self.block(X);
            X = kron(A, eye(id)) + kron(B, kron(eye(id/2), [0 -1; 1 0]));
        end
        
    end
    
end
