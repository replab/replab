classdef IsotypicSimpleCommutant < replab.IsotypicCommutant
    
    methods
        
        function self = IsotypicSimpleCommutant(isotypic)
            self = self@replab.IsotypicCommutant(isotypic);
            self.divisionAlgebraDimension = 1;
        end
        
        function X1 = block(self, X)
        % Returns the block of a matrix projected in the commutant algebra
        %
        % Args:
        %   X (double): Matrix to project on this commutant algebra
        %
        % Returns:
        %   double: The projected block
            m = self.rep.multiplicity;
            id = self.rep.irrepDimension;
            X1 = zeros(m, m);
            for i = 1:id
                X1 = X1 + X(i:id:m*id, i:id:m*id);
            end
            X1 = X1/id;
        end
        
        function X1 = blockFromParent(self, X) 
        % Changes the basis and projects a block on this isotypic component
        %
        % Args:
        %   X (double): Matrix to project on this commutant algebra in the basis of the original representation
        %
        % Returns:
        %   double: Block corresponding to the isotypic component
            m = self.rep.multiplicity;
            id = self.rep.irrepDimension;
            U = self.rep.U;
            X1 = zeros(m, m);
            for i = 1:id
                U1 = U(i:id:m*id, :);
                X1 = X1 + U1*X*U1';
            end
            X1 = X1/id;
        end
        
        function X1 = projectAndReduce(self, X)
            X1 = self.block(X);
        end
        
        function X1 = projectAndReduceFromParent(self, X)
            X1 = self.blockFromParent(X);
        end
        
        function X1 = project(self, X)
            X1 = kron(self.block(X), eye(self.rep.irrepDimension));
        end
        
    end
    
end
