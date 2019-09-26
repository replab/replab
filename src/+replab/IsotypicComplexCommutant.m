classdef IsotypicComplexCommutant < replab.IsotypicCommutant
    
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
            block = zeros(m, m);
            for i = 1:cd
                block = block + X(i:id:m*id, i:id:m*id);
            end
            block = block/id;
        end
        
        function X = projectAndReduce(self, X)
            X = self.block(X);
        end
        
        function X = project(self, X)
            X = kron(self.block(X), eye(self.rep.irrepDimension));
        end
        
    end
    
end
