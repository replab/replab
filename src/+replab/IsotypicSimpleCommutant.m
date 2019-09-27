classdef IsotypicSimpleCommutant < replab.IsotypicCommutant
    
    methods
        
        function self = IsotypicSimpleCommutant(isotypic)
            self = self@replab.IsotypicCommutant(isotypic);
            self.divisionAlgebraDimension = 1;
        end

        function block = projectAndReduceFromParent(self, X)
        % Projects the given matrix given in the parent representation space and removes its redundancy
        %
        % Args:
        %   X (double): Matrix in the parent representation space
        %
        % Returns:
        %   double: The projected block of size `self.reducedBlockSize` corresponding
        %           to this isotypic component, having removed the redundancy due to the irrep dimension
            U = self.rep.U;
            block = self.projectAndReduce(U*X*U');
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
