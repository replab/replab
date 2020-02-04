classdef IsotypicCommutant < replab.Equivariant

    properties
        divisionAlgebraDimension % integer: Size of a block in the division algebra encoding
    end

    methods

        function self = IsotypicCommutant(isotypic)
            self = self@replab.Equivariant(isotypic, isotypic, 'commutant');
        end

        function s = reducedBlockSize(self)
        % Returns the size of a commutant algebra element block, without repetition
            s = self.rep.multiplicity * self.divisionAlgebraDimension;
        end

        function block = projectAndReduce(self, X)
        % Projects the given matrix in the commutant algebra and removes its inherent redundancy
        %
        % Args:
        %   X (double): Matrix in the isotypic component space
        %
        % Returns:
        %   double: The corresponding block of size `reducedBlockSize`,
        %           removing the redundancy due to the irrep dimension
            error('Abstract');
        end

        function block = projectAndReduceFromParent(self, X)
        % Projects the given matrix given in the parent representation space and removes its redundancy
        %
        % Args:
        %   X (double): Matrix in the parent representation space
        %
        % Returns:
        %   double: The projected block of size `reducedBlockSize` corresponding
        %           to this isotypic component, having removed the redundancy due to the irrep dimension
            U = self.rep.U;
            block = self.projectAndReduce(U*X*U');
        end

    end

end
