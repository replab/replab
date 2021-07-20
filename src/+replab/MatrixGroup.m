classdef FiniteMatrixGroup < replab.FiniteGroup

    properties (SetAccess = protected)
        matrixSize % (integer): Matrix size
        elements % (cyclo.CycloSet): List of elements (Java object)
    end

    methods

        function self = FiniteMatrixGroup(matrixSize, generators)
            self.identity = replab.cyclotomic.eye(matrixSize);
        end

    end

    methods % Implementations

        function xInv = inverse(self, x)
            xInv = inv(x);
        end

        function z = compose(self, x, y)
            z = x * y;
        end

    end

end
