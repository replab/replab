classdef PermutationFiniteSet < replab.FiniteSet

    methods % Properties

        function d = domainSize(self)
            d = self.type.domainSize;
        end

    end

end
