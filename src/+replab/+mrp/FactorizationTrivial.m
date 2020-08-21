classdef FactorizationTrivial < replab.mrp.Factorization
% Computes the factorization of a permutation group elements in its generators

    methods

        function self = FactorizationTrivial(group)
            self.group = group;
        end

        function letters = preimageElement(self, g)
            letters = zeros(1, 0);
        end

    end

end
