classdef IrreducibleLaws < replab.Laws

    properties
        irreducible; % Irreducible decomposition
        M;
    end

    methods

        function self = IrreducibleLaws(irreducible)
            self.irreducible = irreducible;
            d = self.irreducible.dimension;
            self.M = replab.domain.Matrices(irreducible.parent.field, d, d);
        end

        function law_decomposes_entire_space(self)
            U = self.irreducible.U;
            self.M.assertEqv(U' * U, eye(self.irreducible.parent.dimension));
        end

        function isotypicLaws = laws_isotypic_components(self)
            children = cellfun(@(x) replab.IsotypicLaws(x), self.irreducible.components, 'uniform', 0);
            isotypicLaws = replab.LawsCollection(children);
        end

    end
end
