classdef AbstractGroupLaws < replab.laws.FiniteGroupLaws
% Law checks for abstract groups: the operations below are pretty expensive

    methods

        function self = AbstractGroupLaws(T)
            self@replab.laws.FiniteGroupLaws(T);
        end

    end

    methods

        function law_relators_are_satisfied_by_permutation_realization_(self)
            assert(self.T.imagesDefineMorphism(self.T.permutationGroup, self.T.permutationGroup.generators));
        end

    end

end
