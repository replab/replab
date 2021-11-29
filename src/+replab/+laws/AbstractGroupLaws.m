classdef AbstractGroupLaws < replab.laws.FiniteGroupLaws
% Law checks for abstract groups: the operations below are pretty expensive

    methods

        function self = AbstractGroupLaws(S)
            assert(false); % TODO
            self@replab.laws.FiniteGroupLaws(S);
        end

    end

    methods

        function law_relators_are_satisfied_by_permutation_realization_(self)
            assert(self.S.isMorphismByImages(self.S.permutationGroup, 'images', self.S.permutationGroup.generators));
        end

    end

end
