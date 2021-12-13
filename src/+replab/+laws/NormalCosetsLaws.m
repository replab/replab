classdef NormalCosetsLaws < replab.laws.LeftCosetLaws & replab.laws.RightCosetLaws

    methods

        function self = NormalCosetsLaws(N)
            self@replab.laws.LeftCosetLaws(N);
            self@replab.laws.RightCosetLaws(N);
        end

    end

    methods % Laws

        function law_subgroup_must_be_normal_(self)
            self.assert(self.S.subgroup.isNormalSubgroupOf(self.S));
        end
        
        function law_same_coset_G(self, g)
            assert(self.S.contains(self.S.group.compose(self.S.sample, g)));
            assert(self.S.contains(self.S.group.compose(g, self.S.sample)));
        end
        
    end

end
