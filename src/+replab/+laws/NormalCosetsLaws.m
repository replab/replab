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
    end

end
