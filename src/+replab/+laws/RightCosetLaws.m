classdef RightCosetLaws < replab.laws.CosetLaws

    methods

        function self = RightCosetLaws(S)
            self@replab.laws.CosetLaws(S);
        end

    end

    methods

        function law_same_coset_G(self, g)
            assert(self.S.contains(self.S.group.compose(g, self.S.sample)));
        end

    end

end
