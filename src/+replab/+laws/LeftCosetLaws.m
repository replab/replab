classdef LeftCosetLaws < replab.laws.CosetLaws

    methods

        function self = LeftCosetLaws(S)
            self@replab.laws.CosetLaws(S);
        end

    end

    methods

        function law_same_coset_G(self, g)
            assert(self.S.contains(self.S.group.compose(self.S.sample, g)));
        end

    end

end
