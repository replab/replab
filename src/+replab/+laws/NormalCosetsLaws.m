classdef NormalCosetsLaws < replab.laws.LeftCosetLaws & replab.laws.RightCosetLaws

    methods

        function self = NormalCosetsLaws(N)
            self@replab.laws.LeftCosetLaws(N);
            self@replab.laws.RightCosetLaws(N);
        end

    end

end
