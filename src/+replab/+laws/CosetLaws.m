classdef CosetLaws < replab.laws.FiniteSetLaws

    properties (SetAccess = protected)
        G % (`+replab.FiniteGroup`): Subgroup
    end

    methods

        function self = CosetLaws(S)
            self@replab.laws.FiniteSetLaws(S);
            self.G = S.subgroup;
        end

    end

end
