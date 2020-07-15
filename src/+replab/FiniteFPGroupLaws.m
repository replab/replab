classdef FiniteFPGroupLaws < replab.FiniteGroupLaws

    methods

        function self = FiniteFPGroupLaws(T)
            self@replab.FiniteGroupLaws(T);
        end

        function law_parse_presentation_roundtrip_(self)
            P1 = self.T.presentation;
            fpg = replab.FiniteFPGroup.parsePresentation(P1);
            P2 = fpg.presentation;
            assertEqual(P1, P2);
        end

    end

end
