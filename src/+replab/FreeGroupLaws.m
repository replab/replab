classdef FreeGroupLaws < replab.GroupLaws
% Group axioms

    methods

        function self = FreeGroupLaws(T)
            self@replab.GroupLaws(T);
        end

        function law_parse_roundtrip_T(self, x)
            self.T.assertEqv(x, self.T.parse(x.toString));
        end

    end

end
