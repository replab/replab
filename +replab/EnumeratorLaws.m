classdef EnumeratorLaws < replab.Laws
    properties (SetAccess = protected)
        E;
        B;
    end
    methods
        function self = EnumeratorLaws(E)
            self.E = E;
            self.B = replab.domain.vpi(1, E.size);
        end
        function law_sample_is_contained_N1(self, dummy)
            e = self.E.sample;
            self.assert(~isempty(self.E.find(e)));
        end
        function law_at_find_roundtrip_B(self, ind)
            a = self.E.at(ind);
            ind1 = self.E.find(a);
            self.B.assertEqv(ind, ind1);
        end
    end
end
