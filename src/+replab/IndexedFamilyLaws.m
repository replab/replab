classdef IndexedFamilyLaws < replab.Laws
    properties (SetAccess = protected)
        E;
        B;
        D;
    end
    methods
        function self = IndexedFamilyLaws(E)
            self.E = E;
            self.B = replab.domain.vpi(1, E.size);
            if E.size > 2^53 - 1
                maxD = 2^53 - 1;
            else
                maxD = double(E.size);
            end
            self.D = replab.domain.intAsDouble(1, maxD);
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
        function law_at_find_roundtrip_double_D(self, ind)
            a = self.E.at(ind);
            ind1 = double(self.E.find(a));
            self.D.assertEqv(ind, ind1);
        end
        function law_at_find_roundtrip_string_D(self, ind)
            a = self.E.at(num2str(ind));
            ind1 = double(self.E.find(a));
            self.D.assertEqv(ind, ind1);
        end
    end
end
