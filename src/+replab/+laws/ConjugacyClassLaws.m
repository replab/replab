classdef ConjugacyClassLaws < replab.Laws
% Law checks for a conjugacy class

    properties (SetAccess = protected)
        C % (`+replab.ConjugacyClass`): Conjugacy class
    end

    methods

        function self = ConjugacyClassLaws(C)
            self.C = C;
        end

        function law_representative_C(self, c)
            r = replab.ConjugacyClasses.representative(self.C.group, c);
            self.C.group.assertEqv(self.C.representative, r);
        end

        function law_contains_C(self, c)
            self.assert(self.C.contains(c));
        end

    end

end
