classdef ConjugacyClassesLaws < replab.Laws
% Law checks for conjugacy classes

    properties (SetAccess = protected)
        G % (`+replab.FiniteGroup`): Group
        C % (`+replab.ConjugacyClasses`): Classes to test
    end

    methods

        function self = ConjugacyClassesLaws(C)
            self.C = C;
            self.G = C.group;
        end

        function law_classes_are_sorted_(self)
            I = self.G.type.sort(self.C.classRepresentatives);
            self.assert(isequal(I, 1:length(I)));
        end

        function law_classIndexOf_G(self, g)
            ind = self.C.classIndexOf(g);
            c = self.C.classes{ind};
            self.assert(c.contains(g));
        end

        function law_groupOrder_(self)
            order = sum([self.C.classSizes{:}]);
            self.assert(order == self.G.order);
        end

        function laws = laws_classes(self)
            laws = replab.laws.Collection(cellfun(@(c) c.laws, self.C.classes, 'uniform', 0));
        end

    end

end
