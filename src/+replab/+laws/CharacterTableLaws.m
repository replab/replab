classdef CharacterTableLaws < replab.Laws
% Law checks for character tables

    properties (SetAccess = protected)
        C % (`+replab.CharacterTable`): Character table
    end

    methods

        function self = CharacterTableLaws(C)
            self.C = C;
        end

        function laws = laws_classes(self)
            laws = replab.laws.Collection({self.C.classes.laws});
        end

        function law_irrep_(self)
            for i = 1:self.C.nIrreps
                irrep = self.C.irreps{i};
                if ~isempty(irrep)
                    for j = 1:length(irrep.preimages)
                        g = irrep.preimages{j};
                        ci = self.C.classes.classIndex(g);
                        self.assert(trace(irrep.images_internal{j}) == self.C.characters(i, ci));
                    end
                end
            end
        end
    end

end
