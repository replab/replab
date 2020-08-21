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

    end

end
