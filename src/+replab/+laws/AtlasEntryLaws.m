classdef AtlasEntryLaws < replab.Laws
% Law checks for character tables

    properties (SetAccess = protected)
        A % (`+replab.AtlasEntry`): Atlas entry
    end

    methods

        function self = AtlasEntryLaws(A)
            self.A = A;
        end

        function laws = laws_characterTable(self)
            laws = replab.laws.Collection({self.A.characterTable.laws});
        end

        function laws = laws_group(self)
            laws = replab.laws.Collection({self.A.group.laws});
        end

    end

end
