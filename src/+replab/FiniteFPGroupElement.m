classdef FiniteFPGroupElement < replab.Str

    properties (SetAccess = protected)
        group % (`.FiniteFPGroup`): Finite finitely presented group this is an element of
        representative % (`.FreeGroupWord`): Word in the free group defining the equivalence class in the FP group
    end

    methods

        function self = FiniteFPGroupElement(group, representative)
            self.group = group;
            self.representative = representative;
        end

    end

    methods % Implementations

        function s = headerStr(self)
            s = self.toString;
        end

        function s = shortStr(self, maxColumns)
            s = self.toString;
        end

        function s = longStr(self, maxRows, maxColumns)
            s = {self.toString};
        end

    end

    methods

        function str = toString(self)
            str = ['[ ' self.representative.toString ' ]'];
        end

    end

end
