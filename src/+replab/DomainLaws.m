classdef DomainLaws < replab.Laws
% Law checks for a generic Domain

    properties (SetAccess = protected)
        T % (`+replab.Domain`): Domain under test
    end

    methods

        function self = DomainLaws(T)
            self.T = T;
        end

        function law_eqv_T(self, x)
        % Checks that an element is equivalent to itself
            self.assert(self.T.eqv(x, x));
        end

    end

end
