classdef DomainLaws < replab.Laws
% Law checks for a generic Domain

    properties (SetAccess = protected)
        S % (`+replab.Domain`): Domain under test
    end

    methods

        function self = DomainLaws(S)
            self.S = S;
        end

        function law_eqv_S(self, x)
        % Checks that an element is equivalent to itself
            self.assert(self.S.eqv(x, x));
        end

    end

end
