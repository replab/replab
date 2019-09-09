classdef DomainLaws < replab.Laws
    properties (SetAccess = protected)
        T; % Domain under test
    end
    methods
        function self = DomainLaws(T)
            self.T = T;
        end
        function law_eqv_T(self, x)
            self.assert(self.T.eqv(x, x));
        end
    end
end
