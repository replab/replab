classdef Coset < replab.FiniteSet

    properties (SetAccess = protected)
        group % (`.FiniteGroup`): Group containing this coset
        subgroup % (`.FiniteGroup`): Subgroup used to decompose `.group`
    end

    methods % Implementations

        % Domain

        function l = laws(self)
            l = replab.laws.CosetLaws(self);
        end

        function b = eqv(self, lhs, rhs)
            b = self.type.eqv(lhs, rhs);
        end

        % FiniteSet

        function s = nElements(self)
            s = self.subgroup.order;
        end

    end

end
