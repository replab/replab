classdef Coset < replab.Coset & replab.gen.FiniteSet

    methods % Implementations

        % gen.FiniteSet

        function l = compatibleWithNiceIsomorphism(self, iso)
            l = self.group.compatibleWithNiceIsomorphism(iso) && self.subgroup.compatibleWithNiceIsomorphism(iso);
        end

    end

end
