classdef DefiningRep < replab.Rep
% Defining representation of a unitary group

    methods

        function self = DefiningRep(group)
            assert(isa(group, 'replab.UnitaryGroup'));
            self.group = group;
            self.field = 'C';
            assert(group.n >= 1);
            self.dimension = group.n;
            self.isUnitary = true;
            self.isTrivial = false;
            self.isIrreducible = true;
        end

        function rho = image_internal(self, g)
            rho = g;
        end

        function rho = imageInverse_internal(self, g)
            rho = g';
        end

    end

end
