classdef DefiningRep < replab.Rep
% Defining representation of a unitary or orthogonal group

    methods

        function self = DefiningRep(group)
            assert(isa(group, 'replab.UnitaryGroup') || isa(group, 'replab.OrthogonalGroup'));
            if isa(group, 'replab.UnitaryGroup')
                field = 'C';
            else % replab.OrthogonalGroup
                field = 'R';
            end
            self@replab.Rep(group, field, group.n, 'isUnitary', true, 'trivialDimension', 0, 'isIrreducible', true);
        end

    end

    methods (Access = protected) % Implementations

        function d = computeDouble(self)
            d = self;
        end

        function e = computeErrorBound(self)
            e = inf;
        end

        function rho = image_intval(self, g)
            rho = replab.rep.findEnclosingUnitary(g);
        end

        function rho = image_double_sparse(self, g)
            rho = g;
        end

    end

    methods % Implementations

        function b = isExact(self)
            b = false; % Note: redundant, as parent is already false
        end

    end

end
