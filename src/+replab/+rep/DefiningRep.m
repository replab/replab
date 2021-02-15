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

        function e = computeErrorBound(self)
            e = inf;
        end

        function rho = image_double_sparse(self, g)
            rho = g;
        end

    end

    methods % Implementations

        % Rep

        function b = isExact(self)
            b = false; % Note: redundant, as parent is already false
        end

        function p = invariantBlocks(self)
            p = replab.Partition.fromBlocks({1:self.dimension});
        end

    end

end
