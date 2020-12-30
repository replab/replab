classdef RestrictionRep < replab.Rep
% Composition of a representation with a morphism

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Representation for the parent group, such that ``self.group`` is a subgroup of ``parent.group``
    end

    methods

        function self = RestrictionRep(parent, group)
            assert(isa(parent, 'replab.Rep'));
            assert(isa(group, 'replab.CompactGroup'));
            if second.knownUnitary
                args = {'isUnitary' true};
            else
                args = {};
            end
            % TODO: propagate more properties
            self@replab.Rep(group, parent.field, parent.dimension, args{:});
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function res = double(self)
            res = replab.rep.RestrictionRep(double(self.parent), self.group);
        end

        function c = decomposeTerm(self)
            c = {self.parent};
        end

        function r = composeTerm(self, newParts)
            r = replab.rep.RestrictionRep(newParts{1}, self.group);
        end

        function rho = image_exact(self, g)
            rho = self.parent.image(g, 'exact');
        end

        function rho = image_double_sparse(self, g)
            rho = self.parent.image(g, 'double/sparse');
        end

    end

    methods % Implementations

        function b = isExact(self)
            b = self.parent.isExact;
        end

    end

end
