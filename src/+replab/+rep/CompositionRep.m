classdef CompositionRep < replab.Rep
% Composition of a representation with a morphism

    properties (SetAccess = protected)
        first % (`+replab.Morphism`): Morphism
        second % (`+replab.Rep`): Representation
    end

    methods

        function self = CompositionRep(first, second)
            assert(isa(first, 'replab.Morphism'));
            assert(isa(second, 'replab.Rep'));
            if second.knownUnitary
                args = {'isUnitary' true};
            else
                args = {};
            end
            % TODO: propagate more properties
            self@replab.Rep(first.source, second.field, second.dimension, args{:});
            self.first = first;
            self.second = second;
            self.group = first.source;
            self.field = second.field;
            self.dimension = second.dimension;
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function c = decomposeTerm(self)
            c = {self.second};
        end

        function r = composeTerm(self, newParts)
            r = replab.rep.Composition(self.first, newParts{1});
        end

        function rho = image_exact(self, g)
            rho = self.second.image(self.first.imageElement(g), 'exact');
        end

        function rho = image_double_sparse(self, g)
            rho = self.second.image(self.first.imageElement(g), 'double/sparse');
        end

    end

    methods % Implementations

        function b = isExact(self)
            b = self.second.isExact;
        end

    end

end
