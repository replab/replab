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
            % TODO: propagate more properties
            self@replab.Rep(first.source, second.field, second.dimension, 'isUnitary', second.isUnitary);
            self.first = first;
            self.second = second;
            self.group = first.source;
            self.field = second.field;
            self.dimension = second.dimension;
        end

    end

    methods % Simplification rules

        function res = rewriteTerm_CompositionRepOfCompositionRep(self, options)
            if isa(self.second, 'replab.Rep.CompositionRep')
                newMorphism = self.first.andThen(self.second.first);
                newRep = self.second.second;
                res = replab.rep.CompositionRep(newMorphish, newRep);
            else
                res = [];
            end
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function b = computeIsUnitary(self)
            if isa(self.first, 'replab.Isomorphism')
                b = self.parent.isUnitary;
            else
                b = computeIsUnitary@replab.Rep(self);
            end
        end

        function e = computeErrorBound(self)
            e = self.second.errorBound;
        end

        function c = decomposeTerm(self)
            c = {self.second};
        end

        function r = composeTerm(self, newParts)
            r = replab.rep.CompositionRep(self.first, newParts{1});
        end

        function rho = image_exact(self, g)
            rho = self.second.image(self.first.imageElement(g), 'exact');
        end

        function rho = image_double_sparse(self, g)
            rho = self.second.image(self.first.imageElement(g), 'double/sparse');
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function M = matrixRowAction_double_sparse(self, g, M)
            h = self.first.imageElement(g);
            M = self.second.matrixRowAction(h, M);
        end

        function M = matrixColAction_double_sparse(self, g, M)
            h = self.first.imageElement(g);
            M = self.second.matrixColAction(h, M);
        end

    end

    methods % Implementations

        % Rep

        function b = isExact(self)
            b = self.second.isExact;
        end

    end

end
