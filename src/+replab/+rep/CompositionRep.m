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

        function b = hasTorusImage(self)
            b = self.group.hasReconstruction && (~isempty(self.first.torusMap) || self.group.maximalTorusDimension == 0);
            b = b && self.second.hasTorusImage;
        end

        function [torusMap, torusInjection, torusProjection] = torusImage(self)
            if isempty(self.first.torusMap)
                n1 = self.group.maximalTorusDimension;
                assert(n1 == 0);
                n2 = self.second.group.maximalTorusDimension;
                torusMap1 = zeros(n2, n1);
            else
                torusMap1 = self.first.torusMap;
            end
            [torusMap2, torusInjection, torusProjection] = self.second.torusImage;
            torusMap = torusMap2 * torusMap1;
        end

    end

end
