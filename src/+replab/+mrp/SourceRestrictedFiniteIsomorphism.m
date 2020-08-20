classdef SourceRestrictedFiniteIsomorphism < replab.FiniteIsomorphism

    properties (SetAccess = protected)
        original % (`+replab.FiniteIsomorphism`): Original isomorphism
    end

    methods

        function self = SourceRestrictedFiniteIsomorphism(original, newSource)
            assert(isa(original, 'replab.FiniteIsomorphism'));
            assert(isa(newSource, 'replab.FiniteGroup') && newSource.isSubgroupOf(original.source));
            targetGenerators = cellfun(@(g) original.imageElement(g), newSource.generators, 'uniform', 0);
            newTarget = original.target.subgroupWithGenerators(targetGenerators);
            self.original = original;
            self.source = newSource;
            self.target = newTarget;
        end

    end

    methods % Implementations

        function t = imageElement(self, s)
            t = self.original.imageElement(s);
        end

        function s = preimageElement(self, t)
            s = self.original.preimageElement(t);
        end

        function T = imageGroup(self, S)
            T = self.original.imageGroup(S);
        end

        function S = preimageGroup(self, T)
            S = self.original.preimageGroup(T);
        end

        function s = preimageRepresentative(self, t)
            s = self.original.preimageRepresentative(t);
        end

        function S = preimagesElement(self, t)
            S = self.original.preimagesElement(t);
        end

        function r = restrictedSource(self, newSource)
            r = replab.mrp.SourceRestrictedFiniteIsomorphism(self.original, newSource);
        end

    end

end
