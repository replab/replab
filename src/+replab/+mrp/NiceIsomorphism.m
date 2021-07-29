classdef NiceIsomorphism < replab.FiniteIsomorphism

    methods % Implementations

        function T = imageGroup(self, S)
            T = S.niceGroup;
        end

        function t = imageElement(self, s)
            t = self.source.niceImage(s);
        end

        function S = preimageGroup(self, T)
            gens = cellfun(@(t) self.preimageElement(t), T.generators, 'uniform', 0);
            S = self.source.niceSubgroup(gens, T.order, T);
        end

        function m = restrictedSource(self, newSource)
            m = newSource.niceMorphism;
        end

    end

end
