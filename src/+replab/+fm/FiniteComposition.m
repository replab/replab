classdef FiniteComposition < replab.FiniteMorphism & replab.fm.Composition

    methods

        function self = FiniteComposition(second, first)
            self@replab.fm.Composition(second, first);
        end

        function s = preimageRepresentative(self, t)
            s = self.first.preimageRepresentative(self.second.preimageRepresentative(t));
        end

        function S = preimageGroup(self, T)
            S = self.first.preimageGroup(self.second.preimageGroup(T));
        end

        function t = imageElement(self, s)
            t = self.second.imageElement(self.first.imageElement(s));
        end

        function T = imageGroup(self, S)
            T = self.second.imageGroup(self.first.imageGroup(S));
        end

    end

    methods (Access = protected)

        function I = computeImage(self)
            I = self.second.imageGroup(self.first.image);
        end

        function K = computeKernel(self)
            K = self.first.preimageGroup(self.second.kernel);
        end

    end

end
