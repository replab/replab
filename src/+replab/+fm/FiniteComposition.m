classdef FiniteComposition < replab.FiniteMorphism

    properties (SetAccess = protected)
        first % (`+replab.FiniteMorphism`): First morphism
        second % (`+replab.FiniteMorphism`): Second morphism
    end

    methods

        function self = FiniteComposition(second, first)
            self.target = second.target;
            self.source = first.source;
            self.first = first;
            self.second = second;
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
