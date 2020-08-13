classdef FiniteIsomorphismWrapper < replab.FiniteIsomorphism
% Wraps a finite morphism which happens to be injective, and restricts its target

    properties (SetAccess = true)
        finiteMorphism
    end

    methods

        function self = FiniteIsomorphismWrapper(finiteMorphism)
            assert(finiteMorphism.kernel.isTrivial);
            self.source = finiteMorphism.source;
            self.target = finiteMorphism.image;
            self.finiteMorphism = finiteMorphism;
        end

    end

    methods (Access = protected)

    % computeImage/computeKernel/computeInverse in FiniteIsomorphism
        function I = computeImageSourceGenerators(self)
            I = self.finiteMorphism.imageSourceGenerators;
        end

    end

    methods % Implementations

        % Morphism

        function t = imageElement(self, s)
            t = self.finiteMorphism.imageElement(s);
        end

        % Isomorphism

        function s = preimageElement(t)
            s = self.finiteMorphism.preimageRepresentative(self, t);
        end

        % FiniteMorphism

        function s = preimageRepresentative(self, t)
            s = self.finiteMorphism.preimageRepresentative(t);
        end

        function S = preimageGroup(self, T)
            S = self.finiteMorphism.preimageGroup(T);
        end

        % imageGroup
        function T = imageGroup(self, S)
            T = self.finiteMorphism.imageGroup(S);
        end

    end

end
