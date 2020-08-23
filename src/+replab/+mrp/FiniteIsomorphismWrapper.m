classdef FiniteIsomorphismWrapper < replab.FiniteIsomorphism
% Wraps a finite morphism which happens to be injective, and restricts its target

    properties (SetAccess = protected)
        finiteMorphism % (`+replab.FiniteMorphism`): Wrapped finite morphism that happens to have trivial kernel
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

    end

    methods % Implementations

        % Morphism

        function t = imageElement(self, s)
            t = self.finiteMorphism.imageElement(s);
        end

        % Isomorphism

        function s = preimageElement(self, t)
            s = self.finiteMorphism.preimageRepresentative(t);
        end

        % FiniteMorphism

        function s = preimageRepresentative(self, t)
            s = self.finiteMorphism.preimageRepresentative(t);
        end

        function T = imageGroup(self, S)
            T = self.finiteMorphism.imageGroup(S);
        end

        function S = preimageGroup(self, T)
            S = self.finiteMorphism.preimageGroup(T);
        end

    end

end
