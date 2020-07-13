classdef FiniteMorphism < replab.Obj
% Describes a morphism between finite groups

    properties (SetAccess = protected)
        source % (`.FiniteGroup`): Source group
        target % (`.FiniteGroup`): Target group
    end

    methods

        function K = kernel(self)
        % Returns the kernel of this morphism
        %
        % Returns:
        %   `+replab.FiniteGroup`: Maximal subgroup of `.source` with trivial image
            K = self.cached('kernel', @() self.computeKernel);
        end

        function K = computeKernel(self)
            K = self.preimageGroup(self.target.trivialSubgroup);
        end

        function s = preimagesElement(self, t)
            error('Abstract');
        end

        function S = preimageGroup(self, T)
            error('Abstract');
        end

        function I = image(self)
        % Returns the image of this morphism
            I = self.cached('image', @() self.imageGroup(self.source));
        end

        function t = imageElement(self, s)
        % Returns the image of the given source element
        %
        % Args:
        %   s (element of `.source`): Element to compute the image of
            error('Abstract');
        end

        function T = imageGroup(self, S)
        % Computes the image of a group
        %
        % Args:
        %   S (subgroup of `.source`): Group to compute the image of
            images = cellfun(@(g) self.imageElement(g), S.generators, 'uniform', 0);
            T = self.target.subgroup(images);
        end

        function res = compose(self, applyFirst)
        % Composition of morphisms, the right hand side applied first
        %
        % Note: if both morphisms are isomorphisms, the resulting morphism is also an isomorphism.
        %
        % Args:
        %   applyFirst (`.FiniteMorphism`): Morphism to apply first
        %
        % Returns:
        %   `.Morphism`: The composition of the given morphism applied first, followed by this morphism.
            m = replab.fm.Composition(self, applyFirst);
        end

        function res = andThen(self, applyLast)
        % Composition of morphisms, the right hand side applied last
        %
        % Args:
        %   applyLast (`.FiniteMorphism`): Morphism to apply last
        %
        % Returns:
        %   `.FiniteMorphism`: The composition of this morphism applied first, followed by the given morphism
            res = applyLast.compose(self);
        end

        function res = mtimes(self, applyFirst)
        % Shorthand for `.compose`
            res = self.compose(applyFirst);
        end

    end

end
