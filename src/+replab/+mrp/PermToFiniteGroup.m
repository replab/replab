classdef PermToFiniteGroup < replab.FiniteMorphism
% Morphism from a permutation group to a finite group
%
% Uses a BSGS chain with associated images to compute images of elements, delegating the operation to `+replab.+mrp.PermToGroup`, stored
% in `.fastMorphism`; the other operations are computed by a round-trip to permutation groups through the nice morphisms.

    properties (SetAccess = protected)
        fastMorphism % (`+replab.+mrp.PermToGroup`): Fast morphism that only computes images of elements
    end

    methods

        function self = PermToFiniteGroup(source, target, images)
        % Constructs a morphism from a permutation group to a finite group
        %
        % Args:
        %   source (`+replab.PermutationGroup`): Source of the morphism
        %   target (`+replab.Group`): Target of the morphism
        %   images (cell(1,\*) of `.target` elements): Images of the source generators
            self.source = source;
            self.target = target;
            self.fastMorphism = replab.mrp.PermToGroup(source, target, images);
        end

        function m = slowFiniteMorphism(self)
        % Returns the finite morphism that implements the complete set of operations
        %
        % Returns:
        %   `+replab.FiniteMorphism`: The finite morphism equivalent to `.fastMorphism`, but with additional methods
            m = self.cached('slowFiniteMorphism', @() self.computeSlowFiniteMorphism);
        end

    end

    methods (Access = protected) % Implementations

        % FiniteMorphism

        function K = computeKernel(self)
            K = self.slowFiniteMorphism.kernel;
        end


        function I = computeImage(self)
            I = self.slowFiniteMorphism.image;
        end

        function I = computeImageSourceGenerators(self)
            I = self.fastMorphism.imageSourceGenerators;
        end

    end

    methods (Access = protected)

        function sfm = computeSlowFiniteMorphism(self)
            images = self.fastMorphism.imageSourceGenerators;
            tnm = self.target.niceMorphism;
            tng = tnm.target;
            prmImages = cellfun(@(g) tnm.imageElement(g), images, 'uniform', 0);
            pm = self.source.morphismByImages(tng, prmImages);
            sfm = pm.andThen(tnm.inverse);
        end

    end

    methods % Implementations

        % FiniteMorphism

        function s = preimageRepresentative(self, t)
            s = self.slowFiniteMorphism.preimageRepresentative(t);
        end

        function t = imageElement(self, s)
            t = self.fastMorphism.imageElement(s);
        end

    end

end
