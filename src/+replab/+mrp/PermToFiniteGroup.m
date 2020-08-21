classdef PermToFiniteGroup < replab.FiniteMorphism
% Morphism from a permutation group to a finite group
%
% Uses a BSGS chain with associated images to compute images of elements, delegating the operation to `+replab.+mrp.PermToGroup`, stored
% in `.fastMorphism`; the other operations are computed by a round-trip to permutation groups through the nice morphisms.

    properties (SetAccess = protected)
        fastMorphism % (`+replab.+mrp.PermToGroup`): Fast morphism that only computes images of elements
    end

    methods

        function self = PermToFiniteGroup(source, target, preimages, images)
        % Constructs a morphism from a permutation group to a finite group
        %
        % Args:
        %   source (`+replab.PermutationGroup`): Source of the morphism
        %   target (`+replab.Group`): Target of the morphism
        %   preimages (cell(1,n) of `.source` elements): Preimages
        %   images (cell(1,n) of `.target` elements): Images
            self.source = source;
            self.target = target;
            self.fastMorphism = replab.mrp.PermToGroup(source, target, preimages, images);
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
            images = self.fastMorphism.images;
            I = self.target.subgroup(images);
        end

    end

    methods (Access = protected)

        function sfm = computeSlowFiniteMorphism(self)
            preimages1 = self.fastMorphism.preimages;
            preimages2 = self.fastMorphism.images;
            iso2 = self.target.niceMorphism;
            images12 = cellfun(@(g) iso2.imageElement(g), preimages2, 'uniform', 0);
            inter = iso2.target.subgroup(images12);
            f1 = self.source.morphismByImages(inter, 'preimages', preimages1, 'images', images12);
            I = self.image;
            f2 = iso2.restrictedSource(I);
            % source -- f1 --> inter <-- f2 -- target
            sfm = f1.andThen(f2.inverse);
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
