classdef PermToFiniteGroup < replab.FiniteMorphism

    properties (SetAccess = protected)
        fastMorphism % (`+replab.fm.PermToGroup`): Fast morphism that cannot compute preimage representatives
    end

    methods

        function self = PermToFinite(source, target, images)
        % Constructs a morphism from a permutation group to a finite group
        %
        % Args:
        %   source (`+replab.PermutationGroup`): Source of the morphism
        %   target (`+replab.Group`): Target of the morphism
        %   images (cell(1,\*) of `.target` elements): Images of the source generators
            self.source = source;
            self.target = target;
            self.fastMorphism = replab.fm.PermToGroup(source, target, images);
        end

        function m = slowFiniteMorphism(self)
            m = self.cached('slowFiniteMorphism', @() self.computeSlowFiniteMorphism);

    end

    methods (Access = protected)

        function sfm = computeSlowFiniteMorphism(self)
            images = self.fastMorphism.imagesSourceGenerators;
            % we have the following diagram, that we use to compute the slow finite morphism
            %
            % self.source -- snm --> sng
            %     |                   |
            %    sfm                 ppm
            %     |                   |
            %     v                   v
            % self.target -- tnm --> tng
            snm = self.source.niceMorphism;
            sng = snm.image;
            %
            tnm = self.target.niceMorphism;
            prmImages = cellfun(@(g) tnm.imageElement(g), images, 'uniform', 0);
            prmSource = self.source.niceMorphism.image;
            prmTarget = self.target.niceMorphism.image;
            m = self.source.morphismByImages(prmImages
        end

    end

end

end
