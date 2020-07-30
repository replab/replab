classdef PermToFinite < replab.FiniteMorphism

    properties (Access = protected)
        imageSourceGenerators_ % (cell(1,\*) of `.target` elements): Images of the source generators
    end

    methods

        function self = PermToFinite(source, target, images)
        % Constructs a morphism from a permutation group to a finite group
        %
        % Args:
        %   source (`+replab.PermutationGroup`): Source of the morphism
        %   target (`+replab.FiniteGroup`): Target of the morphism
        %   images (cell(1,\*) of `.target` elements): Images of the source generators
            self.source = source;
            self.target = target;
            self.imageSourceGenerators_ = images;
        end

        function c = chain(self)
        % Returns the BSGS chain that enables the computation of images
        %
        % Returns:
        %   `+replab.+bsgs.ChainWithImages`: The BSGS chain
            c = self.cached('chain', @() self.computeChain);
        end

    end

    methods (Access = protected)

        function c = computeChain(self)
            n = self.source.domainSize;
            % TODO: optimize ChainWithImages by using deterministic Schreier-Sims while comparing orbits
            c = replab.bsgs.ChainWithImages.make(n, self.target, self.source.generators, self.imageSourceGenerators_, [], ...
                                                 self.source.chain.base, self.source.order);
        end

    end

    methods % Implementations

        function t = imageElement(self, s)
            t = self.chain.image(s);
        end

        function I = imageSourceGenerators(self)
            I = self.imageSourceGenerators_;
        end

    end

end
