classdef PermToGroup < replab.Morphism

    properties (SetAccess = protected)
        preimages % (cell(1,n) of `.source` elements): Preimages
        images % (cell(1,n) of `.target` elements): Images
        imageElementFun % (function_handle or ``[]``): Optional function_handle that implements `.imageElement`
    end

    methods

        function self = PermToGroup(source, target, preimages, images, imageElementFun)
        % Constructs a morphism from a permutation group to a finite group
        %
        % Args:
        %   source (`+replab.PermutationGroup`): Source of the morphism
        %   target (`+replab.Group`): Target of the morphism
        %   preimages (cell(1,n) of `.source` elements): Preimages
        %   images (cell(1,n) of `.target` elements): Images
            self.source = source;
            self.target = target;
            self.preimages = preimages;
            self.images = images;
            self.imageElementFun = imageElementFun;
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
            c = replab.bsgs.ChainWithImages.make(n, self.target, self.preimages, self.images, ...
                                                 self.source.chain.base, self.source.order);
        end

    end

    methods % Implementations

        function t = imageElement(self, s)
            if ~isempty(self.imageElementFun)
                f = self.imageElementFun;
                t = f(s);
                return
            end
            t = self.chain.image(s);
        end

    end

end
