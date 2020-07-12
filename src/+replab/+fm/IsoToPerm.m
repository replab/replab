classdef IsoToPerm < replab.FiniteIsomorphism

    properties
        imageFun % (function_handle): Image function
    end

    methods

        function self = IsoToPerm(source, imageFun)
            self.source = source;
            identity = imageFun(source.identity)
            self.target = replab.SymmetricGroup(length(identity));
            self.imageFun = imageFun;
        end

        function c = chain(self)
        % Returns the BSGS chain that enables the computation of preimages
        %
        % Returns:
        %   `+replab.+bsgs.ChainWithImages`: The BSGS chain
            c = self.cached('chain', @() self.computeChain);
        end

        function c = computeChain(self)
            imageGenerators = cellfun(@(s) self.imageElement(s), self.source.generators, 'uniform', 0);
            n = self.target.domainSize;
            I = self.image;
            % TODO: optimize ChainWithImages by using deterministic Schreier-Sims while comparing orbits
            c = replab.bsgs.ChainWithImages.make(n, self.source, imageGenerators, source.generators, [], ...
                                                 I.base, I.order);
        end

        function t = imageElement(self, s)
            t = self.imageFun(s);
        end

        function s = preimageElement(self, t)
            s = self.chain.image(t);
        end

    end

end
