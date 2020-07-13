classdef NiceFiniteGroupIsomorphism < replab.FiniteIsomorphism

    methods

        function self = NiceFiniteGroupIsomorphism(source)
        % Constructs the nice isomorphism from a nice finite group to its permutation group
        %
        % Args:
        %   source (`+replab.NiceFiniteGroup`): Source of the isomorphism
        %
        % Returns:
        %   `+replab.FiniteIsomorphism`: The constructed isomorphism
            self.source = source;
            domainSize = length(source.niceImage(source.identity));
            self.target = replab.SymmetricGroup(domainSize);
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
            c = replab.bsgs.ChainWithImages.make(n, self.source, imageGenerators, self.source.generators, [], ...
                                                 I.chain.base, I.order);
        end

        function T = imageGroup(self, S)
            T = S.niceGroup;
        end

        function t = imageElement(self, s)
            t = self.source.niceImage(s);
        end

        function s = preimageElement(self, t)
            s = self.chain.image(t);
        end

        function S = preimageGroup(self, T)
            gens = cellfun(@(t) self.chain.image(t), T.generators, 'uniform', 0);
            S = self.source.niceSubgroup(gens, T.order, T);
        end

    end

end
