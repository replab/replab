classdef AbstractGroupIsomorphism < replab.FiniteIsomorphism

    methods

        function self = AbstractGroupIsomorphism(source)
        % Constructs the nice isomorphism from an abstract finite group to its permutation group
        %
        % Args:
        %   source (`+replab.AbstractGroup`): Source of the isomorphism
        %   target (`+replab.PermutationGroup`): Target of the isomorphism
        %
        % Returns:
        %   `+replab.FiniteIsomorphism`: The constructed isomorphism
            self.source = source;
            self.target = source.permutationGroup;
        end

        function c = chain(self)
        % Returns the BSGS chain that enables the computation of preimages
        %
        % Returns:
        %   `+replab.+bsgs.ChainWithWords`: The BSGS chain
            c = self.cached('chain', @() self.computeChain);
        end

        function c = computeChain(self)
            target = self.target;
            c = replab.bsgs.ChainWithWords.make(target.domainSize, target.generators, target.chain.base, target.chain.order);
        end

        function T = imageGroup(self, S)
            T = S.niceGroup;
        end

        function t = imageElement(self, s)
            t = self.source.niceImage(s);
        end

        function s = preimageElement(self, t)
            s = self.source.fromLetters(self.chain.image(t));
        end

        function S = preimageGroup(self, T)
            gens = cellfun(@(t) self.source.fromLetters(self.chain.image(t)), T.generators, 'uniform', 0);
            S = self.source.niceSubgroup(gens, T.order, T);
        end

    end

end
