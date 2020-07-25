classdef AbstractGroupIsomorphismChain < replab.FiniteIsomorphism
% Describes an isomorphism from an abstract group to its realization (permutation group) using a stabilizer chain

    methods

        function self = AbstractGroupIsomorphismChain(source)
        % Constructs the nice isomorphism from an abstract finite group to its permutation group
        %
        % Args:
        %   source (`+replab.AbstractGroup`): Source of the isomorphism
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
            c = replab.bsgs.ChainWithWords(self.target);
            c.sgsWordQuick(1000);
            c.setCompleted;
        end

        function T = imageGroup(self, S)
            T = S.niceGroup;
        end

        function t = imageElement(self, s)
            t = self.source.niceImage(s);
        end

        function s = preimageElement(self, t)
            s = self.source.fromLetters(self.chain.word(t));
        end

        function S = preimageGroup(self, T)
            gens = cellfun(@(t) self.source.fromLetters(self.chain.word(t)), T.generators, 'uniform', 0);
            S = self.source.niceSubgroup(gens, T.order, T);
        end

    end

end
