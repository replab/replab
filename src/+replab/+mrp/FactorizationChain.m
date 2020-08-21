classdef FactorizationChain < replab.mrp.Factorization
% Factorizes permutations using Minkwitz algorithm based on a BSGS chain with short transversal words

    methods

        function self = FactorizationChain(group, chain)
        % Constructs a FactorizationChain
        %
        % Args:
        %   group (`+replab.PermutationGroup`): Group whose elements are factorized
        %   chain (`+replab.+bsgs.ChainWithWords`, optional): Computed stabilizer chain with reduced words
            self.group = group;
            if nargin >= 2 && ~isempty(chain)
                self.cache('chain', chain, 'ignore');
            end
        end

        function c = chain(self)
        % Returns the BSGS chain that enables the computation of preimages
        %
        % Returns:
        %   `+replab.+bsgs.ChainWithWords`: The BSGS chain
            c = self.cached('chain', @() self.computeChain);
        end

        function c = computeChain(self)
            c = replab.bsgs.ChainWithWords(self.group);
            c.sgsWordQuick(1000);
            c.setCompleted;
        end

        function letters = preimageElement(self, g)
            letters = self.chain.word(g);
        end

    end

end
