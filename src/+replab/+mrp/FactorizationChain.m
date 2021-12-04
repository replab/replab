classdef FactorizationChain < replab.mrp.Factorization
% Factorizes permutations using Minkwitz algorithm based on a BSGS chain with short transversal words

    methods

        function self = FactorizationChain(group, generators, useInverses, chain)
        % Constructs a FactorizationChain
        %
        % Args:
        %   group (`+replab.PermutationGroup`): Group whose elements are factorized
        %   generators (cell(1,\*) of elements of ``group``): Generators on which to perform the factorization
        %   useInverses (logical, optional): Whether to use inverses in the decomposition (default: true)
        %   chain (`+replab.+bsgs.ChainWithWords`, optional): Computed stabilizer chain with reduced words
            if nargin < 3 || isempty(useInverses)
                useInverses = true;
            end
            if nargin < 4
                chain = [];
            end
            self.group = group;
            self.generators = generators;
            self.useInverses = useInverses;
            if ~isempty(chain)
                self.cache('chain', chain, 'ignore');
            end
        end

        function ind = generatorInverseIndices(self)
        % Returns the inverse index of each generator when it exists
        %
        % See `+replab.+mrp.inverseIndices`
        %
        % Returns:
        %   integer(1,\*): Generator inverse indices
            ind = self.cached('generatorInverseIndices', @() replab.mrp.inverseIndices(self.group, self.generators));
        end

        function eo = generatorElementOrders(self)
        % Returns the element order of each generator
        %
        % Returns:
        %   integer(1,\*): Element order of each generator
            eo = self.cached('generatorElementOrders', @() cellfun(@(g) self.group.elementOrder(g), self.generators));
        end

        function c = chain(self)
        % Returns the BSGS chain that enables the computation of preimages
        %
        % Returns:
        %   `+replab.+bsgs.ChainWithWords`: The BSGS chain
            c = self.cached('chain', @() self.computeChain);
        end

        function c = computeChain(self)
            c = replab.bsgs.ChainWithWords(self.group, self.generators);
            c.sgsWordQuick(1000);
            c.setCompleted;
        end

        function letters = substituteInverses(self, letters)
            ind = find(letters < 0);
            if ~isempty(ind)
                II = self.generatorInverseIndices;
                inv = II(-letters(ind));
                mask = inv > 0;
                letters(ind(mask)) = inv(mask);
                if any(~mask)
                    EO = self.generatorElementOrders;
                    for i = fliplr(find(letters < 0)) % go right to left, as we are increasing the word length
                        l = -letters(i);
                        letters = [letters(1:i-1) repmat(l, 1, EO(l)-1) letters(i+1:end)];
                    end
                end
            end
        end

        function letters = factorize(self, g)
            letters = self.chain.word(g);
            if ~self.useInverses
                letters = self.substituteInverses(letters);
            end
        end

        function [l, r] = factorizeRepresentativeOfLeftCoset(self, leftCoset)
        % Returns a tentatively short word corresponding to an element of the given coset
        %
        % Args:
        %   leftCoset (`+replab.LeftCoset`): Left coset subset of `.group`
        %
        % Returns
        % -------
        %   l: integer(1,\*)
        %     Word expressed in letters
        %   r: group element
        %     Chosen coset representative
            [l, r] = self.chain.wordLeftCoset(leftCoset.representative, leftCoset.subgroup.chain);
        end

        function n = maximumWordLength(self)
            c = self.chain;
            n = 0;
            for i = 1:c.k
                W = c.nuw{i};
                l = cellfun(@(w) length(self.substituteInverses(w)), c.nuw{i});
                n = n + max(l);
            end
        end

    end

end
