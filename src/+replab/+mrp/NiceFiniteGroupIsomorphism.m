classdef NiceFiniteGroupIsomorphism < replab.FiniteIsomorphism

    methods

        function self = NiceFiniteGroupIsomorphism(source, target)
        % Constructs the nice isomorphism from a nice finite group to its permutation group
        %
        % Args:
        %   source (`+replab.NiceFiniteGroup`): Source of the isomorphism
        %   target (`+replab.PermutationGroup`): Target of the isomorphism
        %
        % Returns:
        %   `+replab.FiniteIsomorphism`: The constructed isomorphism
            self.source = source;
            self.target = target;
        end

        function c = chain(self)
        % Returns the BSGS chain that enables the computation of preimages
        %
        % Returns:
        %   `+replab.+bsgs.ChainWithImages`: The BSGS chain
            c = self.cached('chain', @() self.computeChain);
        end

        function c = computeChain(self)
            niceGroup = self.source.niceGroup;
            targetGenerators = cellfun(@(g) self.source.niceImage(g), self.source.generators, 'uniform', 0);
            c = replab.bsgs.ChainWithImages.make(self.target.domainSize, self.source, ...
                                                 targetGenerators, self.source.generators, [], ...
                                                 niceGroup.chain.base, niceGroup.order);
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
