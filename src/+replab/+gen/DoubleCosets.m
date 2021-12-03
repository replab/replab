classdef DoubleCosets < replab.DoubleCosets

    methods

        function self = DoubleCosets(nice, niceIsomorphism)
            self.group = niceIsomorphism.preimageElement
            assert(group.hasSameTypeAs(H));
            assert(group.hasSameTypeAs(K));
            self.isomorphism = group.niceMorphism;
            self.group = group;
            self.H = H;
            self.K = K;
            self.Hprmgrp = group.niceMorphism.imageGroup(H);
            self.Kprmgrp = group.niceMorphism.imageGroup(K);
            self.leftCosets = group / K;
        end

        function s = size(self)
        % Returns the number of double cosets
        %
        % Returns:
        %   integer: Number of double cosets
            s = length(self.transversal);
        end

        function t = cosetRepresentative(self, g)
        % Returns the canonical coset representative corresponding to the given element
        %
        % Args:
        %   g (element of `.group`): Group element
        %
        % Returns:
        %   t (element of `.group`): Double coset canonical representative
            S = replab.perm.Set(self.Hprmgrp.domainSize);
            elP = self.isomorphism.imageElement(g);
            repP = replab.bsgs.Cosets.leftRepresentative(self.Kprmgrp.lexChain, elP);
            minP = repP;
            S.insert(repP');
            toCheck = 1;
            % compute the orbit of the left cosets ``element K``
            while ~isempty(toCheck)
                i = toCheck(end);
                toCheck = toCheck(1:end-1);
                rep = S.at(i)';
                for j = 1:self.Hprmgrp.nGenerators
                    h = self.Hprmgrp.generator(j);
                    elP = h(rep);
                    repP = replab.bsgs.Cosets.leftRepresentative(self.Kprmgrp.lexChain, elP);
                    ind = S.find(repP');
                    if ind == 0
                        if replab.DoubleCoset.lexCompare(repP, minP) < 0
                            minP = repP;
                        end
                        ind = S.insert(repP');
                        toCheck = [toCheck ind];
                    end
                end
            end
            t = self.isomorphism.preimageElement(minP);
        end

        function T = transversal(self)
        % Returns all the canonical representatives of cosets
        %
        % Returns:
        %   cell(1, \*) of `.group` elements: Transversal
            T = self.cached('transversal', @() self.computeTransversal);
        end

        function T = computeTransversal(self)
        % See `.transversal`
            M = self.leftCosets.transversalAsMatrix;
            HActionOnLeftCosets = self.leftCosets.leftAction.imageGroup(self.H);
            orbits = HActionOnLeftCosets.orbits.blocks;
            Tperms = zeros(0, self.Hprmgrp.domainSize);
            for i = 1:length(orbits)
                orbit = orbits{i};
                sorted = sortrows(M(:,orbit)');
                Tperms(end+1, :) = sorted(1,:);
            end
            Tperms = sortrows(Tperms);
            T = arrayfun(@(i) self.isomorphism.preimageElement(Tperms(i,:)), 1:size(Tperms, 1), 'uniform', 0);
        end

        function C = elements(self)
        % Returns the set of double cosets as a cell array
        %
        % Returns:
        %   cell(1,\*) of `+replab.DoubleCoset`: Set of double cosets
            T = self.transversal;
            C = cellfun(@(t) replab.DoubleCoset(self.H, t, self.K, self.group), T, 'uniform', 0);
        end

    end

end
