classdef Centralizer
% Computes the centralizer of a subgroup

    properties
        group % `+replab.PermutationGroup`: Group
        other % `+replab.PermutationGroup`: Subgroup to centralize
    end

    methods

        function self = Centralizer(group, other)
            self.group = group;
            self.other = other;
        end


        function s = subgroup(self)
            if self.group.isTrivial || self.other.isTrivial
                s = self.group;
                return
            end
            degree = self.group.domainSize;
            orbits = self.other.orbits.blocks;
            numOrbits = length(orbits);
            [~, ind] = sort(-cellfun(@length, orbits));
            orbits = orbits(ind); % sort by decreasing length
            longBase = zeros(1, 0);
            orbitReps = zeros(1, numOrbits);
            orbitRepsIndices = zeros(1, numOrbits);
            orbitDescr = zeros(1, degree);
            for i = 1:numOrbits
                orbit = orbits{i};
                orbitReps(i) = orbit(1);
                orbitRepsIndices(i) = length(longBase) + 1;
                orbitDescr(orbit) = i;
                longBase = [longBase orbit];
            end
            chain = self.group.chain.mutableCopy;
            chain.baseChange(longBase);
            chain.removeRedundantBasePointsAtTheEnd;
            base = chain.base;
            L = length(base);
            j = find(cellfun(@(orbit) any(orbit == base(end)), orbits), 1, 'last');
            relOrbits = orbits(1:min(j+1, length(orbits)));
            numRelOrbits = length(relOrbits);
            otherOrbits = cell(1, numRelOrbits);
            otherTransversals = cell(1, numRelOrbits);
            for j = 1:numRelOrbits
                rep = orbitReps(j);
                [O T] = replab.bsgs.orbitTransversal(degree, self.other.generatorsAsMatrix', rep);
                otherOrbits{j} = O;
                otherTransversals{j} = T;
            end
            tests = cell(1, L);
            for l = 1:L
                if any(base(l) == orbitReps)
                    tests{l} = @(g, data) deal(true, []);
                else
                    tests{l} = @(g, data) deal(replab.bsgs.Centralizer.test(g, base(l), orbitDescr, orbitReps, otherOrbits, otherTransversals), []);
                end
            end

            startData = [];
            initSubgroup = [];
            bt = replab.bsgs.Backtrack(chain, @(g) self.prop(g), tests, startData, initSubgroup, true);
            s = replab.PermutationGroup.fromChain(bt.subgroup, self.group.type);
        end

        function b = prop(self, g)
            for i = 1:self.other.nGenerators
                gen = self.other.generator(i);
                if ~isequal(g(gen), gen(g))
                    b = false;
                    return
                end
            end
            b = true;
        end

    end

    methods (Static)

        function [ok, data] = test(g, beta, orbitDescr, orbitReps, otherOrbits, otherTransversals)
            repOrbIndex = orbitDescr(beta);
            rep = orbitReps(repOrbIndex);
            im = g(beta);
            imRep = g(rep);
            trEl = otherTransversals{repOrbIndex}(:, find(otherOrbits{repOrbIndex} == beta))';
            ok = (im == trEl(imRep));
        end

    end



end
