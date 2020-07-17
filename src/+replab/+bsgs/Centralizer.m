classdef Centralizer < replab.bsgs.Backtrack
% Computes the centralizer of a subgroup

    properties
        other % `+replab.PermutationGroup`: Subgroup to centralize

        orbitDescr
        orbitReps
        otherOrbits
        otherTransversals

        cutAfter
    end

    methods

        function self = Centralizer(group, other, knownSubgroup, debug)
            if nargin < 4 || isempty(debug)
                debug = false;
            end
            if nargin < 3 || isempty(knownSubgroup)
                knownSubgroup = group.trivialSubgroup;
            end
            degree = group.domainSize;
            orbits = other.orbits.blocks;
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
            self@replab.bsgs.Backtrack(group, longBase, knownSubgroup, knownSubgroup, debug);
            base = self.base;
            cutAfter = find(self.orbitSizes ~= 1, 1, 'last');
            if ~isempty(cutAfter)
                base = base(1:cutAfter);
            end
            self.cutAfter = cutAfter;
            L = length(base);
            j = find(cellfun(@(orbit) any(orbit == base(end)), orbits), 1, 'last');
            relOrbits = orbits(1:min(j+1, length(orbits)));
            numRelOrbits = length(relOrbits);
            otherOrbits = cell(1, numRelOrbits);
            otherTransversals = cell(1, numRelOrbits);
            for j = 1:numRelOrbits
                rep = orbitReps(j);
                [O T] = replab.bsgs.orbitTransversal(degree, other.generatorsAsMatrix', rep);
                otherOrbits{j} = O;
                otherTransversals{j} = T;
            end
            self.other = other;
            self.orbitDescr = orbitDescr;
            self.orbitReps = orbitReps;
            self.otherOrbits = otherOrbits;
            self.otherTransversals = otherTransversals;
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

        function ok = test(self, l, gPrev, ul)
            if l > self.cutAfter
                ok = true;
                return
            end
            beta = self.base(l);
            if any(self.orbitReps == self.base(l))
                g = gPrev(ul);
                repOrbIndex = self.orbitDescr(beta);
                rep = self.orbitReps(repOrbIndex);
                im = gPrev(ul(beta));
                imRep = gPrev(ul(rep));
                trEl = self.otherTransversals{repOrbIndex}(:, find(self.otherOrbits{repOrbIndex} == beta))';
                ok = (im == trEl(imRep));
            else
                ok = true;
            end
        end

    end

end
