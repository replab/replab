classdef DoubleCosets < replab.DoubleCosets
% Describes the set of double cosets in a group
%
% A double coset is a set of the form ``{ H g K } = { h g k : h \in H, k \in K}`` for subgroups
% ``H`` and ``K`` of a group ``G``

    methods

        function self = DoubleCosets(group, leftSubgroup, rightSubgroup)
            assert(group.hasSameTypeAs(leftSubgroup));
            assert(group.hasSameTypeAs(rightSubgroup));
            self.group = group;
            self.leftSubgroup = leftSubgroup;
            self.rightSubgroup = rightSubgroup;
        end

    end

    methods

        function t = cosetRepresentative(self, g)
            type = self.group.type;
            S = replab.perm.Set(type.domainSize);
            rep = replab.bsgs.Cosets.leftRepresentative(self.rightSubgroup.lexChain, g);
            min = rep;
            S.insert(rep');
            toCheck = 1;
            % compute the orbit of the left cosets ``element rightSubgroup``
            while ~isempty(toCheck)
                i = toCheck(end);
                toCheck = toCheck(1:end-1);
                rep = S.at(i)';
                for j = 1:self.leftSubgroup.nGenerators
                    h = self.leftSubgroup.generator(j);
                    g = h(rep);
                    rep = replab.bsgs.Cosets.leftRepresentative(self.rightSubgroup.lexChain, g);
                    ind = S.find(rep');
                    if ind == 0
                        if type.compare(rep, min) < 0
                            min = rep;
                        end
                        ind = S.insert(rep');
                        toCheck = [toCheck ind];
                    end
                end
            end
            t = min;
        end

        function T = transversal(self)
            leftCosets = self.group / self.rightSubgroup;
            M = leftCosets.transversalAsMatrix;
            HActionOnLeftCosets = leftCosets.leftAction.imageGroup(self.leftSubgroup);
            orbits = HActionOnLeftCosets.orbits.blocks;
            Tperms = zeros(0, self.leftSubgroup.domainSize);
            for i = 1:length(orbits)
                orbit = orbits{i};
                sorted = sortrows(M(:,orbit)');
                Tperms(end+1, :) = sorted(1,:);
            end
            Tperms = sortrows(Tperms);
            T = arrayfun(@(i) Tperms(i,:), 1:size(Tperms, 1), 'uniform', 0);
        end

    end

end
