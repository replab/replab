classdef PermutationGroupRightCosets < replab.Str
% Enables computations with the right cosets of a permutation group
%
% The transversal elements are given by the minimal elements of the coset under lexicographic ordering;
% the right cosets themselves are ordered by their transversal elements, in lexicographic ordering too.
%
% Thus, computations involving cosets are deterministic, and do not depend on the details of
% the construction of the permutation group and subgroup.

    properties (SetAccess = protected)
        group % (`.PermutationGroup`): Group
        subgroup % (`.PermutationGroup`): Subgroup of `.group`
    end

    properties %(Access = protected)
        groupChain % (`+replab.+bsgs.Chain`): Stabilizer chain of `.group` with monotonically increasing base
        subgroupChain % (`+replab.+bsgs.Chain`): Stabilzier chain of `.subgroup` with monotonically increasing base
    end

    methods

        function self = PermutationGroupRightCosets(group, subgroup)
            self.group = group;
            self.subgroup = subgroup;
            groupChain = group.niceGroup.chain;
            subgroupChain = subgroup.niceGroup.chain;
            if isempty(groupChain.B)
                groupChain = groupChain.mutableCopy;
                groupChain.baseChange(1);
                groupChain.makeImmutable;
            end
            subgroupChain = subgroupChain.mutableCopy;
            subgroupChain.baseChange(groupChain.B);
            subgroupChain.makeImmutable;
            self.groupChain = groupChain;
            self.subgroupChain = subgroupChain;
        end

        function t = findTransversalElement(self, g)
            group = self.groupChain;
            sub = self.subgroupChain.mutableCopy;
            n = group.n;
            L = group.length;
            t = g;
            for l = 1:L
                beta = group.B(l);
                tbeta = t(beta);
                b = replab.bsgs.minimalInOrbit(n, sub.strongGeneratorsForLevel(l), tbeta);
                sub.baseChange([sub.B(1:l-1) b]);
                uinv = sub.uinv(l, tbeta);
                t = uinv(t);
            end
        end

        function T = transversalMatrix(self)
            group = self.groupChain;
            subgroup = self.subgroupChain;
            n = group.n;
            L = group.length;
            stackG = zeros(n, L+1); % current group element
            stackG(:,1) = 1:n;
            stackS = cell(1, L+1);
            stackS{1} = subgroup; % subgroup chain
            stackB = cell(1, L+1);
            stackB{1} = sort(group.Delta{1}); % remaining points in orbit to check
            stackM = cell(1, L+1);
            stackM{1} = replab.bsgs.minimalMaskInOrbit(n, stackS{1}.strongGeneratorsForLevel(1));
            l = 1;
            T = [];
            while l >= 1
                g = stackG(:,l)';
                while 1
                    B = stackB{l};
                    if isempty(B)
                        l = l - 1;
                        break
                    end
                    b = B(1);
                    stackB{l} = B(2:end);
                    bg = g(b);
                    if stackM{l}(bg)
                        if l == L
                            T(:,end+1) = g(group.u(l, b));
                        else
                            stackS{l+1} = stackS{l}.stabilizer(bg);
                            gnext = g(group.u(l, b));
                            stackG(:,l+1) = gnext; % compose transversal element
                            orbit = sort(group.Delta{l+1});
                            [~, I] = sort(gnext(orbit));
                            stackB{l+1} = orbit(I);
                            stackM{l+1} = replab.bsgs.minimalMaskInOrbit(n, stackS{l+1}.strongGeneratorsForLevel(1));
                            l = l + 1;
                            break
                        end
                    end
                end
            end
            T = T';
        end

    end

end
