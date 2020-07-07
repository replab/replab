classdef PermutationGroupLeftCosets < replab.Str

    properties (SetAccess = protected)
        group % (`.PermutationGroup`): Group
        subgroup % (`.PermutationGroup`): Subgroup of `.group`
    end

    properties (Access = protected)
        inverse % (`.PermutationGroupRightCosets`): Right cosets
    end

    properties %(Access = protected)
        groupChain % (`+replab.+bsgs.Chain`): Stabilizer chain of `.group` with monotonically increasing base
        subgroupChain % (`+replab.+bsgs.Chain`): Stabilzier chain of `.subgroup` with monotonically increasing base
    end

    methods

        function self = PermutationGroupLeftCosets(group, subgroup, inverse)
            self.group = group;
            self.subgroup = subgroup;
            if nargin < 3 || isempty(inverse)
                self.inverse = replab.PermutationGroupRightCosets(group, subgroup);
            end
        end

        function t = findTransversalElement(self, g)
        % Finds the transversal element corresponding to a given group element
            t = self.group.inverse(self.inverse.findTransversalElement(self.group.inverse(g)));
            % we have the left coset decomposition g = h t, where h \in subgroup and t is a transversal
            % then g^-1 = t^-1 h^-1, and we ask for the right transversal element t^-1 corresponding to g^-1
        end

        function T = transversalMatrix(self)
            Tinv = self.inverse.transversalMatrix;
            T = zeros(size(Tinv));
            for i = 1:size(Tinv, 1)
                T(i,Tinv(i,:)) = 1:self.group.domainSize;
            end
        end

    end

end
