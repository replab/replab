classdef PermutationGroupLeftCosets < replab.LeftCosets

    properties (Access = protected)
        inverse % (`.PermutationGroupRightCosets`): Right cosets
    end

    methods

        function self = PermutationGroupLeftCosets(group, subgroup, inverse)
            self.group = group;
            self.subgroup = subgroup;
            if nargin < 3 || isempty(inverse)
                self.inverse = replab.PermutationGroupRightCosets(group, subgroup);
            end
        end

        function t = canonicalRepresentative(self, g)
            t = self.group.inverse(self.inverse.findTransversalElement(self.group.inverse(g)));
            % we have the left coset decomposition g = h t, where h \in subgroup and t is a transversal
            % then g^-1 = t^-1 h^-1, and we ask for the right transversal element t^-1 corresponding to g^-1
        end

        function T = transversalAsMatrix(self)
            Tinv = self.inverse.transversalAsMatrix;
            T = zeros(size(Tinv));
            for i = 1:size(Tinv, 1)
                T(i,Tinv(i,:)) = 1:self.group.domainSize;
            end
        end

        function C = cosetAsMatrix(self, g)
        % Returns the given coset as a matrix, with elements in rows
        %
        % Args:
        %   g (permutation): Coset representative
        %
        % Returns:
        %   integer(\*,\*): Coset matrix
            H = self.subgroupChain.allElements;
            C = g(H');
            C = sortrows(C');
        end

        function T = transversal(self)
            M = self.transversalAsMatrix;
            T = arrayfun(@(i) M(i,:), 1:size(M,1), 'uniform', 0);
        end

        function C = coset(self, g)
            M = self.cosetAsMatrix(g);
            N = size(M, 1);
            C = arrayfun(@(i) M(i,:), 1:size(M,1), 'uniform', 0);
        end

    end

end
