classdef PermutationGroupLeftCosets < replab.LeftCosets
% Enables computations with the left cosets of a permutation group
%
% The transversal elements are given by the minimal elements of the coset under lexicographic ordering;
% the left cosets themselves are ordered by their transversal elements, in lexicographic ordering too.
%
% Thus, computations involving cosets are deterministic, and do not depend on the details of
% the construction of the permutation group and subgroup.
%
% See `.LeftCosets` for the basic information.

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
            sub = self.inverse.subgroupChain;
            L = sub.length;
            t = g;
            for l = 1:L
                orbit = sub.Delta{l};
                [~, i] = min(g(orbit));
                t = t(sub.u(l, orbit(i)));
            end
        end

        function T = transversalAsMatrix(self)
            Tinv = self.inverse.transversalAsMatrix;
            T = zeros(size(Tinv));
            n = self.group.domainSize;
            g = zeros(1, n);
            for i = 1:size(Tinv, 1)
                g(Tinv(i,:)) = 1:n; % invert
                T(i,:) = self.canonicalRepresentative(g);
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
            H = self.inverse.subgroupChain.allElements;
            C = g(H');
            C = sortrows(C);
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

        function mu = action(self)
            nG = self.group.nGenerators;
            T = self.transversalAsMatrix;
            n = size(T, 1);
            images = cell(1, nG);
            for i = 1:nG
                g = self.group.generator(i);
                img = zeros(1, n);
                for j = 1:n
                    gt = self.canonicalRepresentative(g(T(j,:)));
                    [ok, loc] = ismember(gt, T, 'rows');
                    assert(ok);
                    img(j) = loc;
                end
                images{i} = img;
            end
            Sn = replab.S(n);
            mu = self.group.morphismByImages(Sn, images);
        end

    end

end
