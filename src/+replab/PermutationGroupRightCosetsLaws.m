classdef PermutationGroupRightCosetsLaws < replab.Laws

    properties (SetAccess = protected)
        L % (`.PermutationGroupLeftCosets`): Left cosets structure
        G % (`.PermutationGroup`): Group
        H % (`.PermutationGroup`): Subgroup
    end
    methods
        function self = PermutationGroupRightCosetsLaws(L)
            self.L = L;
            self.G = L.group;
            self.H = L.subgroup;
        end
    end
    methods
        function law_transversal_is_ordered_(self)
            T = self.L.transversalMatrix;
            assertEqual(T, sortrows(T));
        end
        function law_representative_is_minimal_(self)
            T = self.L.transversalMatrix;
            U = self.L.subgroup.elements.toCell;
            for t = 1:size(T, 2)
                C = [];
                tel = T(t,:);
                for u = 1:length(U)
                    uel = U{u};
                    C(:,end+1) = self.G.compose(uel, tel);
                end
                C = sortrows(C');
                assertEqual(tel, C(1,:));
            end
        end
        function law_transversal_is_complete_(self)
            n = self.L.group.domainSize;
            T = self.L.transversalMatrix';
            U = self.L.subgroup.elements.toCell;
            els = zeros(n, 0);
            for t = 1:size(T, 2)
                tel = T(:,t)';
                for u = 1:length(U)
                    uel = U{u};
                    els(:,end+1) = self.G.compose(uel, tel);
                end
            end
            assert(size(els', 1) == self.G.order);
            assert(size(unique(els', 'rows'), 1) == self.G.order);
        end
    end
end
