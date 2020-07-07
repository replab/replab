classdef PermutationGroupRightCosetsLaws < replab.Laws

    properties (SetAccess = protected)
        R % (`.PermutationGroupRightCosets`): Right cosets structure
        G % (`.PermutationGroup`): Group
        H % (`.PermutationGroup`): Subgroup
    end
    methods
        function self = PermutationGroupRightCosetsLaws(R)
            self.R = R;
            self.G = R.group;
            self.H = R.subgroup;
        end
    end
    methods
        function law_transversal_is_ordered_(self)
            T = self.R.transversalAsMatrix;
            assertEqual(T, sortrows(T));
        end
        function law_representative_is_minimal_(self)
            T = self.R.transversalAsMatrix;
            Hels = self.R.subgroup.elements.toCell;
            for t = 1:size(T, 1)
                C = [];
                tel = T(t,:);
                for h = 1:length(Hels)
                    hel = Hels{h};
                    C(:,end+1) = self.G.compose(hel, tel);
                end
                C = sortrows(C');
                assertEqual(tel, C(1,:));
            end
        end
        function law_transversal_is_complete_(self)
            n = self.R.group.domainSize;
            T = self.R.transversalAsMatrix;
            Hels = self.R.subgroup.elements.toCell;
            els = zeros(0, n);
            for t = 1:size(T, 1)
                tel = T(t,:);
                for h = 1:length(Hels)
                    hel = Hels{h};
                    els(end+1,:) = self.G.compose(hel, tel);
                end
            end
            assert(size(els, 1) == self.G.order);
            els = unique(els, 'rows');
            assert(size(els, 1) == self.G.order);
        end
    end
end
