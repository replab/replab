classdef PermutationGroupLeftCosetsLaws < replab.LeftCosetsLaws
% Laws for `.PermutationGroupLeftCosets`
%
% Property types
%
% L is now a `.PermutationGroupLeftCosets`, G and H are `.PermutationGroup`
    methods
        function self = PermutationGroupLeftCosetsLaws(L)
            assert(isa(L, 'replab.PermutationGroupLeftCosets'));
            self@replab.LeftCosetsLaws(L);
        end
    end
    methods
        function law_transversal_is_ordered_(self)
            T = self.L.transversalAsMatrix;
            assertEqual(T, sortrows(T));
        end
        function law_representative_is_minimal_(self)
            T = self.L.transversalAsMatrix;
            Hels = self.L.subgroup.elements.toCell;
            for t = 1:size(T, 1)
                C = [];
                tel = T(t,:);
                for h = 1:length(Hels)
                    hel = Hels{h};
                    C(end+1,:) = self.G.compose(tel, hel);
                end
                C = sortrows(C);
                assertEqual(tel, C(1,:));
            end
        end
        function law_transversal_is_complete_(self)
            n = self.L.group.domainSize;
            T = self.L.transversalAsMatrix;
            Hels = self.L.subgroup.elements.toCell;
            els = zeros(0, n);
            for t = 1:size(T, 1)
                tel = T(t,:);
                for h = 1:length(Hels)
                    hel = Hels{h};
                    els(end+1,:) = self.G.compose(tel, hel);
                end
            end
            assert(size(els, 1) == self.G.order);
            els = unique(els, 'rows');
            assert(size(els, 1) == self.G.order);
        end
    end
end
