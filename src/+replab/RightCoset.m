classdef RightCoset < replab.FiniteSet & replab.CosetBase
% Describes a right coset of a finite group
%
% Let $g \in G$ be the coset representative and $H \le G$ the subgroup. The right coset is $\{ h g : h \in H\}$.
%
% Note:
%   This structure is implemented by referencing the permutation realizations of the group and subgroup through
%   a nice isomorphism.

    methods (Access = protected)

        function self = RightCoset(group, subgroup, representative)
            self@replab.CosetBase(group, subgroup);
            self.type = group.type;
            self.representative = representative;
        end
    end

    methods

        function s = size(self)
        % Returns the size of this coset
        %
        % Returns:
        %   vpi: Coset size
            s = self.group.order/self.subgroup.order;
        end

        function b = contains(self, el)
        % Returns if this coset contains the given element
        %
        % Args:
        %   el (element of `.type`): Element to check
        %
        % Returns:
        %   logical: True if this coset contains the element
            if ~isempty(self.isomorphism)
                el = self.isomorphism.imageElement(el);
            end
            el = replab.bsgs.Cosets.rightRepresentative(self.subgroupChain, el);
            b = isequal(self.representative, el);
        end

        function E = computeElements(self)
        % Returns an indexed family of the elements of this coset
        %
        % Returns:
        %   `+replab.IndexedFamily`: Elements
            H = self.subgroupChain.allElements;
            g = self.representative;
            if ~isempty(self.isomorphism)
                g = self.isomorphism.imageElement(g);
            end
            matrix = zeros(size(H));
            for i = 1:size(H, 2)
                matrix(:,i) = H(g,i);
            end
            matrix = sortrows(matrix')';
            E = replab.indf.FiniteGroupIndexedFamily(matrix, self.isomorphism);
        end

    end

    methods (Static)


        function L = make(group, subgroup, element)
            assert(isa(group, 'replab.PermutationGroup')); % TODO
            L = replab.RightCoset(group, subgroup, replab.bsgs.Cosets.rightRepresentative(subgroup.lexChain, element));
        end

        function L = fromRepresentative(group, subgroup, representative)
            L = replab.RightCoset(group, subgroup, representative);
        end

    end

end
