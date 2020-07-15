classdef LeftCoset < replab.CosetBase & replab.FiniteSet
% Describes a left coset of a finite group
%
% Let $g \in G$ be the coset representative and $H \le G$ the subgroup. The left coset is $\{ g h : h \in H\}$.
%
% Note:
%   This structure is implemented by referencing the permutation realizations of the group and subgroup through
%   a nice isomorphism.

    methods

        function self = LeftCoset(group, subgroup, canonicalRepresentative)
            self@replab.CosetBase(group, subgroup);
            self.type = group.type;
            self.representative = canonicalRepresentative;
        end

        function s = cardinality(self)
        % Returns the size of this coset
        %
        % Returns:
        %   vpi: Coset size
            s = self.subgroup.order;
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
            el = replab.bsgs.Cosets.leftRepresentative(self.subgroupChain, el);
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
            matrix = sortrows(g(H'))';
            E = replab.indf.FiniteGroupIndexedFamily(matrix, self.isomorphism);
        end

    end

end
