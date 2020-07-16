classdef RightCoset < replab.Coset
% Describes a right coset of a finite group
%
% Let $g \in G$ be a coset representative and $H \le G$ a group. The right coset is $H g = \{ h g : h \in H\}$.
%
% Note:
%   This structure is implemented by referencing the permutation realizations of the group through a nice isomorphism.

    methods

        function self = RightCoset(group, canonicalRepresentative)
            self.group = group;
            self.representative = canonicalRepresentative;
            self.isomorphism = group.niceMorphism;
            self.groupChain = self.isomorphism.image.lexChain;
            self.type = group.type;
        end

    end

    methods (Static)


        function l = make(group, element)
            iso = group.niceMorphism;
            chain = iso.image.lexChain;
            permRep = replab.bsgs.Cosets.rightRepresentative(iso.image.lexChain, element);
            l = replab.RightCoset(group, iso.preimageElement(permRep));
        end

    end

    methods % Implementations

        % FiniteSet

        function s = cardinality(self)
        % Returns the size of this coset
        %
        % Returns:
        %   vpi: Coset size
            s = self.group.order;
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
            el = replab.bsgs.Cosets.rightRepresentative(self.groupChain, el);
            b = isequal(self.representative, el);
        end

        function E = computeElements(self)
        % Returns an indexed family of the elements of this coset
        %
        % Returns:
        %   `+replab.IndexedFamily`: Elements
            H = self.groupChain.allElements;
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

end
