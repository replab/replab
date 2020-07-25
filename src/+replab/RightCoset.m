classdef RightCoset < replab.Coset
% Describes a right coset of a finite group
%
% Let $g \in G$ be a coset representative and $H \le G$ a group. The right coset is $H g = \{ h g : h \in H\}$.
%
% Note:
%   This structure is implemented by referencing the permutation realizations of the group through a nice isomorphism.

    methods

        function self = RightCoset(group, canonicalRepresentative, parent)
            self.type = parent.type;
            self.parent = parent;
            self.isomorphism = parent.niceMorphism;
            self.group = group;
            self.representative = canonicalRepresentative;
            self.groupChain = parent.niceMorphism.imageGroup(group).lexChain;
        end

    end

    methods (Static)

        function r = make(group, element, parent)
            if nargin < 3 || isempty(parent)
                parent = group.closure(element);
            end
            if group.isNormalizedBy(element)
                r = replab.NormalCoset(group, element, parent);
                return
            end
            chain = parent.niceMorphism.imageGroup(group).lexChain;
            permRep = replab.bsgs.Cosets.rightRepresentative(chain, parent.niceMorphism.imageElement(element));
            r = replab.RightCoset(group, parent.niceMorphism.preimageElement(permRep), parent);
        end

    end

    methods (Access = protected)

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

    methods % Implementations

        % Domain

        function s = sample(self)
            s = self.type.compose(self.group.sample, self.representative);
        end

        % FiniteSet

        function s = size(self)
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

    end

end
