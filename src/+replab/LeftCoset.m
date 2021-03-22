classdef LeftCoset < replab.Coset
% Describes a left coset of a finite group
%
% Let $g \in G$ be a coset representative and $H \le G$ a group. The left coset is $g H = \{ g h : h \in H\}$.
%
% Note:
%   This structure is implemented by referencing the permutation realizations of the group through a nice isomorphism.

    methods

        function self = LeftCoset(group, canonicalRepresentative, parent)
            self.type = parent.type;
            self.parent = parent;
            self.isomorphism = parent.niceMorphism;
            self.group = group;
            self.representative = canonicalRepresentative;
            self.groupChain = parent.niceMorphism.imageGroup(group).lexChain;
        end

    end

    methods (Static)

        function l = make(group, element, parent)
            if nargin < 3 || isempty(parent)
                parent = group.closure(element);
            end
            if group.isNormalizedBy(element)
                l = replab.NormalCoset(group, element, parent);
                return
            end
            chain = parent.niceMorphism.imageGroup(group).lexChain;
            permRep = replab.bsgs.Cosets.leftRepresentative(chain, parent.niceMorphism.imageElement(element));
            l = replab.LeftCoset(group, parent.niceMorphism.preimageElement(permRep), parent);
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
            matrix = sortrows(g(H'))';
            E = replab.indf.FiniteGroupIndexedFamily(matrix, self.isomorphism);
        end

    end

    methods % Implementations

        % Domain

        function s = sample(self)
            s = self.type.compose(self.representative, self.group.sample);
        end

        % FiniteSet

        function s = nElements(self)
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
            el = replab.bsgs.Cosets.leftRepresentative(self.groupChain, el);
            b = isequal(self.representative, el);
        end

        % Coset

        function [l, r] = factorizeShortRepresentativeLetters(self)
        % Returns a tentatively short word corresponding to an element of this coset
        %
        % An effort is made to identify a short word, but without optimality guarantees.
        %
        % Returns
        % -------
        %   l: integer(1,\*)
        %     Letters of the word representing an element of ``self``
        %   r: element of `.group`
        %     Represented coset element
            error('WIP');
        end
    end

end
