classdef ConjugacyClass < replab.FiniteSet
% Describes a conjugacy class of a finite group
%
% A conjugacy class containing the representative $r \in G$ is the set $\{g r g^{-1} : g \in G \}$.
%
% The centralizer of $r$ in $G$ is the subgroup $C_{G}(r) = \{ g r g^{-1} == r : g \in G \}$.
%
% Thus, the left cosets $G/C_{G}(r) = \{ g C_{G}(r) : g \in G \}$ are in one to one correspondence with
% the elements of the conjugacy class.

    properties (SetAccess = protected)
        group % (`+replab.FiniteGroup`): Group containing this conjugacy class
        representativeCentralizer % (`+replab.FiniteGroup`): Centralizer of `.representative` in `.group`
    end

    methods

        function self = ConjugacyClass(group, representative, representativeCentralizer)
            self.type = group.type;
            self.group = group;
            self.representative = representative;
            if nargin < 3 || isempty(representativeCentralizer)
                representativeCentralizer = self.group.centralizer(representative);
            end
            self.representativeCentralizer = representativeCentralizer;
        end

    end

    methods (Access = protected)

        function o = computeElementOrder(self)
            o = self.group.elementOrder(self.representative);
        end

    end

    methods

        function o = elementOrder(self)
        % Returns the order of the elements in this conjugacy class
        %
        % Returns:
        %   integer: Element order
            o = self.cached('elementOrder', @() self.computeElementOrder);
        end

        function c1 = imap(self, f, imageGroup, preserveLexOrder)
        % Maps this conjugacy class under an isomorphism
        %
        % Args:
        %   f (`.FiniteIsomorphism`): Isomorphism with ``self.group.isSubgroupOf(f.source)``
        %   imageGroup (`.FiniteGroup`, optional): Image of `.group` under ``f``, default ``[]`` (recompute)
        %   preserveLexOrder (logical, optional): Whether the isomorphism preserves the lexicographic order of group elements, default false
        %
        % Returns:
        %   `.ConjugacyClass`: The conjugacy class mapped under ``f``, expressed as a subset of ``f.image``
            if nargin < 3 || isempty(imageGroup)
                imageGroup = f.imageGroup(self.group);
            end
            if nargin < 4 || isempty(preserveLexOrder)
                preserveLexOrder = false;
            end
            if preserveLexOrder
                c1 = replab.ConjugacyClass(imageGroup, f.imageElement(self.representative), f.imageGroup(self.representativeCentralizer));
            else
                c1 = replab.ConjugacyClass.make(imageGroup, f.imageElement(self.representative), f.imageGroup(self.representativeCentralizer));
            end
        end

    end

    methods (Static)

        function c = make(group, element, elementCentralizer)
            if nargin < 3
                elementCentralizer = group.centralizer(element);
            end
            prmGroup = group.niceMorphism.image;
            prmElement = group.niceMorphism.imageElement(element);
            [h1 g] = replab.bsgs.ConjugacyClasses.representative(prmGroup, prmElement);
            representative = group.niceMorphism.preimageElement(h1);
            representativeCentralizer = elementCentralizer.leftConjugateGroup(group.niceMorphism.preimageElement(g));
            c = replab.ConjugacyClass(group, representative, representativeCentralizer);
        end

    end

    methods (Access = protected)

        function E = computeElements(self)
            T = self.group.leftCosetsOf(self.representativeCentralizer).transversal;
            E = cellfun(@(t) self.group.leftConjugate(t, self.representative), T, 'uniform', 0);
        end

        function n = computeNElements(self)
            n = self.group.order / self.representativeCentralizer.order;
        end

    end

    methods % Implementations

        % Str

        function s = shortStr(self, maxColumns)
            s = sprintf('ConjugacyClass of %s in %s', replab.shortStr(self.representative, maxColumns), replab.shortStr(self.group, maxColumns));
            if length(s) > maxColumns
                s = sprintf('ConjugacyClass of %s', replab.shortStr(self.representative, maxColumns));
            end
        end

        % Obj

        function l = laws(self)
            l = replab.laws.ConjugacyClassLaws(self);
        end

        % Domain

        function s = sample(self)
            t = self.group.sample;
            s = self.group.leftConjugate(t, self.representative);
        end

        % FiniteSet

        function s = nElements(self)
            s = self.cached('nElements', @() self.computeNElements);
        end

        function b = contains(self, t)
        % Returns whether the given element is part of this conjugacy class
        %
        % Args:
        %   t (element of `.group`): Group element
        %
        % Returns:
        %   logical: True if ``t`` is a member of this conjugacy classx
            s = self.representative;
            sCentralizer = self.representativeCentralizer;
            % We want to solve ``t == b s b^-1`` with ``s`` the representative
            B = self.group.findLeftConjugations(s, t, sCentralizer);
            b = ~isempty(B);
        end

    end

end
