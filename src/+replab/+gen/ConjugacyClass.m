classdef ConjugacyClass < replab.ConjugacyClass & replab.gen.FiniteSet

    methods

        function self = ConjugacyClass(type, nice, niceIsomorphism)
            self@replab.gen.FiniteSet(type, nice, niceIsomorphism);
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

        function c1 = imap(self, f)
        % Maps this conjugacy class under an isomorphism
        %
        % Args:
        %   f (`.FiniteIsomorphism`): Isomorphism with ``self.group.isSubgroupOf(f.source)``
        %
        % Returns:
        %   `.ConjugacyClass`: The conjugacy class mapped under ``f``, expressed as a subset of ``f.image``
            if self.group.order < f.source.order
                f = f.restrictedSource(self.group);
            end
            c1 = replab.ConjugacyClass.make(f.target, f.imageElement(self.representative), f.imageGroup(self.representativeCentralizer));
        end

    end

    methods (Static)

        function c = make(group, element, elementCentralizer)
            if nargin < 3
                elementCentralizer = group.centralizer(element);
            end
            prmGroup = group.permutationGroup;
            prmElement = group.permutationIsomorphism.imageElement(element);
            [h1 g] = replab.bsgs.ConjugacyClasses.representative(prmGroup, prmElement);
            representative = group.permutationIsomorphism.preimageElement(h1);
            representativeCentralizer = elementCentralizer.leftConjugateGroup(group.permutationIsomorphism.preimageElement(g));
            c = replab.ConjugacyClass(group, representative, representativeCentralizer);
        end

    end

    methods (Access = protected)

        function E = computeElements(self)
            T = self.group.leftCosets(self.representativeCentralizer).transversal;
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
            B = self.group.findLeftConjugations(s, t, 'sCentralizer', sCentralizer);
            b = ~isempty(B);
        end

    end

end
