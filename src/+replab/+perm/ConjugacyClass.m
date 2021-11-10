classdef ConjugacyClass < replab.ConjugacyClass & replab.PermutationFiniteSet

    methods

        function self = ConjugacyClass(group, representative)
            self.type = group.type;
            self.representative_ = representative;
            self.group = group;
        end

    end

    methods (Access = protected)

        function o = computeElementOrder(self)
            o = self.group.elementOrder(self.representative);
        end

        function n = computeNElements(self)
            n = self.group.order / self.representativeCentralizer.order;
        end

        function G = computeRepresentativeCentralizer(self)
            G = self.group.centralizer(self.representative);
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

        function s = nElements(self)
            s = self.cached('nElements', @() self.computeNElements);
        end

        % ConjugacyClass

        function o = elementOrder(self)
            o = self.cached('elementOrder', @() self.computeElementOrder);
        end

        function E = elementsSequence(self)
            mat = self.group.leftCosets(self.representativeCentralizer).transversalAsMatrix;
            y = self.representative;
            for i = 1:size(mat, 2)
                x = mat(:,i);
                mat(x,i) = x(y); % x y xInv
            end
            mat = sortrows(mat')';
            E = replab.perm.Sequence(mat);
        end

        function c = representativeCentralizer(self)
            c = self.cached('representativeCentralizer', @() self.computeRepresentativeCentralizer);
        end

    end

end
