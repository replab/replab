classdef ConjugacyClass < replab.ConjugacyClass & replab.PermutationFiniteSet

    methods

        function self = ConjugacyClass(group, representative, representativeCentralizer)
            self.type = group.type;
            self.representative = representative;
            self.group = group;
            if nargin >= 3 && ~isempty(representativeCentralizer)
                self.cache('representativeCentralizer', representativeCentralizer);
            end
        end

    end

    methods (Static)

        function c = make(group, element, elementCentralizer)
            [representative, g] = replab.bsgs.ConjugacyClasses.representative(group, element);
            if nargin >= 3
                representativeCentralizer = elementCentralizer.leftConjugateGroup(g);
                c = replab.perm.ConjugacyClass(group, representative, representativeCentralizer);
            else
                c = replab.perm.ConjugacyClass(group, representative);
            end
        end

    end

    methods (Access = protected)

        function o = computeElementOrder(self)
            o = self.group.elementOrder(self.representative);
        end

% $$$         function E = computeElementsSequence(self)
% $$$             T = self.group.leftCosets(self.representativeCentralizer).transversal;
% $$$             E = cellfun(@(t) self.group.leftConjugate(t, self.representative), T, 'uniform', 0);
% $$$         end

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

        function c = representativeCentralizer(self)
            c = self.cached('representativeCentralizer', @() self.computeRepresentativeCentralizer);
        end

    end

end
