classdef RightCosets < replab.Cosets
% Describes the set of right cosets of a nice finite group
%
% Let $H$ be a subgroup of a group $G$. Then the right cosets are the sets $H g = \{ h g : h \in H \}$.
% The set of such right cosets is often written $H \\ G = \{ H g : g \in G \}$.

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.RightCosetsLaws(self);
        end

    end

    methods

        function t = cosetRepresentative(self, g)
        % Returns the coset representative corresponding to the given element
        %
        % If ``t = cosetRepresentative(g)``, then $g t^{-1} = h \in H$, with the decomposition $g = h t$.
        %
        % Moreover ``R.cosetRepresentative(g) == R.cosetRepresentative(compose(h, g))`` for any $h \in H$.
        %
        % Finally, ``R.cosetRepresentative(g) == R.coset(g).representative``.
        %
        % Args:
        %   g (element of `.group`): Group element
        %
        % Returns:
        %   t (element of `.group`): Coset canonical representative
            error('Abstract');
        end

        function C = elements(self)
        % Returns the set of right cosets as a cell array
        %
        % Returns:
        %   cell(1,\*) of `+replab.RightCoset`: Set of right cosets
            C = cellfun(@(t) replab.RightCoset(self.subgroup, t, self.group), self.transversal, 'uniform', 0);
        end

        function d = leftCosetsBy(self, subgroup1)
            d = self.group.doubleCosets(self.subgroup, subgroup1);
        end

        function d = mrdivide(self, subgroup1)
            d = self.leftCosetsBy(subgroup1);
        end

        function s = nElements(self)
        % Returns the number of right cosets
        %
        % Returns:
        %   integer: Number of right cosets
            s = self.group.order / self.subgroup.order;
            assert(s < 2^53 - 1);
            s = double(s);
        end

        function T = transversal(self)
        % Returns all the canonical representatives of cosets
        %
        % Returns:
        %   cell(1, \*) of `.group` elements: Transversal
            error('Abstract');
        end

    end

end
