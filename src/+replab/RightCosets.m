classdef RightCosets < replab.CosetBase
% Describes the set of right cosets of a nice finite group
%
% Let $H$ be a subgroup of a group $G$. Then the right cosets are the sets $H g = \{ h g : h \in H \}$.
% The set of such right cosets is often written $H \\ G = \{ H g : g \in G \}$.

    methods

        function self = RightCosets(group, subgroup)
            self@replab.CosetBase(group, subgroup);
        end

        function s = cardinality(self)
        % Returns the number of right cosets
        %
        % Returns:
        %   integer: Number of right cosets
            s = self.group.order / self.subgroup.order;
            assert(s < 2^53 - 1);
            s = double(s);
        end

        function C = coset(self, g)
        % Returns the right coset containing the given element
        %
        % Args:
        %   g (element of `.group`): Group element
        %
        % Returns:
        %   `+replab.RightCoset`: Right coset
            C = replab.RightCoset(self.group, self.subgroup, self.cosetRepresentative(g));
        end

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
            g = self.isomorphism.imageElement(g);
            t = replab.bsgs.Cosets.rightRepresentative(self.subgroupChain, g);
            t = self.isomorphism.preimageElement(t);
        end

        function T = transversal(self)
        % Returns all the canonical representatives of cosets
        %
        % Returns:
        %   cell(1, \*) of `.group` elements: Transversal
            M = replab.bsgs.Cosets.rightTransversalMatrix(self.groupChain, self.subgroupChain);
            T = arrayfun(@(i) self.isomorphism.preimageElement(M(:,i)'), 1:self.cardinality, 'uniform', 0);
        end

        function C = elements(self)
        % Returns the set of right cosets as a cell array
        %
        % Returns:
        %   cell(1,\*) of `+replab.RightCoset`: Set of right cosets
            C = cellfun(@(t) replab.RightCoset(self.group, self.subgroup, t), self.transversal, 'uniform', 0);
        end

    end

end
