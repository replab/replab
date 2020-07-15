classdef LeftCosets < replab.CosetBase
% Describes the set of left cosets of a finite group
%
% Let $H$ be a subgroup of a group $G$. Then the left cosets are the sets $g H = \{ g h : h \in H \}$.
% The set of such left cosets is often written $G / H = \{ g H : g \in G \}$.

    methods

        function self = LeftCosets(group, subgroup)
            self@replab.CosetBase(group, subgroup);
        end

        function s = cardinality(self)
        % Returns the number of left cosets
        %
        % Returns:
        %   integer: Number of left cosets
            s = self.group.order / self.subgroup.order;
            assert(s <= 2^53 - 1);
            s = double(s);
        end

        function C = coset(self, g)
        % Returns the left coset containing the given element
        %
        % Args:
        %   g (element of `.group`): Group element
        %
        % Returns:
        %   `+replab.LeftCoset`: Left coset
            C = replab.LeftCoset(self.group, self.subgroup, self.cosetRepresentative(g));
        end

        function t = cosetRepresentative(self, g)
        % Returns the canonical coset representative corresponding to the given element
        %
        % If ``t = cosetRepresentative(g)``, then $t^{-1} g = h \in H$, with decomposition $g = t h$.
        %
        % Moreover ``L.cosetRepresentative(g) == L.cosetRepresentative(compose(g, h))`` for any $h \in H$.
        %
        % Finally, ``L.cosetRepresentative(g) == L.coset(g).representative``.
        %
        % Args:
        %   g (element of `.group`): Group element
        %
        % Returns:
        %   t (element of `.group`): Coset canonical representative
            g = self.isomorphism.imageElement(g);
            t = replab.bsgs.Cosets.leftRepresentative(self.subgroupChain, g);
            t = self.isomorphism.preimageElement(t);
        end

        function T = transversal(self)
        % Returns all the canonical representatives of cosets
        %
        % Returns:
        %   cell(1, \*) of `.group` elements: Transversal
            M = replab.bsgs.Cosets.leftTransversalAsMatrix(self.groupChain, self.subgroupChain);
            T = arrayfun(@(i) self.isomorphism.preimageElement(M(:,i)'), 1:self.cardinality, 'uniform', 0);
        end

        function C = elements(self)
        % Returns the set of left cosets as a cell array
        %
        % Returns:
        %   cell(1,\*) of `+replab.LeftCoset`: Set of left cosets
            T = self.transversal;
            C = cellfun(@(t) replab.LeftCoset(self.group, self.subgroup, t), T, 'uniform', 0);
            assert(isscalar(C{1}));
        end

% $$$         function mu = leftAction(group, subgroup)
% $$$         % Returns, as a morphism, the action of the given group of its left cosets
% $$$             nG = group.nGenerators;
% $$$             T = self.transversalAsMatrix;
% $$$             n = size(T, 1);
% $$$             images = cell(1, nG);
% $$$             for i = 1:nG
% $$$                 g = self.group.generator(i);
% $$$                 img = zeros(1, n);
% $$$                 for j = 1:n
% $$$                     gt = self.canonicalRepresentative(g(T(j,:)));
% $$$                     loc = replab.util.findRowInMatrix(gt, T);
% $$$                     % [ok, loc] = ismember(gt, T, 'rows');
% $$$                     assert(length(loc) == 1);
% $$$                     img(j) = loc;
% $$$                 end
% $$$                 images{i} = img;
% $$$             end
% $$$             Sn = replab.S(n);
% $$$             mu = self.group.morphismByImages(Sn, images);
% $$$         end

    end

end
