classdef Backtrack < replab.Obj
% Performs backtracking search among the group elements
%
% It is a base class used for different but related tasks:
%
% - the computation of a subgroup of a group, whose elements satisfy a property,
% - the search of an element of a group satisfying a property.
%
% It attempts to reduce the search by considering only minimal elements in the double coset ``H g K`` of a group `.G`.
%
% We also use a object-oriented interface, instead of passing function handles for property and base image tests.
%
% In the subgroup search case, note that we assume that the elements that satisfy ``prop`` form a group. The current
% subgroup found is stored in a stabilizer chain `.subChain0`, and copies of it are made in `.HchainInBase0` and
% `.KchainInBase0`. Those chains are updated whenever a new subgroup generator is found.
% It is assumed that the user-given groups `.H` and `.K` are subgroups of the group we search for.
% Thus at the start of the search, `.subChain0` is initialized to the stabilizer chain of the closure of `.H` and `.K`.
%
% In the (coset) backtrack search case, `.subChain0` is not used, and `.HchainInBase0`, `.KchainInBase0` are initialized
% to the stabilizer chains of given groups `.H` and `.K`; and it is assumed that for any element ``g`` satisfying the property,
% then ``h g k`` also satisfies the property for any ``h`` in `.H` and ``k`` in `.K`.
%
% The property being tested should be implemented in `.prop`.
%
% The user can specify partial image tests. Those image tests are conducted relative to the user given partial base
% `.partialBase`. Note that the actual base starting with `.partialBase` may be longer; and that the code actually
% operates on a stabilizer chain whose redundant base points have been removed.
%
% However, the partial image tests are called as if the stabilizer chain started with `.partialBase`. Notably,
% `.partialBase`, `.base` and `.orbitSizes` correspond to the base with possible redundancy.
%
% Three references will be helpful in making sense of the code below:
%
% D. Holt, B. Eick, E. O’Brien, Handbook of Computational Group Theory. (CRC Press, 2005).
%
% which presents the most complete pruning criteria we know. An implementation following closely the pseudo-code
% in the book is available in Sympy 1.6, ``PermutationGroup.subgroup_search``; however, Sympy does not implement
% backtracking coset search.
%
% Our implementation follows more closely the overall structure given in
%
% G. Butler, Fundamental Algorithms for Permutation Groups. vol. 559 (Springer Berlin Heidelberg, 1991).
%
% which is more versatile, as their ``generate`` procedure directly provides search in double cosets.
%
% A `.debug` parameter can be set, which computes by brute enumeration the set of elements satisfying the property,
% which are then used to validate the partial base image tests.

    properties

        degree % (integer): Permutation group degree
        G % (`.PermutationGroup`): Group to search in
        H % (`.PermutationGroup`): Right cosets group
        K % (`.PermutationGroup`): Left cosets group

        partialBase % (integer(1,\*)): Prescribed base
        base % (integer(1,\*)): Actual base, starts with `.partialBase`
        orbitSizes % (integer(1,\*)): Orbit sizes corresponding to `.base`

        Gchain0 % (`.Chain`): Stabilizer chain of `.G` following `.base` with redundant points removed
        base0 % (integer(1,\*)): Prescribed base with redundant points removed
        baseLen0 % (integer): Equal to ``length(base0)``
        baseOrdering0 % (integer(1,domainSize+2)): Base ordering including two guard elements
        numRed0 % (integer(1,\*)): Number of redundant points in the original base preceding each base point of base0; number at the end is the number of points *after* the last nonredundant point

        HchainInBase0 % (`.Chain`): Current mutable BSGS chain for `.H` in base `.base0`
        KchainInBase0 % (`.Chain`): Current mutable BSHS chain for `.K` in base `.base0`
        subChain0 % (`.Chain`): Mutable chain of the subgroup being constructed in the base `.base0`, only used when searching for subgroups

        debug % (logical): Whether to do extra checks

    end

    methods

        function self = Backtrack(G, partialBase, H, K, debug)

            if nargin < 3 || isempty(H)
                H = G.trivialSubgroup;
            end
            if nargin < 4 || isempty(K)
                K = G.trivialSubgroup;
            end
            if nargin < 5 || isempty(debug)
                debug = false;
            else
                debug = debug;
            end

            degree = G.domainSize;

            % First block of properties
            self.degree = degree;
            self.G = G;
            self.H = H;
            self.K = K;

            % compute reduced base and reconstruct properties of the prescribed base
            Gchain0 = replab.bsgs.Backtrack.groupChainInBase(G, partialBase, true, true);
            base0 = Gchain0.base;
            if isempty(partialBase)
                base = base0;
            else
                firstNot = find(~ismember([base0 0], partialBase), 1);
                base = [partialBase base0(firstNot:end)];
            end
            assert(length(unique(base)) == length(base));
            assert(all(ismember(base0, base)));
            [numRed0 orbitSizes] = replab.bsgs.Backtrack.manageRedundantPoints(base, base0, Gchain0.orbitSizes);

            % Second block of properties

            self.partialBase = partialBase;
            self.base = base;
            self.orbitSizes = orbitSizes;
            self.Gchain0 = Gchain0;
            self.base0 = base0;
            self.baseLen0 = length(base0);
            self.baseOrdering0 = [replab.bsgs.Backtrack.computeBaseOrdering(degree, base0) degree+1 0];
            self.numRed0 = numRed0;

            self.HchainInBase0 = [];

            self.debug = debug;
        end

        function s = resultSet(self)
        % Computes explicitly the result set
        %
        % Returns:
        %   `+replab.+perm.Set`: Set of searched elements, used for debugging
            s = self.cached('resultSet', @() self.computeResultSet);
        end

        function s = computeResultSet(self)
            s = replab.perm.Set(self.degree);
            E = self.Gchain0.allElements;
            for i = 1:size(E, 2)
                if self.prop(E(:,i)')
                    s.insert(E(:,i));
                end
            end
        end

        function c = comparePermutations(self, lhs, rhs)
            diff = self.baseOrdering0(lhs(self.base0)) - self.baseOrdering0(rhs(self.base0))
            ind = find(diff == 0, 1);
            if isempty(ind)
                c = 0;
            else
                c = sign(diff(ind));
            end
        end

        function matrix = sortPermutations(self, matrix)
            baseImages = matrix(self.base0, :);
            sz = size(baseImages);
            baseImages = reshape(self.baseOrdering0(baseImages), sz); % in case vectors are funny
            [~, ind] = sortrows(baseImages');
            matrix = matrix(:, ind);
        end

        function res = find(self)
            self.subChain0 = [];
            self.HchainInBase0 = replab.bsgs.Backtrack.groupChainInBase(self.H, self.base0, false, true);
            self.KchainInBase0 = replab.bsgs.Backtrack.groupChainInBase(self.K, self.base0, false, true);
            res = self.generate(0, self.G.identity, self.HchainInBase0.mutableCopy, 1);
        end

        function res = subgroup(self)
            self.search(0);
            self.subChain0.makeImmutable;
            res = replab.PermutationGroup.fromChain(self.subChain0);
            self.subChain0 = [];
            self.HchainInBase0 = [];
            self.KchainInBase0 = [];
        end

        function search(self, s)
        % Search ``G^s`` for the subgroup of elements satisfying the property
        %
        % The call should start at level ``s = 0`` so the tests for the first redundant prescribed base points
        % can run.
        %
        % Args:
        %   s (integer): Level to search
            identity = 1:self.degree;
            if s == 0
                self.test0(0, identity, identity);
                self.search(s + 1);
            elseif s == self.baseLen0 + 1
                knownSubgroup = self.H.closure(self.K);
                self.subChain0 = replab.bsgs.Backtrack.groupChainInBase(knownSubgroup, self.base0, false, false);
                assert(isequal(self.subChain0.base, self.base0), 'Known subgroup is not a subgroup');
                self.HchainInBase0 = self.subChain0;
                self.KchainInBase0 = self.subChain0;
            else
                ok = self.test0(s, identity, identity);
                assert(ok); % for now, when looking for a subgroup
                self.search(s + 1);
                % note: compared to Butler, page 101, Algorithm 1, we need to remove the base point itself
                % there is no use in having the transversal u_s fixing the current base point u_s(base0(s)) = base0(s)
                orbit = self.Gchain0.Delta{s};
                orbit = orbit(2:end);
                [~, ind] = sort(self.baseOrdering0(orbit));
                orbit = orbit(ind);
                mask = replab.bsgs.minimalMaskInOrbit(self.degree, self.HchainInBase0.S, self.baseOrdering0);
                mu = self.computeMu(s, identity);
                ind = min(self.computeNu(s)+1-1, length(orbit)); % +1 because we removed the first point
                for gamma_s = orbit(1:ind)
                    if mask(gamma_s) && self.greaterThan(gamma_s, mu)
                        % this performs u = self.Gchain0.u(s, gamma_s)
                        ind = self.Gchain0.iDelta(gamma_s, s);
                        u = self.Gchain0.U{s}(:, ind)';

                        ok = self.test0(s, identity, u);
                        if ok
                            % Creates the stabilizer chain; cut at the current level, force
                            % stabilization of the current base image, and then pass the result
                            % to `.generate`, without cutting the first base point if it has
                            % changed (this is why `.generate` takes a ``HstabLevel`` argument).
                            Hstab = self.HchainInBase0.chainFromLevel(s, false);
                            if Hstab.stabilizes(1, gamma_s) % if the chain stabilizes already the point
                                HstabLevel = 1; % we don't need to do a base change
                            else
                                Hstab.changeBasePointAt(1, gamma_s); % otherwise, stabilize that point
                                HstabLevel = 2;
                            end
                            os = Hstab.orbitSizes; % test if the chain is trivial
                            if all(os(HstabLevel:end) == 1)
                                Hstab = []; % then pass the dummy value
                            end
                            found = self.generate(s+1, self.Gchain0.u(s, gamma_s), Hstab, HstabLevel);
                            if ~isempty(found)
                                if self.debug
                                    assert(self.resultSet.find(found(:)) > 0);
                                end
                                self.subChain0.stripAndAddStrongGenerator(found);
                                self.subChain0.randomizedSchreierSims([]);
                                self.HchainInBase0 = self.subChain0;
                                self.KchainInBase0 = self.subChain0;
                                mask = replab.bsgs.minimalMaskInOrbit(self.degree, self.HchainInBase0.S, self.baseOrdering0);
                                assert(~mask(gamma_s)); % verifies that the point is no longer minimal in its H-orbit
                            end
                        end
                    end
                end
            end

        end

        function found = generate(self, i, prevG, Hstab, HstabLevel)
        % Generate the elements of ``G^s``
        %
        % Those elements have base image ``[gamma(1) ... gamma(i-1)] = [prevG(base0(1)) ... prevG(base0(i-1))]``
        % and they may have the required property.
        % If one is found then we extend the subgroup of ``G^(s)`` of elements with property ``P`` that have
        % already been found, and return to search.
        %
        % Args:
        %   i (integer): Level to search
        %   prevG (permutation): Product ``u_1 ... u_{i-1}``
        %   Hstab (`.Chain` or ``[]``): Chain for `.H` with ``gamma(1) ... gamma(i-1)`` stabilized
        %                               However, one has to consider the subchain starting at the level
        %                               given below. Also, if the chain describes the trivial group,
        %                               we pass the value ``[]`` instead.
        %   HstabLevel (integer): Level at which to start consider ``Hstab``
        %
        % Returns:
        %   permutation or ``[]``: Element satisfying `.prop` if found, otherwise ``[]``
            if i == 0
                identity = 1:self.degree;
                self.test0(0, identity, identity);
                found = self.generate(i + 1, identity, Hstab, 1);
            elseif i == self.baseLen0 + 1
                if self.prop(prevG);
                    found = prevG;
                else
                    found = [];
                end
            else
                orbit_g = prevG(self.Gchain0.Delta{i});
                [~, ind] = sort(self.baseOrdering0(orbit_g));
                orbit_g = orbit_g(ind);
                if isempty(Hstab)
                    mask = true(1, self.degree);
                else
                    mask = replab.bsgs.minimalMaskInOrbit(self.degree, Hstab.S(:,Hstab.Sind(HstabLevel):end), self.baseOrdering0);
                end
                mu = self.computeMu(i, prevG);
                ind = min(self.computeNu(i)-1, length(orbit_g));
                for gamma_i = orbit_g(1:ind)
                    if mask(gamma_i) && self.baseOrdering0(gamma_i) > self.baseOrdering0(mu)
                        b = find(prevG == gamma_i);
                        ind = self.Gchain0.iDelta(b, i);
                        u = self.Gchain0.U{i}(:, ind)';
                        if self.test0(i, prevG, u)
                            if isempty(Hstab)
                                Hstab1 = Hstab;
                                % if trivial do nothing
                                HstabLevel1 = 0;
                            elseif Hstab.stabilizes(HstabLevel, gamma_i)
                                % if it already stabilizes, pass it on unchanged
                                Hstab1 = Hstab;
                                HstabLevel1 = HstabLevel;
                            else
                                % otherwise stabilize the new point
                                Hstab.changeBasePointAt(HstabLevel, gamma_i);
                                HstabLevel1 = HstabLevel + 1;
                                % if trivial, pass on the empty value
                                os = Hstab.orbitSizes;
                                if all(os(HstabLevel1:end) == 1)
                                    Hstab1 = [];
                                else
                                    Hstab1 = Hstab;
                                end
                            end
                            found = self.generate(i + 1, prevG(u), Hstab1, HstabLevel1);
                            if ~isempty(found)
                                return
                            end
                        end
                    end
                end
                found = [];
            end
        end

        function mu = computeMu(self, l, g)
            mu = self.degree + 2;
            beta = self.base0(l);
            for j = 1:l-1
                if any(self.KchainInBase0.Delta{j} == beta)
                    candidate = g(self.base0(j));
                    if self.greaterThan(candidate, mu)
                        mu = candidate;
                    end
                end
            end
        end

        function nu = computeNu(self, l)
            idx = length(self.Gchain0.Delta{l}) + 2 - length(self.KchainInBase0.Delta{l});
            if idx > length(self.Gchain0.Delta{l})
                nu = self.degree + 1;
            else
                nu = idx;
            end
        end

        function verifyPartialBaseUnsatisfied(self, l, gPrev, ul)
            pb = self.base(1:l);
            gamma = gPrev(ul(pb));
            matrix = self.resultSet.matrix;
            for i = 1:l
                mask = matrix(self.base(i), :) == gamma(i);
                matrix = matrix(:, mask);
            end
            assert(isempty(matrix));
        end

        function ok = test0(self, l0, gPrev, ul)
        % Tests if base images are possible at a given level
        %
        % Same as `.test` except it performs the test on the reduced base `.base0`
        %
        % ... should be called for ``l0 == 0`` as well
            l = sum(self.numRed0(1:l0)) + l0; % index in the original base
            if l > length(self.partialBase)
                ok = true;
                return
            end
            ok = true;
            if l0 >= 1
                ok = self.test(l, gPrev, ul);
            else
                ok = true;
            end
            if ~ok && self.debug
                self.verifyPartialBaseUnsatisfied(l, gPrev, ul);
            end
            nR = self.numRed0(l0+1);
            if nR > 0 && ok
                g = gPrev(ul);
                u = 1:length(ul); % identity
                for i = 1:nR
                    ok = self.test(l + i, g, u);
                    if ~ok
                        if self.debug
                            self.verifyPartialBaseUnsatisfied(l + i, g, u);
                        end
                        break
                    end
                end
            end
        end

        function ok = test(self, l, gPrev, ul)
        % Tests if base images are possible at a given level
        %
        % Note that the element ``g`` in Holt is given by ``g = compose(gPrev, ul)``.
        % Hence, the permutation being tested is ``gPrev(ul)``.
        %
        % This tests if there is any element in the searched set that has the partial base image
        % ``[g(base(1)) ... g(base(l))]`` of the preimage ``[base(1) ... base(l)]``.
        % False positives are possible, but not false negatives.
        %
        % Note that those tests are performed in the context of the original base `.base`, and thus
        % ``l`` has to be understood in that context. The rest of the code transparently handles the removed base
        % points by calling this method according to the original convention.
        %
        % Args:
        %   l (integer): Level of the test
        %   gPrev (permutation): Partial product ``u1 ... u{l-1}``
        %   ul (permutation): Current transversal element
        %
        % Returns:
        %   logical: False if there no element in the searched set that has the current partial base image
            ok = true; % default implementation
        end

        function l = greaterThan(self, x, y)
        % Point comparison using the base ordering
            l = self.baseOrdering0(x) > self.baseOrdering0(y);
        end

        function l = lessThan(self, x, y)
        % Point comparison using the base ordering
            l = self.baseOrdering0(x) < self.baseOrdering0(y);
        end

        function ok = prop(self, g)
        % Tests the property that elements of the searched subgroup satisfy
            error('Abstract');
        end

    end

    methods (Static)

        function [numRed0 orbitSizes] = manageRedundantPoints(base, base0, orbitSizes0)
        % Extract properties about a full base, from the data about a base with redundant points removed
        %
        % Args:
        %   base (integer(1,k)): Full base
        %   base0 (integer(1,k0)): Subset of ``base`` with redundant points removed
        %   orbitSizes0 (integer(1,k0)): Orbit sizes corresponding to base0
        %
        % Returns
        % -------
        %   numRed0:
        %     integer(1,k0+1): Number of redundant points preceding each base point, position ``k0+1`` counts the final redundant points
        %   orbitSizes:
        %     integer(1,k): Orbit sizes corresponding to ``base``
            k0 = length(base0);
            k = length(base);
            numRed0 = zeros(1, k0 + 1);
            i = 1;
            orbitSizes = zeros(1, k);
            for i0 = 1:k0
                if i > k
                    break
                end
                while base(i) ~= base0(i0)
                    orbitSizes(i) = 1;
                    numRed0(i0) = numRed0(i0) + 1;
                    i = i + 1;
                    if i > k
                        break
                    end
                end
                orbitSizes(i) = orbitSizes(i0);
                i = i + 1;
            end
            while i <= k
                numRed0(end) = numRed0(end) + 1;
                orbitSizes(i) = 1;
                i = i + 1;
            end
        end

        function bo = computeBaseOrdering(degree, base)
        % Returns a base ordering from a degree and given base
            rest = sort(setdiff(1:degree, base));
            bo = zeros(1, degree);
            bo(base) = 1:length(base);
            bo(rest) = length(base) + (1:length(rest));
        end

        function b = baseHasLexOrder(base)
        % Returns true if the base points are increasing
            b = all(base(2:end) > base(1:end-1));
        end

        function c = groupChainInBase(group, base, removeRedundant, immutable)
        % Get a stabilizer chain for the given group in the given base
        %
        % Args:
        %   group (`+replab.PermutationGroup`): Permutation group to get the stabilizer chain from
        %   base (integer(1,\*)): Prescribed base
        %   removeRedundant (logical): Whether to remove all redundant base points
        %   immutable (logical): Whether to return an immutable chain or a fresh mutable copy
            if replab.bsgs.Backtrack.baseHasLexOrder(base)
                c = group.lexChain.mutableCopy;
            else
                c = group.chain.mutableCopy;
            end
            c.baseChange(base, removeRedundant);
            if immutable
                c.makeImmutable;
            end
        end

    end

end
