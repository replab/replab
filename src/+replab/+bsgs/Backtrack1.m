classdef Backtrack1 < replab.Str

    properties

        degree % (integer): Permutation group degree
        group % (`.Chain`): Group to search in
        knownSubgroup % (`.Chain`): Known subgroup of the searched for subgroup

        base % (integer(1,\*)): Prescribed base
        debug % (logical): Whether to do extra checks

        group0 % (`.Chain`): Stabilizer chain of `.group` following `.base` with redundant points removed
        base0 % (integer(1,\*)): Prescribed base with redundant points removed
        baseOrdering0 % (integer(1,domainSize+2)): Base ordering including two guard elements
        numRed0 % (integer(1,\*)): Number of redundant points in the original base preceding each base point of base0; number at the end is the number of points *after* the last nonredundant point

        KinBase0 % (`.Chain`): Current mutable BSGS chain for the subgroup in base `.base0`

        debugSet % (`+replab.+perm.Set`): Set of searched elements, used for debugging
    end

    methods

        function self = Backtrack1(group, base, knownSubgroup, debug)
            degree = group.n;
            self.degree = degree;
            self.group = group;
            if nargin < 3 || isempty(knownSubgroup)
                knownSubgroup = replab.bsgs.Chain(degree);
                knownSubgroup.makeImmutable;
            end
            if nargin < 4
                debug = false;
            end
            self.knownSubgroup = knownSubgroup;
            self.base = base;
            group0 = group.mutableCopy;
            group0.baseChange(base, true);
            group0.makeImmutable;
            self.group0 = group0;
            base0 = group0.base;
            self.base0 = base0;
            rest = sort(setdiff(1:degree, base0));
            baseOrdering0 = zeros(1, degree);
            baseOrdering0(base0) = 1:length(base0);
            baseOrdering0(rest) = length(base0) + (1:length(rest));
            self.baseOrdering0 = [baseOrdering0 degree+1 0];
            if group0.length == 0
                return
            end
            numRed0 = zeros(1, length(base0) + 1);
            i = 1;
            for i0 = 1:length(base0)
                while base(i) ~= base0(i0)
                    numRed0(i0) = numRed0(i0) + 1;
                    i = i + 1;
                end
                i = i + 1;
            end
            numRed0(end) = length(base) - i + 1;
            self.numRed0 = numRed0;
            self.KinBase0 = [];
            if debug
                self.debug = true;
                debugSet = replab.perm.Set(degree);
                E = group0.allElements;
                for i = 1:size(E, 2)
                    if self.prop(E(:,i)')
                        debugSet.insert(E(:,i)');
                    end
                end
                self.debugSet = debugSet;
            else
                self.debug = false;
                self.debugSet = [];
            end
            % verify ordering of groups DEBUG
            ch = group0;
            while ch.length > 1
                ch1 = ch.stabilizer(ch.B(1));
                el1 = self.sortPermutations(ch1.allElements);
                el = self.sortPermutations(ch.allElements);
                assert(isequal(el(:,1:size(el1,2)), el1));
                ch = ch1;
            end
            % END DEBUG
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

        function res = subgroup(self)
            disp('start');
            self.test0(0, 1:self.degree, 1:self.degree);
            self.search2(1);
            self.KinBase0.makeImmutable;
            res = replab.PermutationGroup.fromChain(self.KinBase0);
            self.KinBase0 = [];
        end

        function search1(self, s)
        % Search ``G^s`` for the subgroup K of elements satisfying the property
        %
        % Args:
        %   s (integer): Level to search
            if s == length(self.base0) + 1
                c = self.knownSubgroup.mutableCopy;
                c.baseChange(self.base0);
                assert(isequal(c.base, self.base0), 'Known subgroup is not a subgroup');
                self.KinBase0 = c;
            else
                self.search1(s + 1);
                identity = 1:self.degree;
                % note: compared to Butler, page 101, Algorithm 1, we need to remove the base point itself
                % there is no use in having the transversal u_s fixing the current base point u_s(base0(s)) = base0(s)
                orbit = self.sort(self.group0.Delta{s}(2:end));
                mask = replab.bsgs.minimalMaskInOrbit(self.degree, self.KinBase0.S, self.baseOrdering0);
                for gamma_s = orbit
                    cond1 = mask(gamma_s);
                    cond2 = replab.bsgs.minimalInOrbit(self.degree, self.KinBase0.S, gamma_s, self.baseOrdering0) == gamma_s;
                    assert(cond1 == cond2);
                    if mask(gamma_s)
                        u = self.group0.u(s, gamma_s);
                        found = self.generate1(s, s+1, self.group0.u(s, gamma_s));
                        if found
                            mask = replab.bsgs.minimalMaskInOrbit(self.degree, self.KinBase0.S, self.baseOrdering0);
                            assert(~mask(gamma_s)); % DEBUG asserts the point is no longer minimal in its K-orbit
                        end
                    end
                end
            end
        end

        function found = generate1(self, s, i, prevG)
        % Generate the elements of ``G^s``
        %
        % Those elements have base image ``[gamma(1) ... gamma(i-1)] = [prevG(base0(1)) ... prevG(base0(i-1))]``
        % and they may have the required property.
        % If one is found then we extend ``K``, the subgroup of ``G^(s)`` of elements with property ``P`` that have
        % already been found, and return to search.
            if i == length(self.base0) + 1
                found = self.prop(prevG);
                if found
                    self.KinBase0.stripAndAddStrongGenerator(prevG);
                    self.KinBase0.randomizedSchreierSims([]);
                end
            else
                orbit_g = self.sort(prevG(self.group0.Delta{i}));
                for gamma_i = orbit_g
                    b = find(prevG == gamma_i);
                    u = self.group0.u(i, b);
                    assert(prevG(u(self.base0(i))) == gamma_i);
                    found = self.generate1(s, i + 1, prevG(u));
                    if found
                        return
                    end
                end
            end
        end

        function search2(self, s)
        % Search ``G^s`` for the subgroup K of elements satisfying the property
        %
        % Args:
        %   s (integer): Level to search
            if s == length(self.base0) + 1
                c = self.knownSubgroup.mutableCopy;
                c.baseChange(self.base0);
                assert(isequal(c.base, self.base0), 'Known subgroup is not a subgroup');
                self.KinBase0 = c;
            else
                identity = 1:self.degree;
                ok = self.test0(s, identity, identity);
                assert(ok); % for now, when looking for a subgroup
                self.search2(s + 1);
                % note: compared to Butler, page 101, Algorithm 1, we need to remove the base point itself
                % there is no use in having the transversal u_s fixing the current base point u_s(base0(s)) = base0(s)
                orbit = self.sort(self.group0.Delta{s}(2:end));
                mask = replab.bsgs.minimalMaskInOrbit(self.degree, self.KinBase0.S, self.baseOrdering0);
                for gamma_s = orbit
                    cond1 = mask(gamma_s);
                    cond2 = replab.bsgs.minimalInOrbit(self.degree, self.KinBase0.S, gamma_s, self.baseOrdering0) == gamma_s;
                    assert(cond1 == cond2);
                    if mask(gamma_s)
                        u = self.group0.u(s, gamma_s);
                        ok = self.test0(s, identity, u);
                        if ok
                            found = self.generate2(s, s+1, self.group0.u(s, gamma_s));
                            if found
                                mask = replab.bsgs.minimalMaskInOrbit(self.degree, self.KinBase0.S, self.baseOrdering0);
                                assert(~mask(gamma_s)); % DEBUG asserts the point is no longer minimal in its K-orbit
                            end
                        end
                    end
                end
            end
        end

        function found = generate2(self, s, i, prevG)
        % Generate the elements of ``G^s``
        %
        % Those elements have base image ``[gamma(1) ... gamma(i-1)] = [prevG(base0(1)) ... prevG(base0(i-1))]``
        % and they may have the required property.
        % If one is found then we extend ``K``, the subgroup of ``G^(s)`` of elements with property ``P`` that have
        % already been found, and return to search.
            if i == length(self.base0) + 1
                found = self.prop(prevG);
                if found
                    self.KinBase0.stripAndAddStrongGenerator(prevG);
                    self.KinBase0.randomizedSchreierSims([]);
                end
            else
                orbit_g = self.sort(prevG(self.group0.Delta{i}));
                for gamma_i = orbit_g
                    b = find(prevG == gamma_i);
                    u = self.group0.u(i, b);
                    if self.test0(i, prevG, u)
                        assert(prevG(u(self.base0(i))) == gamma_i);
                        found = self.generate2(s, i + 1, prevG(u));
                        if found
                            return
                        end
                    end
                end
                found = false;
            end
        end

        function ok = test0(self, l0, gPrev, ul)
        % Tests if base images are possible at a given level
        %
        % Same as `.test` except it performs the test on the reduced base `.base0`
        %
        % ... should be called for ``l0 == 0`` as well

            l = sum(self.numRed0(1:l0)) + l0; % index in the original base
            assert(l0 == 0 || self.base(l) == self.base0(l0));
            ok = true;
            if l0 >= 1
                ok = self.test(l, gPrev, ul);
            else
                ok = true;
            end
            nR = self.numRed0(l0+1);
            if ~ok || nR == 0
                return
            end
            g = gPrev(ul);
            u = 1:length(ul); % identity
            for i = 1:nR
                ok = self.test(l + i, g, u);
                if ~ok
                    return
                end
            end
        end


        function ok = test(self, l, gPrev, ul)
        % Tests if base images are possible at a given level
        %
        % Note that the element ``g`` in Holt is given by ``g = compose(gPrev, ul)``.
        %
        % This tests if there is any element in the searched set that has the partial base image
        % ``[g(base(1)) ... g(base(l))]``. False negatives are possible, but not false positives.
        %
        % Note that those tests are performed in the context of the original base `.base`, and thus
        % ``l`` has to be understood in that context. The code transparently handles the removed base
        % points by calling the tests according to the original conventionx.
        %
        % Args:
        %   l (integer): Level of the test
        %   gPrev (permutation): Partial product ``u1 ... u{l-1}``
        %   ul (permutation): Current transversal element
        %
        % Returns:
        %   logical: False if there no element in the searched set that has the current partial base image
            error('Abstract');
        end

        function l = greaterThan(self, x, y)
        % Point comparison using the base ordering
            l = self.baseOrdering0(x) > self.baseOrdering0(y);
        end

        function l = lessThan(self, x, y)
        % Point comparison using the base ordering
            l = self.baseOrdering0(x) < self.baseOrdering0(y);
        end

        function t = sort(self, s)
        % Sorts a sequence of points under the ordering
            [~, ind] = sort(self.baseOrdering0(s));
            t = s(ind);
        end

        function ok = prop(self, g)
        % Tests the property that elements of the searched subgroup satisfy
            error('Abstract');
        end

    end

end
