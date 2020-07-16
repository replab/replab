classdef Backtrack1 < replab.Str

    properties

        degree % (integer): Permutation group degree
        group % (`.Chain`): Group to search in
        knownSubgroup % (`.Chain`): Known subgroup of the searched for subgroup
        baseOrdering % (integer(1,domainSize+2)): Base ordering including two guard elements

        base % (integer(1,\*)): Prescribed base
        debug % (logical): Whether to do extra checks

        group0 % (`.Chain`): Stabilizer chain of `.group` following `.base` with redundant points removed
        base0 % (integer(1,\*)): Prescribed base with redundant points removed
        numRed0 % (integer(1,\*)): Number of redundant points in the original base following each base point of base0

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
            rest = sort(setdiff(1:degree, base));
            baseOrdering = zeros(1, degree);
            baseOrdering(base) = 1:length(base);
            baseOrdering(rest) = length(base) + (1:length(rest));
            self.baseOrdering = [baseOrdering degree+1 0];
            self.base = base;
            group0 = group.mutableCopy;
            group0.baseChange(base, true);
            group0.makeImmutable;
            self.group0 = group0;
            base0 = group0.base;
            self.base0 = base0;
            if group0.length == 0
                return
            end
            baseStart = 1;
            while base(baseStart) ~= base0(1)
                baseStart = baseStart + 1;
            end
            i = baseStart;
            numRed0 = zeros(1, length(base0));
            for i0 = 1:length(base0)-1
                i = i + 1;
                while base(i) ~= base0(i0+1)
                    numRed0(i0) = numRed0(i0) + 1;
                    i = i + 1;
                end
            end
            numRed0(end) = length(base) - i;
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
        end

        function res = subgroup(self)
            self.search(1);
            self.KinBase0.makeImmutable;
            res = replab.PermutationGroup.fromChain(self.KinBase0);
            self.KinBase0 = [];
        end

        function search(self, s)
        % Search ``G^s`` for the subgroup K of elements satisfying the property
        %
        % Args:
        %   s (integer): Level to search
            if s == length(self.base0 + 1)
                c = self.knownSubgroup.mutableCopy;
                c.baseChange(self.base0);
                assert(isequal(c.base, self.base0), 'Known subgroup is not a subgroup');
                self.KinBase0 = c;
            else
                self.search(s + 1);
                identity = 1:self.degree;
                orbit = self.sort(self.group0.Delta{s});
                %assert(isequal(self.KinBase0.base, self.base0));
                mask = replab.bsgs.minimalMaskInOrbit(self.degree, self.KinBase0.S, self.baseOrdering);
                for gamma_s = orbit
                    cond1 = mask(gamma_s);
                    cond2 = replab.bsgs.minimalInOrbit(self.degree, self.KinBase0.S, gamma_s, self.baseOrdering) == gamma_s;
                    assert(cond1 == cond2);
                    if mask(gamma_s)
                        found = self.generate(s, s+1, self.group0.u(s, gamma_s));
                        if found
                            mask = replab.bsgs.minimalMaskInOrbit(self.degree, self.KinBase0.S, self.baseOrdering);
                        end
                    end
                end
            end
        end

        function found = generate(self, s, i, prevG)
        % Generate the elements of ``G^s``
        %
        % Those elements have base image ``[gamma(1) ... gamma(i-1)] = [prevG(base0(1)) ... prevG(base0(i-1))]``
        % and they may have the required property.
        % If one is found then we extend ``K``, the subgroup of ``G^(s)`` of elements with property ``P`` that have
        % already been found, and return to search.
            if i == length(self.base0 + 1)
                found = self.prop(prevG);
                if found
                    self.KinBase0.stripAndAddStrongGenerator(prevG);
                    order = self.KinBase0.order;
                    self.KinBase0.randomizedSchreierSims([]);
                    assert(self.KinBase0.order == order);
                end
            else
                orbit = self.group0.Delta{i};
                for b = orbit
                    u = self.group0.u(i, b);
                    found = self.generate(s, i + 1, prevG(u));
                    if found
                        return
                    end
                end
            end
        end

        function l = greaterThan(self, x, y)
        % Point comparison using the base ordering
            l = self.baseOrdering(x) > self.baseOrdering(y);
        end

        function l = lessThan(self, x, y)
        % Point comparison using the base ordering
            l = self.baseOrdering(x) < self.baseOrdering(y);
        end

        function t = sort(self, s)
        % Sorts a sequence of points under the ordering
            [~, ind] = sort(self.baseOrdering(s));
            t = s(ind);
        end

        function ok = prop(self, g)
        % Tests the property that elements of the searched subgroup satisfy
            error('Abstract');
        end

        function ok = test(self, l, prev, ul)
        % Tests if base images are possible at a given level
        %
        % Note that the element ``g`` in Holt is given by ``g = compose(prev, ul)``.
        %
        % This tests if there is any element in the searched set that has the partial base image
        % ``[g(base(1)) ... g(base(l))]``. False negatives are possible, but not false positives.

        % Args:
        %   l (integer): Level of the test
        %   prev (permutation): Partial product ``u1 ... u{l-1}``
        %   ul (permutation): Current transversal element
        %
        % Returns:
        %   logical: False if there no element in the searched set that has the current partial base image
            error('Abstract');
        end

    end

end
