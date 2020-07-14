classdef Backtrack < replab.Str

    properties

        degree % (integer): Permutation group degree
        group % (`+replab.PermutationGroup`): Group to search in
        initSubgroup

        baseOrdering % Base ordering including two guard elements

        base % Non-redundant base
        baseLen % local
        resBasicOrbits
        sortedOrbits
        orbits
        slowCosetTest
        tests
        testData
        groupedTests
        startData
        prop
        res
        g
        f
        l
        c
        mu
        nu
        minimalMaskInOrbit
    end

    methods

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

        function res = subgroup(self)
            % line 8: main loop
            while 1
                while self.l < self.baseLen
                    % line 10: apply all tests
                    img = self.g{self.l}(self.base(self.l));
                    if ~self.greaterThan(img, self.mu(self.l)) || ~self.lessThan(img, self.nu(self.l)) || ~self.minimalMaskInOrbit{self.l}(img)
                        break
                    end
                    ok = true;
                    data = self.testData{self.l};
                    seq = self.groupedTests{self.l};
                    for j = 1:length(seq)
                        [ok, data] = seq{j}(self.g{self.l}, data);
                        if ~ok
                            break
                        end
                    end
                    self.testData{self.l+1} = data;
                    if ~ok
                        break
                    end
                    if self.slowCosetTest
                        % line 11: change the (partial) base of K
                        self.res.baseChange([self.res.B(1:self.l-1) img]);
                        % line 12: calculate the minimal orbit representative mask
                        self.minimalMaskInOrbit{self.l+1} = replab.bsgs.minimalMaskInOrbit(self.degree, self.res.strongGeneratorsForLevel(self.l+1), self.baseOrdering);
                    else
                        self.minimalMaskInOrbit{self.l+1} = true(1, self.degree);
                    end
                    % line 13: recompute sorted orbits
                    self.l = self.l + 1;
                    self.sortedOrbits{self.l} = self.sort(self.g{self.l-1}(self.group.Delta{self.l}));
                    % lines 14 and 15: update variables used in minimality tests
                    self.computeMu(self.l);
                    self.computeNu(self.l);
                    % line 16: determine the new transversal element
                    self.c(self.l) = 1;
                    idx = self.sortedOrbits{self.l}(self.c(self.l));
                    gamma = find(self.g{self.l-1} == idx);
                    ul = self.group.u(self.l, gamma);
                    self.g{self.l} = self.g{self.l-1}(ul);
                end
                % lines 17: apply the tests to the group element found
                if self.l == self.baseLen
                    img = self.g{self.l}(self.base(self.l));
                    if self.minimalMaskInOrbit{self.l}(img) && self.greaterThan(img, self.mu(self.l)) && self.lessThan(img, self.nu(self.l))
                        ok = true;
                        data = self.testData{self.l};
                        seq = self.groupedTests{self.l};
                        for j = 1:length(seq)
                            [ok, data] = seq{j}(self.g{self.l}, data);
                            if ~ok
                                break
                            end
                        end
                        if ok && self.prop(self.g{self.l})
                            % line 18: add new strong generator for K
                            if self.slowCosetTest
                                % line 19-20: reset the base of K
                                self.res.baseChange(self.base);
                            end
                            self.res.stripAndAddStrongGenerator(self.g{self.l});
                            self.resBasicOrbits = self.res.Delta;
                            if self.slowCosetTest
                                % line 21: recalculate orbit representatives
                                self.minimalMaskInOrbit{self.f} = replab.bsgs.minimalMaskInOrbit(self.degree, self.res.strongGeneratorsForLevel(self.f), self.baseOrdering);
                            else
                                self.minimalMaskInOrbit{self.f} = true(1, degree);
                            end
                            % line 22: reset the search depth
                            self.l = self.f;
                        end
                    end
                end
                % line 23: go up the tree until in the first branch not fully seached
                while self.l > 0 && self.c(self.l) == self.group.orbitSize(self.l)
                    self.l = self.l - 1;
                end
                % line 24: if the entire tree is traversed, return K
                if self.l == 0
                    self.res.makeImmutable;
                    res = self.res;
                    return
                end
                % line 25-27: update orbit representatives
                if self.l < self.f
                    % line 26
                    self.f = self.l;
                    self.c(self.l) = 1;
                    if self.slowCosetTest
                        % line 27
                        self.minimalMaskInOrbit{self.f} = replab.bsgs.minimalMaskInOrbit(self.degree, self.res.strongGeneratorsForLevel(self.f), self.baseOrdering);
                    else
                        self.minimalMaskInOrbit{f} = true(1, self.degree);
                    end
                    % line 28: update variables used for minimality testing
                    self.mu(self.l) = self.degree + 2; % = 0
                    self.computeNu(self.l);
                end
                % line 29: set the next element from the current branch and update accordingly
                self.c(self.l) = self.c(self.l) + 1;
                idx = self.sortedOrbits{self.l}(self.c(self.l));
                if self.l == 1
                    gamma = idx;
                else
                    gamma = find(self.g{self.l-1} == idx);
                end
                ul = self.group.u(self.l, gamma);
                if self.l == 1
                    self.g{self.l} = ul;
                else
                    self.g{self.l} = self.g{self.l-1}(ul);
                end
            end
        end

        function self = Backtrack(group, prop, tests, startData, initSubgroup, slowCosetTest)
            degree = group.n;
            if nargin < 6 || isempty(slowCosetTest)
                self.slowCosetTest = true;
            else
                self.slowCosetTest = slowCosetTest;
            end
            if nargin < 5 || isempty(initSubgroup)
                self.initSubgroup = replab.bsgs.Chain(degree);
            else
                self.initSubgroup = initSubgroup;
            end
            if nargin < 4 || isempty(tests)
                self.tests = {};
                self.startData = [];
            else
                self.tests = tests;
                self.startData = startData;
            end
            self.prop = prop;
            self.baseOrdering = [replab.bsgs.baseOrdering(degree, group.base) degree+1 0];
            self.group = group;
            self.degree = group.n;
            [self.group, self.groupedTests, self.startData] = replab.bsgs.cleanUpBaseAndTests(self.group, self.tests, self.startData);
            self.base = self.group.base;
            self.baseLen = length(self.base);
            if self.baseLen == 0
                self.res = replab.bsgs.Chain(self.degree);
                self.res.makeImmutable;
                return
            end
            identity = 1:self.degree;
            self.testData = cell(1, self.baseLen);
            self.testData{1} = self.startData;
            for i = 1:self.baseLen
                seq = self.groupedTests{i};
                data = self.testData{i};
                for j = 1:length(seq)
                    [ok, data] = seq{j}(identity, data);
                end
                self.testData{i+1} = data;
            end
            % line 1: more initializations
            self.res = self.initSubgroup.mutableCopy;
            self.f = self.baseLen;
            self.l = self.baseLen;
            % line 2: set the base for K to the base of G
            % line 3: compute BSGS and related structure for K
            self.res.baseChange(self.base);
            self.resBasicOrbits = self.res.Delta; % Delta_K
                                                  % line 4: orbit representatives for f-th basic stabilizer of K
                                                  % instead of storing the orbit representatives and using ismember, we store
                                                  % a logical mask which is true if the element is minimal
            self.minimalMaskInOrbit{self.f} = replab.bsgs.minimalMaskInOrbit(self.degree, self.res.strongGeneratorsForLevel(self.f), self.baseOrdering);
            % line 5: remove the base point from the representatives to avoid getting the identity element as a generator for K
            self.minimalMaskInOrbit{self.f}(self.base(self.f)) = false;
            % line 6: more initializations
            self.c = zeros(1, self.baseLen);
            self.sortedOrbits = cell(1, self.baseLen); % = \Lambda
            for i = 1:self.baseLen
                self.sortedOrbits{i} = self.sort(self.group.Delta{i});
            end
            % line 7: initializations
            self.mu = zeros(1, self.baseLen);
            self.nu = zeros(1, self.baseLen);
            % this corresponds to the element smaller than all points
            self.mu(self.l) = self.degree + 2;
            self.computeNu(self.l);
            % initialized computed words
            self.g = repmat({identity}, 1, self.baseLen);
        end

        function computeMu(self, l)
        % Computes the ``mu`` bound in the backtracking tests
            mu_l = self.degree + 2; % place holder for element < all others in base ordering
            for j = 1:l
                if any(self.base(l) == self.resBasicOrbits{j})
                    candidate = self.g{j}(self.base(j));
                    if self.baseOrdering(candidate) > mu_l
                        mu_l = candidate;
                    end
                end
            end
            self.mu(l) = mu_l;
        end

        function nu = computeNu(self, l);
            idx = length(self.group.Delta{l}) + 2 - length(self.resBasicOrbits{l});
            if idx > length(self.sortedOrbits{l})
                self.nu(l) = self.degree + 1; % place holder for element > all others in base ordering
            else
                self.nu(l) = self.sortedOrbits{l}(idx);
            end
        end

    end

end
