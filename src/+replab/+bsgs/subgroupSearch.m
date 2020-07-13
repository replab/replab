function res = subgroupSearch(group, prop, tests, startData, initSubgroup)
% Find the subgroup of all elements satisfying the property ``prop``
%
% See also a similar implementation in Sympy 1.6, PermutationGroup.subgroup_search
%
% It is assumed that the set of elements of ``group`` that satisfy ``prop`` form a group. We call
% that subgroup ``K``.
%
% The tests are given by function handles in the cell array ``tests``. They refer to the base ``base`` of the
% stabilizer chain ``group`` defining the group.
% The tests are used to rule out group elements by partial base images.
%
% We consider the ``l``-th level, and the test ``tests{l}``. This test has the calling convention
%
% ``[ok, outdata] = tests{l}(g, indata)``.
%
% The return value ``ok`` is logical, and can be false only when the subset ``P`` of ``K`` with the
% same partial base image as ``g`` is empty:
%
% ``P = {k \in K : k(base(i)) = g(base(i)) for all i=1,...,l}`` is empty.
%
% In the simple case, ``indata`` and ``outdata`` are always set to ``[]``. Otherwise, we use the following convention.
%
% Note that the tested elements ``g = u{1} u{2} ... u{l}`` where ``u{1} ... u{l}`` are transversal elements from the stabilizer chain ``group``.
%
% - At level ``l = 1``, ``indata{1}`` is given as a parameter, and we store the result ``outdata{2}`` from ``[ok, outdata{2}] = tests{1}(u{1}, indata{1})``.
% - At level ``l = 2``, we call ``[ok, outdata{3}] = tests{2}(compose(u{1}, u{2}), outdata{2})``.
% - At level ``l``, we call ``[ok, ~] = tests{l}(compose(u{1}, ..., u{l-1}, u{l}), outdata{l})``.
%
% Args:
%   group (`+replab.+bsgs.Chain`): BSGS chain representing the group
%   prop (function_handle): Tested property. Takes a permutation as an argument and returns a logical value.
%   tests (cell(1,\*) of function_handle, optional): A list of function handles, see description above.
%   startData (user-defined): Data to pass to the first test function
%   initSubgroup (`+replab.+bsgs.Chain`): BSGS chain representing a subgroup, if a subgroup of the sought group is
%                                         known in advance, it can be passed to the function as this parameter.
%
% Returns:
%   `+replab.+bsgs.Chain`: The BSGS chain representing the subgroup of all elements satisfying ``prop``. The strong
%                          generating set of this chain is guaranteed to be a strong generating set relative to the
%                          base ``base``.
%
%
% Example:
%   >>> n = 12;
%   >>> Sn = replab.S(n);
%   >>> isEven = @(g) replab.Permutation.sign(g) == 1;
%   >>> H = replab.PermutationGroup.fromChain(replab.bsgs.subgroupSearch(Sn.chain, isEven));
%   >>> Sn.order/H.order
%        2
%
% Note:
%   This function is tricky to get right, and comes from the pseudo code in pp. 114-117 of
%   D. Holt et al, Handbook of Computational Group Theory (CRC Press, 2005).
%   As in the Sympy implementation, we annotate the code according to the lines
%   of the pseudocode.
    degree = group.n;
    % initialize basic properties
    if nargin < 5
        initSubgroup = replab.bsgs.Chain(degree);
    end
    if nargin < 4 || isempty(tests)
        tests = {};
        startData = [];
    end
    [group, groupedTests, startData] = replab.bsgs.cleanUpBaseAndTests(group, tests, startData);
    tests = [];
    base = group.base;
    baseLen = length(base);
    if baseLen == 0
        res = replab.bsgs.Chain(degree);
        res.makeImmutable;
        return
    end
    identity = 1:degree;
    testData = cell(1, baseLen);
    testData{1} = startData;
    for i = 1:baseLen
        seq = groupedTests{i};
        data = testData{i};
        for j = 1:length(seq)
            [ok, data] = seq{j}(identity, data);
        end
        testData{i+1} = data;
    end
    baseOrdering = [replab.bsgs.baseOrdering(degree, base) degree+1 0];
    % line 1: more initializations
    res = initSubgroup.mutableCopy;
    f = baseLen;
    l = baseLen;
    % line 2: set the base for K to the base of G
    % line 3: compute BSGS and related structure for K
    res.baseChange(base);
    resBasicOrbits = res.Delta; % Delta_K
                                % line 4: orbit representatives for f-th basic stabilizer of K
                                % instead of storing the orbit representatives and using ismember, we store
                                % a logical mask which is true if the element is minimal
    minimalMaskInOrbit{f} = replab.bsgs.minimalMaskInOrbit(degree, res.strongGeneratorsForLevel(f), baseOrdering);
    % line 5: remove the base point from the representatives to avoid getting the identity element as a generator for K
    minimalMaskInOrbit{f}(base(f)) = false;
    % line 6: more initializations
    c = zeros(1, baseLen);
    u = repmat({identity}, 1, baseLen);
    sortedOrbits = cell(1, baseLen); % = \Lambda
    for i = 1:baseLen
        sortedOrbits{i} = replab.bsgs.sortByOrdering(group.Delta{i}, baseOrdering);
    end
    % line 7: initializations
    mu = zeros(1, baseLen);
    nu = zeros(1, baseLen);
    % this corresponds to the element smaller than all points
    mu(l) = degree + 2;
    nu(l) = replab.bsgs.computeNu(degree, l, sortedOrbits, group.Delta, resBasicOrbits);
    % initialized computed words
    g = repmat({identity}, 1, baseLen);
    greaterThan = @(x, y) baseOrdering(x) > baseOrdering(y);
    lessThan = @(x, y) baseOrdering(x) < baseOrdering(y);
    % line 8: main loop
    while 1
        while l < baseLen
            % line 10: apply all tests
            img = g{l}(base(l));
            if ~greaterThan(img, mu(l)) || ~lessThan(img, nu(l)) || ~minimalMaskInOrbit{l}(img)
                break
            end
            ok = true;
            data = testData{l};
            seq = groupedTests{l};
            for j = 1:length(seq)
                [ok, data] = seq{j}(g{l}, data);
                if ~ok
                    break
                end
            end
            testData{l+1} = data;
            if ~ok
                break
            end
            % line 11: change the (partial) base of K
            res.baseChange([res.B(1:l-1) img]);
            % line 12: calculate the minimal orbit representative mask
            minimalMaskInOrbit{l+1} = replab.bsgs.minimalMaskInOrbit(degree, res.strongGeneratorsForLevel(l+1), baseOrdering);
            % line 13: recompute sorted orbits
            l = l + 1;
            sortedOrbits{l} = replab.bsgs.sortByOrdering(g{l-1}(group.Delta{l}), baseOrdering);
            % lines 14 and 15: update variables used in minimality tests
            mu(l) = replab.bsgs.computeMu(degree, l, base, g, resBasicOrbits, baseOrdering);
            nu(l) = replab.bsgs.computeNu(degree, l, sortedOrbits, group.Delta, resBasicOrbits);
            % line 16: determine the new transversal element
            c(l) = 1;
            idx = sortedOrbits{l}(c(l));
            gamma = find(g{l-1} == idx);
            ul = group.u(l, gamma);
            g{l} = g{l-1}(ul);
        end
        % lines 17: apply the tests to the group element found
        if l == baseLen
            img = g{l}(base(l));
            if minimalMaskInOrbit{l}(img) && greaterThan(img, mu(l)) && lessThan(img, nu(l))
                ok = true;
                data = testData{l};
                seq = groupedTests{l};
                for j = 1:length(seq)
                    [ok, data] = seq{j}(g{l}, data);
                    if ~ok
                        break
                    end
                end
                if ok && prop(g{l})
                    % line 18: add new strong generator for K
                    % line 19-20: reset the base of K
                    res.baseChange(base);
                    res.stripAndAddStrongGenerator(g{l});
                    resBasicOrbits = res.Delta;
                    % line 21: recalculate orbit representatives
                    minimalMaskInOrbit{f} = replab.bsgs.minimalMaskInOrbit(degree, res.strongGeneratorsForLevel(f), baseOrdering);
                    % line 22: reset the search depth
                    l = f;
                end
            end
        end
        % line 23: go up the tree until in the first branch not fully seached
        while l > 0 && c(l) == group.orbitSize(l)
            l = l - 1;
        end
        % line 24: if the entire tree is traversed, return K
        if l == 0
            res.makeImmutable;
            return
        end
        % line 25-27: update orbit representatives
        if l < f
            % line 26
            f = l;
            c(l) = 1;
            % line 27
            minimalMaskInOrbit{f} = replab.bsgs.minimalMaskInOrbit(degree, res.strongGeneratorsForLevel(f), baseOrdering);
            % line 28: update variables used for minimality testing
            mu(l) = degree + 2; % = 0
            nu(l) = replab.bsgs.computeNu(degree, l, sortedOrbits, group.Delta, resBasicOrbits);
        end
        % line 29: set the next element from the current branch and update accordingly
        c(l) = c(l) + 1;
        idx = sortedOrbits{l}(c(l));
        if l == 1
            gamma = idx;
        else
            gamma = find(g{l-1} == idx);
        end
        ul = group.u(l, gamma);
        if l == 1
            g{l} = ul;
        else
            g{l} = g{l-1}(ul);
        end
    end
end
