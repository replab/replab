function res = backtrackSearch(group, prop, tests, startData, leftSubgroup, rightSubgroup)
% Find an element in group that satisfies a property
%
% We accelerate the search by the following assumption: for any element ``g`` of ``group`` that
% satisfies the property, then all elements of the coset ``leftSubgroup g rightSubgroup``
% also satisfy the property. Note that one or two of these subgroups may be trivial.
%
% This implementation is based heavily on `+replab.+bsgs.subgroupSearch`.
%
% Args:
%   group (`+replab.+bsgs.Chain`): BSGS chain representing the group
%   prop (function_handle): The property to be used. Has to be callable on group elements
%                           and return a logical value.
%   base (integer(1,\*), optional): A (potentially partial) base for the group, used in conjonction with the tests
%   tests (cell(1,\*) of function_handle, optional): A list of function handles, see description in `.subgroupSearch`
%   startData (user-defined): Data to pass to the first test function
%   leftSubgroup (`+replab.+bsgs.Chain`): BSGS chain representing the left subgroup in the equivalence assumption above
%   rightSubgroup (`+replab.+bsgs.Chain`): BSGS chain representing the right subgroup in the equivalence assumption above
%
% Returns:
%   permutation: An element satisfying the property if it exists
    degree = group.n;
    % initialize basic properties
    if nargin < 6 || isempty(rightSubgroup)
        rightSubgroup = replab.bsgs.Chain(degree);
    end
    if nargin < 5 || isempty(leftSubgroup)
        leftSubgroup = replab.bsgs.Chain(degree);
    end
    if nargin < 4
        tests = {};
        startData = [];
    end
    identity = 1:degree;
    [group, tests, startData] = replab.bsgs.cleanUpBaseAndTests(group, tests, startData);
    base = group.base;
    baseLen = length(base);
    if baseLen == 0
        if prop(identity)
            res = identity;
        else
            res = [];
        end
        return
    end
    testData = cell(1, baseLen);
    testData{1} = startData;
    for i = 1:baseLen
        [ok, testData{i+1}] = tests{i}(identity, testData{i});
    end
    baseOrdering = [replab.bsgs.baseOrdering(degree, base) degree+1 0];
    % line 1-2: initialization
    % in the subgroup search algorithm, we construct a subgroup K by adding the new strong
    % generators found, and the coset minimality tests were performed on K by using it both
    % as a left and right subgroup. Here, we have two different subgroups playing that role.
    left = leftSubgroup.mutableCopy;
    left.baseChange(base);
    right = rightSubgroup.mutableCopy;
    right.baseChange(base);
    f = baseLen;
    l = baseLen;
    % line 3: compute BSGS and related structure for K
    minimalMaskInOrbit{f} = replab.bsgs.minimalMaskInOrbit(degree, right.strongGeneratorsForLevel(f), baseOrdering);
    % line 5: DO NOT remove the base point from the representatives, we want to test the identity!
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
    nu(l) = replab.bsgs.computeNu(degree, l, sortedOrbits, group.Delta, left.Delta);
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
            [ok, testData{l+1}] = tests{l}(g{l}, testData{l});
            if ~ok
                break
            end
            % line 11: change the (partial) base
            right.baseChange([right.B(1:l-1) img]);
            % line 12: calculate the minimal orbit representative mask
            minimalMaskInOrbit{l+1} = replab.bsgs.minimalMaskInOrbit(degree, right.strongGeneratorsForLevel(l+1), baseOrdering);
            % line 13: recompute sorted orbits
            l = l + 1;
            sortedOrbits{l} = replab.bsgs.sortByOrdering(g{l-1}(group.Delta{l}), baseOrdering);
            % lines 14 and 15: update variables used in minimality tests
            mu(l) = replab.bsgs.computeMu(degree, l, base, g, left.Delta, baseOrdering);
            nu(l) = replab.bsgs.computeNu(degree, l, sortedOrbits, group.Delta, left.Delta);
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
                [ok, ~] = tests{l}(g{l}, testData{l});
                if ok && prop(g{l})
                    res = g{l};
                    return
                end
                % skip lines 18-22 as we are not modifying the left and right subgroups
            end
        end
        % line 23: go up the tree until in the first branch not fully seached
        while l > 0 && c(l) == group.orbitSize(l)
            l = l - 1;
        end
        % line 24: if the entire tree is traversed, return ``[]`` as we are not successful
        if l == 0
            res = [];
            return
        end
        % line 25-27: update orbit representatives
        if l < f
            % line 26
            f = l;
            c(l) = 1;
            % line 27
            minimalMaskInOrbit{f} = replab.bsgs.minimalMaskInOrbit(degree, right.strongGeneratorsForLevel(f), baseOrdering);
            % line 28: update variables used for minimality testing
            mu(l) = degree + 2; % = 0
            nu(l) = replab.bsgs.computeNu(degree, l, sortedOrbits, group.Delta, left.Delta);
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
