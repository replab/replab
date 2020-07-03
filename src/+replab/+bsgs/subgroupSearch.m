function res= subgroupSearch(group, prop, base, tests, initSubgroup)
% Find the subgroup of all elements satisfying the property ``prop``
%
% See also a similar implementation in Sympy 1.6, PermutationGroup.subgroup_search
%
% Args:
%   group (`+replab.+bsgs.Chain`): BSGS chain representing the group
%   prop (function_handle): The property to be used. Has to be callable on group elements
%                           and return a logical value. It is assumed that all group elements
%                           satisfying ``prop`` indeed form a subgroup.
%   base (integer(1,\*), optional): A (potentially partial) base for the supergroup, used in conjonction with the tests
%   tests (cell(1,\*) of function_handle, optional): A list of function handles. These are used to rule ou group elements
%                                                    by partial base images, so that ``tests{l}(g)`` returns false if the
%                                                    elements ``g`` is known not to satisfy ``prop`` on where ``g`` sends the
%                                                    first ``l`` base points.
%   initSubgroup (`+replab.+bsgs.Chain`): BSGS chain representing a subgroup, if a subgroup of the sought group is
%                                         known in advance, it can be passed to the function as this parameter.
%
% Returns:
%   `+replab.+bsgs.Chain`: The BSGS chain representing the subgroup of all elements satisfying ``prop``. The strong
%                          generating set of this chain is guaranteed to be a strong generating set relative to the
%                          base ``base``.
%
% Example:
%   >>> n = 12;
%   >>> Sn = replab.S(n);
%   >>> isEven = @(g) Sn.sign(g) == 1;
%   >>> H = replab.PermutationGroup.fromChain(replab.bsgs.subgroupSearch(Sn.chain, isEven));
%   >>> Sn.order/H.order
%        2
%
    degree = group.n;
    % initialize basic properties
    if nargin < 5
        initSubgroup = replab.bsgs.Chain(degree);
    end
    if nargin < 4
        tests = {};
        base = [];
    end
    if ~isempty(base)
        % change the BSGS base to the base provided
        group = group.mutableCopy;
        group.baseChange(base);
        group.makeImmutable;
        % verify we start with the partial base provided
        assert(isequal(group.B(1:length(base)), base));
    end
    % but the actual base may have additional points at the end
    base = group.B;
    baseLen = length(base);
    % handle tests: if tests do not cover base fully, we complete by trivial tests
    for i = length(tests)+1:baseLen
        tests{1,i} = @(x) true;
    end
    identity = 1:degree;
    baseOrdering = [replab.bsgs.baseOrdering(degree, base) degree+1 0];
    % line 1: more initializations
    res = initSubgroup.mutableCopy;
    f = baseLen; % IDX
    l = baseLen; % IDX
                 % line 2: set the base for K to the base of G
                 % line 3: compute BSGS and related structure for K
    res.baseChange(base);
    resBasicOrbits = res.Delta; % Delta_K
                                % line 4: orbit representatives for f-th basic stabilizer of K
                                % instead of storing the orbit representatives and using ismember, we store
                                % a logical mask which is true if the element is minimal
    minimalInOrbit{f} = replab.bsgs.minimalByBaseOrderingInOrbit(degree, res.strongGeneratorsForLevel(f), baseOrdering);
    % line 5: remove the base point from the representatives to avoid getting the identity element as a generator for K
    minimalInOrbit{f}(base(f)) = false;
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
    nu(l) = computeNu(degree, l, sortedOrbits, group.Delta, resBasicOrbits);
    % initialized computed words
    g = repmat({identity}, 1, baseLen);
    greaterThan = @(x, y) baseOrdering(x) > baseOrdering(y);
    lessThan = @(x, y) baseOrdering(x) < baseOrdering(y);
    % line 8: main loop
    while 1
        while l < baseLen
            % line 10: apply all tests
            img = g{l}(base(l));
            if ~greaterThan(img, mu(l)) || ~lessThan(img, nu(l)) || ~minimalInOrbit{l}(img) || ~tests{l}(g{l})
                break
            end
            % line 11: change the (partial) base of K
            res.baseChange([res.B(1:l-1) img]);
            % line 12: calculate the minimal orbit representative mask
            minimalInOrbit{l+1} = replab.bsgs.minimalByBaseOrderingInOrbit(...
                degree, res.strongGeneratorsForLevel(l + 1), baseOrdering);
            % line 13: recompute sorted orbits
            l = l + 1;
            sortedOrbits{l} = replab.bsgs.sortByOrdering(g{l-1}(group.Delta{l}), baseOrdering);
            % lines 14 and 15: update variables used in minimality tests
            mu(l) = computeMu(degree, l, base, g, resBasicOrbits, baseOrdering);
            nu(l) = computeNu(degree, l, sortedOrbits, group.Delta, resBasicOrbits);
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
            if minimalInOrbit{l}(img) && greaterThan(img, mu(l)) && lessThan(img, nu(l)) && tests{l}(g{l}) && prop(g{l})
                % line 18: add new strong generator for K
                % line 19-20: reset the base of K
                res.baseChange(base);
                res.stripAndAddStrongGenerator(g{l});
                resBasicOrbits = res.Delta;
                % line 21: recalculate orbit representatives
                minimalInOrbit{f} = replab.bsgs.minimalByBaseOrderingInOrbit(degree, res.strongGeneratorsForLevel(f), baseOrdering);
                % line 22: reset the search depth
                l = f;
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
            minimalInOrbit{f} = replab.bsgs.minimalByBaseOrderingInOrbit(degree, res.strongGeneratorsForLevel(f), baseOrdering);
            % line 28: update variables used for minimality testing
            mu(l) = degree + 2; % = 0
            nu(l) = computeNu(degree, l, sortedOrbits, group.Delta, resBasicOrbits);
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

function mu = computeMu(degree, l, base, g, resBasicOrbits, baseOrdering)
    mu = degree + 2; % place holder for element < all others in base ordering
    for j = 1:l
        if ismember(base(l), resBasicOrbits{j})
            candidate = g{j}(base(j));
            if baseOrdering(candidate) > mu
                mu = candidate;
            end
        end
    end
end

function nu = computeNu(degree, l, sortedOrbits, orbits, resBasicOrbits)
    idx = length(orbits{l}) + 2 - length(resBasicOrbits{l});
    if idx > length(sortedOrbits{l})
        nu = degree + 1; % place holder for element > all others in base ordering
    else
        nu = sortedOrbits{l}(idx);
    end
end
