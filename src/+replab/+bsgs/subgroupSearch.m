function res = subgroupSearch(group, prop, base, tests, init_subgroup)
% Find the subgroup of all elements satisfying the property ``prop``
%
% Code lifted from Sympy 1.6, PermutationGroup.subgroup_search
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
%   init_subgroup (`+replab.+bsgs.Chain`): BSGS chain representing a subgroup, if a subgroup of the sought group is
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
%   >>> isEven = @(g) replab.bsgs.permsign(g) == 1;
%   >>> H = replab.PermutationGroup.fromChain(replab.bsgs.subgroupSearch(Sn.chain, isEven));
%   >>> Sn.order/H.order
%        2
%
%
% Notes:
%   This function is extremely lengthy and complicated and will require
%   some careful attention. The implementation is described in
%   [1], pp. 114-117, and the comments for the code here follow the lines
%   of the pseudocode in the book for clarity.
%   The complexity is exponential in general, since the search process by
%   itself visits all members of the supergroup. However, there are a lot
%   of tests which are used to prune the search tree, and users can define
%   their own tests via the ``tests`` parameter, so in practice, and for
%   some computations, it's not terrible.
%   A crucial part in the procedure is the frequent base change performed
%   (this is line 11 in the pseudocode) in order to obtain a new basic
%   stabilizer. The book mentiones that this can be done by using a base swap,
%   however the current implementation uses the `+replab.+bsgs.Chain.stabilizer` method
%   on the previous basic stabilizer.
    degree = group.n;
    % initialize basic properties
    if nargin < 5
        init_subgroup = replab.bsgs.Chain(degree);
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
    base_len = length(base);
    % handle tests: if tests do not cover base fully, we complete by trivial tests
    for i = length(tests)+1:base_len
        tests{1,i} = @(x) true;
    end
    identity = 1:degree;
    % add elements larger and smaller than all points
    base_ordering = [base_ordering_(base, degree) degree+1 0];
    % line 1: more initializations
    res = init_subgroup.mutableCopy;
    f = base_len; % IDX
    l = base_len; % IDX
    % line 2: set the base for K to the base of G
    res.baseChange(base);
    % line 3: compute BSGS and related structure for K
    % res_basic_orbits_init_base == res.Delta
    % initialize orbit representatives
    orbit_reps = cell(1, base_len);
    % line 4: orbit representatives for f-th basic stabilizer of K
    orbits = orbits_(degree, res.strongGeneratorsForLevel(f));
    orbit_reps{f} = get_reps(orbits, base_ordering);
    % line 5: remove the base point from the representatives to avoid getting the identity element as a generator for K
    orbit_reps{f} = setdiff(orbit_reps{f}, base(f));
    % line 6: more initializations
    c = zeros(1, base_len);
    u = repmat({identity}, 1, base_len);
    sorted_orbits = cell(1, base_len);
    for i = 1:base_len
        sorted_orbits{i} = group.Delta{i};
        [~, I] = sort(base_ordering(sorted_orbits{i}));
        sorted_orbits{i} = sorted_orbits{i}(I);
    end
    % line 7: initializations
    mu = zeros(1, base_len);
    nu = zeros(1, base_len);
    % this corresponds to the element smaller than all points
    mu(l) = degree + 2; % IDX
    nu = update_nu(nu, l, group.Delta, res.Delta, sorted_orbits, degree, base_ordering);
    % initialized computed words
    computed_words = repmat({identity}, 1, base_len);
    % line 8: main loop
    while 1
        % apply all the tests
         while l < base_len && ...
                 base_ordering(computed_words{l}(base(l))) < base_ordering(nu(l)) && ...
                 ismember(computed_words{l}(base(l)), orbit_reps{l}) && ...
                 tests{l}(computed_words{l})
             %base_ordering(mu(l)) < base_ordering(computed_words{l}(base(l))) && ...
            % change the (partial) base of K
            res.baseChange([res.B(1:l-1) computed_words{l}(base(l))]);
            orbits = orbits_(degree, res.strongGeneratorsForLevel(l+1));
            orbit_reps{l+1} = get_reps(orbits, base_ordering);
            % line 13: amend sorted orbits
            l = l + 1;
            temp_orbit = computed_words{l-1}(group.Delta{l});
            [~, I] = sort(base_ordering(temp_orbit));
            temp_orbit = temp_orbit(I);
            sorted_orbits{l} = temp_orbit;
            % lines 14 and 15: update variables used in minimality tests
            new_mu = degree + 2; % IDX position of the 0 canary
            for i = 1:l
                if res.iDelta(base(l), i) ~= 0
                    candidate = computed_words{i}(base(i));
                    if base_ordering(candidate) > base_ordering(new_mu)
                        new_mu = candidate;
                    end
                end
            end
            mu(l) = new_mu;
            nu = update_nu(nu, l, group.Delta, res.Delta, sorted_orbits, degree, base_ordering);
            % line 16: determine the new transversal element
            c(l) = 1; % IDX
            temp_point = sorted_orbits{l}(c(l));
            gamma = find(computed_words{l-1} == temp_point);
            ul = group.u(l, gamma);
            computed_words{l} = computed_words{l-1}(ul);
        end
        % lines 17 & 18: apply the tests to the group element found
        g = computed_words{l};
        temp_point = g(base(l));
        if l == base_len && ... % IDX
              base_ordering(temp_point) < base_ordering(nu(l)) && ...
              ismember(temp_point, orbit_reps{l}) && ...
              tests{l}(computed_words{l}) && ...
              prop(g)
            % base_ordering(mu(l)) < base_ordering(temp_point) && ...
            % line 18: add new strong generator for K
            % line 19-20: reset the base of K
            res.baseChange(base);
            res.stripAndAddStrongGenerator(g);
            % line 21: recalculate orbit representatives
            orbits = orbits_(degree, res.strongGeneratorsForLevel(f));
            orbit_reps{f} = get_reps(orbits, base_ordering);
            % line 22: reset the search depth
            l = f;
        end
        % line 23: go up the tree until in the first branch not fully seached
        while l > 0 && c(l) == group.orbitSize(l) % IDX IDX
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
            c(l) = 1; % IDX
            % line 27
            temp_orbits = orbits_(degree, res.strongGeneratorsForLevel(f));
            orbit_reps{f} = get_reps(temp_orbits, base_ordering);
            % line 28: update variables used for minimality testing
            mu(l) = degree + 2; % position of the minimal one
            temp_index = group.orbitSize(l) + 2 - res.orbitSize(l);
            if temp_index > length(sorted_orbits{l}) % IDX
                nu(l) = base_ordering(degree + 1);
            else
                nu(l) = sorted_orbits{l}(temp_index);
            end
        end
        % line 29: set the next element from the current branch and update accordingly
        c(l) = c(l) + 1;
        if l == 1
            gamma = sorted_orbits{l}(c(l));
        else
            gamma = find(computed_words{l-1} == sorted_orbits{l}(c(l)));
        end
        ul = group.u(l, gamma);
        if l == 1
            computed_words{l} = ul;
        else
            computed_words{l} = computed_words{l-1}(ul);
        end
    end
end

function ordering = base_ordering_(base, degree)
    base_len = length(base);
    ordering = zeros(1, degree);
    ordering(base) = 1:base_len;
    rest = setdiff(1:degree, base);
    ordering(rest) = base_len + (1:length(rest));
end

function res = get_reps(orbits, base_ordering)
% get the minimal element in the base ordering
    res = zeros(1, length(orbits));
    for i = 1:length(orbits)
        orbit = orbits{i};
        [~,ind] = min(base_ordering(orbit));
        res(i) = orbit(ind);
    end
end

function nu = update_nu(nu, l, basic_orbits, res_basic_orbits_init_base, sorted_orbits, degree, base_ordering)
    temp_index = length(basic_orbits{l}) + 2 - length(res_basic_orbits_init_base{l}); % IDX
    if temp_index > length(sorted_orbits{l})
        nu(l) = base_ordering(degree + 1); % IDX
    else
        nu(l) = sorted_orbits{l}(temp_index);
    end
end

function orbs = orbits_(degree, generators)
    if size(generators, 2) == 0
        generators = zeros(degree, 0);
    end
    orbs = replab.Partition.permutationsOrbits(generators').blocks;
end