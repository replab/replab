function res = equality(lhs, rhs)
    res = false;
    if ~isa(lhs, 'replab.FiniteSet') || ~isa(rhs, 'replab.FiniteSet')
        return
    end
    if ~lhs.hasSameTypeAs(rhs)
        return
    end
    if lhs.nElements ~= rhs.nElements
        return
    end
    c = [replab.FiniteSet.shortType(lhs) replab.FiniteSet.shortType(rhs)];
    type = lhs.type;
    switch c
      case 'GG'
        res = lhs.isSubgroupOf(rhs) && rhs.isSubgroupOf(lhs);
        return
        % the following handles all combinations of cosets
      case {'LL', 'LN', 'NL', 'NN', 'RR', 'RN', 'NR'}
        % g1 H1 == g2 H2
        % g2^-1 g1 H1 == H2
        % if g2^-1 g1 is not in H1, then H1 is not a group and the equality test fails
        % then we must have H1 == H2
        % and then both cosets must have the same canonical representative
        res = (lhs.subgroup == rhs.subgroup) && type.eqv(lhs.representative, rhs.representative);
        return
      case 'CC'
        if lhs.group == rhs.group
            res = type.eqv(lhs.representative, rhs.representative);
            return
        end
    end
    replab.msg(1, 'FiniteSet.eq: using generic equality test, which can be slow')
    E1 = lhs.elements;
    E2 = rhs.elements;
    for i = 1:length(E1)
        if ~type.eqv(E1{i}, E2{i})
            return
        end
    end
    res = true;
end
