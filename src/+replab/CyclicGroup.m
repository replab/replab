function grp = CyclicGroup(n)
    Sn = replab.S(n);
    if n == 1
        grp = Sn;
    else
        grp = Sn.subgroup({[2:n 1]}, vpi(n));
    end
end
