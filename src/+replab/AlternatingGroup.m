function grp = AlternatingGroup(n)
    Sn = replab.S(n);
    if n <= 2
        grp = replab.PermutationGroup.trivial(n);
    else
        t = [2 3 1 4:n];
        if n == 3
            grp = Sn.subgroup({t}, Sn.order/2);
            return
        end
        if mod(n, 2) == 0
            s = [2 1 4:n 3];
        else
            s = [1 2 4:n 3];
        end
        grp = Sn.subgroup({s t}, Sn.order/2);
    end
end
