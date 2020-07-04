function U = randomSubgroup(G)
% Returns a random subgroup of G
%
% Does not provide any guarantee on uniformity, etc.
%
% Note that the subgroup of Sn obtained in this way is often the alternating group
%
% Args:
%   G (`+replab.NiceFiniteGroup`): Group
%
% Returns:
%   `+replab.NiceFiniteGroup`: A random subgroup of G
    assert(~isprime(G.order), 'A cyclic group of prime order has no proper subgroup');
    U = G;
    while U.order == G.order
        s1 = G.sample;
        while G.isIdentity(s1)
            s1 = G.sample;
        end
        s2 = G.sample;
        while G.isIdentity(s2)
            s2 = G.sample;
        end
        U = G.subgroup({s1 s2});
    end
end
