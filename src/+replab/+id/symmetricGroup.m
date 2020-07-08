function [pres s t] = symmetricGroup(G)
% Recognizes if the given group is the symmetric group
%
% Args:
%   G (`+replab.PermutationGroup`): Permutation group to test for isomorphism to a symmetric group
    pres = [];
    s = [];
    t = [];
    [n r] = replab.id.unfactorial(G.order);
    n = double(n);
    if r ~= 0
        return
    end
    if n == G.domainSize && false
        s = [2:n 1];
        t = [2 1 3:n];
    else
        s = replab.id.findElementOfOrder(G, n);
        if isempty(s)
            return
        end
        t = replab.id.findElementOfOrder(G, 2, @(g) testPresentation(G, n, s, g));
        if isempty(t)
            return
        end
    end
end

function ok = testPresentation(G, n, s, t)
% Assumes that s^n = t^2 = id
    ok = false;
    if ~G.isIdentity(G.composeN(s(t), n-1))
        return
    end
    sj = s; % for j = 1
    for j = 2:floor(n/2)
        sj = sj(s); % now sj = s^j
        % compute t sj t^-1 sj^-1 = t sj (sj t)^-1
        comm = G.compose(t(sj), G.inverse(sj(t)));
        comm2 = comm(comm);
        if ~G.isIdentity(comm2)
            return
        end
    end
    ok = true;
end
