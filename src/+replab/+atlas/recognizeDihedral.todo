function R = recognizeDihedral(G)
% Recognizes if the given group is the dihedral group and provides the generators according to the standard presentation
%
% The standard presentation is ``<x, a| a^n = x^2 = id, x a x^-1 = a>``
    assert(G.order > 2);
    R = [];
    x = [];
    a = [];
    if G.order > 2^53-1
        return
    end
    if mod(G.order, 2) == 1
        return
    end
    n = double(G.order/2);
    if mod(n, 2) == 1
        Zn = G.derivedSubgroup;
        if ~(Zn.isCyclic && Zn.order == n)
            return
        end
    else
        G1 = G.derivedSubgroup;
        if ~(G1.isCyclic && G1.order == n/2)
            return
        end
        T = G.rightCosetsOf(G1).transversal;
        good = false;
        for i = 1:4
            Zn = G1.closure(T{i});
            if Zn.isCyclic && Zn.order == n
                good = true;
                break
            end
        end
        if ~good
            return
        end
    end
    x = G.sample;
    while Zn.contains(x)
        x = G.sample;
    end
    if ~(G.elementOrder(x) == 2 && all(cellfun(@(a) G.isIdentity(G.composeAll({x a x a})), Zn.generators)))
        return
    end
    Rcyclic = replab.Atlas.recognizeCyclic(Zn);
    a = Rcyclic.isomorphism.imageElement(Rcyclic.atlasEntry.group.generator(1));
    entry = replab.AtlasEntry.dihedral(n);
    R = replab.AtlasResult(G, entry, entry.group.isomorphismByImages(G, 'images', {x a}));
end
