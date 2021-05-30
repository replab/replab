function R = recognizeCyclic(G)
% Recognizes if the given group is the cyclic group and provides the group generator
%
% The standard presentation is ``<x| x^n = id>``
    R = [];
    if G.isTrivial
        return
    end
    if G.order > 2^53-1
        return
    end
    n = double(G.order);
    x = [];
    if ~G.isCyclic
        return
    end
    x = G.sample;
    while G.elementOrder(x) ~= G.order
        x = G.sample;
    end
    entry = replab.AtlasEntry.cyclic(n);
    R = replab.AtlasResult(G, entry, entry.group.isomorphismByImages(G, 'images', {x}));
end
