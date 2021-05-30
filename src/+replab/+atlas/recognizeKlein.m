function R = recognizeKlein(G)
% Recognizes if the given group is the Klein four-group and provides the generators according to the standard presentation
%
% The standard presentation is ``<x, a| a^2 = x^2 = id, x a x^-1 = a>``
    R = [];
    if G.order ~= 4
        return
    end
    if G.isCyclic
        return
    end
    x = G.generator(1);
    a = [];
    for i = 2:G.nGenerators
        g = G.generator(i);
        if ~G.eqv(x, g)
            a = g;
            break
        end
    end
    assert(~isempty(a));
    entry = replab.AtlasEntry.kleinFourGroup;
    R = replab.AtlasResult(G, entry, entry.group.isomorphismByImages(G, 'images', {x a}));
end
