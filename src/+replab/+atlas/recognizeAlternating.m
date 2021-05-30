function R = recognizeAlternating(G)
% Recognizes the alternating group and returns the generators corresponding to the standard presentation
    R = [];
    [n r] = replab.util.unfactorial(G.order*2);
    if r ~= 0
        return
    end
    n = double(n);
    C = G.conjugacyClasses.classes;
    entry = replab.AtlasEntry.alternating(n);
    for i = 1:length(C)
        S = C{i};
        s = S.representative;
        if G.elementOrder(s) == n-2
            for j = 1:length(C)
                T = C{j};
                if G.elementOrder(T.representative) == 3
                    U = T.elements;
                    for k = 1:length(U)
                        t = U{k};
                        if entry.group.isMorphismByImages(G, 'images', {s t})
                            if G.subgroup({s, t}).order == G.order
                                R = replab.AtlasResult(G, entry, entry.group.isomorphismByImages(G, 'images', {s t}));
                                return
                            end
                        end
                    end
                end
            end
        end
    end
end
