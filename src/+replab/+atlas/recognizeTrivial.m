function R = recognizeTrivial(G)
% Recognizes if the given group is the trivial group and provides the generators according to the standard presentation
%
% The standard presentation is ``<x| x = id>``
    R = [];
    if ~G.isTrivial
        return
    end
    entry = replab.AtlasEntry.trivial;
    R = replab.AtlasResult(G, entry, entry.group.isomorphismByImages(G));
end
