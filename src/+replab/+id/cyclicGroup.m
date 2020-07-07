function [pres x] = cyclicGroup(G)
% Recognizes if the given group is the cyclic group and provides the group generator
%
% The standard presentation is ``<x| x^n = id>``
    assert(~G.isTrivial, 'The case of the trivial group is not implemented');
    if G.order > 2^53-1
        return
    end
    n = double(G.order);
    pres = [];
    x = [];
    if ~G.isCyclic
        return
    end
    x = G.sample;
    while G.elementOrder(x) ~= G.order
        x = G.sample;
    end
    pres = replab.Presentation(1, {ones(1,n)})
end
