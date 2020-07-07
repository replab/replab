function [name, pres, f, g] = identifyGroup(G)
% Attempts to identify the given group
%
% For now, we use the notation conventions on the wiki https://groupprops.subwiki.org/wiki/Main_Page
%
% Example:
%    >>> replab.identifyGroup(replab.DihedralGroup(3))
%      'Dihedral group of order 6 and degree 3'
% Args:
%   G (`+replab.NiceFiniteGroup`): The group to recognize
%
% Returns
% -------
%   name:
%     charstring: A human-readable description of the group
%   pres:
%     `+replab.Presentation`: The standard presentation of the recognized group
%   f:
%     `+replab.Morphism`: An isomorphism from ``G`` to a group that obeys the presentation
%   g:
%     `+replab.Morphism`: Inverse of the isomorphism ``f``

    [pres x] = replab.id.cyclicGroup(G);
    if ~isempty(pres)
        source = replab.PermutationGroup.of(x);
        n = double(source.order);
        name = sprintf('Cyclic group of order %d', n);
        target = replab.CyclicGroup(n);
        f = source.morphismByImages(target, target.generators);
        g = target.morphismByImages(source, source.generators);
        return
    end

    [pres x a] = replab.id.dihedralGroup(G);
    if ~isempty(pres)
        source = replab.PermutationGroup.of(x, a);
        d = double(source.order)/2;
        name = sprintf('Dihedral group of order %d and degree %d', 2*d, d);
        target = replab.DihedralGroup(d);
        f = source.morphismByImages(target, target.generators);
        g = target.morphismByImages(source, source.generators);
        return
    end

end
