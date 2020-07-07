function [name, pres, f, g] = identifyGroup(G)

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
