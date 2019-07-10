function O = orbits(rep)
% Identifies subrepresentations from the sparsity pattern of a representation
    assert(isa(rep.group, 'replab.FinitelyGeneratedGroup'));
    mask = false(rep.dimension, rep.dimension);
    for i = 1:rep.group.nGenerators
        rho = rep.image(rep.group.generator(i));
        mask = mask | (abs(rho) > replab.Settings.doubleEigTol);
    end
    O = replab.Partition.connectedComponents(mask);
end
