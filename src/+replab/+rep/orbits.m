function O = orbits(rep)
% Identifies subrepresentations from the sparsity pattern of a representation
    if isa(rep.group, 'replab.FiniteGroup')
        mask = false(rep.dimension, rep.dimension);
        for i = 1:rep.group.nGenerators
            rho = rep.image(rep.group.generator(i));
            mask = mask | (abs(rho) > replab.Settings.doubleEigTol);
        end
        O = replab.Partition.connectedComponents(mask);
    else
        O = replab.Partition.fromBlockIndices(ones(1, rep.dimension));
    end
end
