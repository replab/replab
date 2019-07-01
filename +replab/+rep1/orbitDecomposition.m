function [orbits bases] = orbitDecomposition(rep)
% Identifies subrepresentations from the sparsity pattern of a representation
%
% Used as a first step in the decomposition into irreducibles
    assert(isa(rep.group, 'replab.FinitelyGeneratedGroup'));
    mask = false(rep.dimension, rep.dimension);
    for i = 1:rep.group.nGenerators
        rho = rep.image(rep.group.generator(i));
        mask = mask | (abs(rho) > replab.Settings.doubleEigTol);
    end
    orbits = replab.Partition.connectedComponents(mask);
    n = orbits.nBlocks;
    bases = cell(1, n);
    for b = 1:n
        block = orbits.block(b);
        d = length(block);
        bases{b} = sparse(block, 1:d, ones(1, d), d, rep.dimension);
    end
end
