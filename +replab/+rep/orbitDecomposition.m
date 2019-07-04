function [sub O] = orbitDecomposition(rep)
% Identifies subrepresentations from the sparsity pattern of a representation
%
% Those subrepresentations are not necessarily irreducible in general.
%
% Used as a first step in the decomposition into irreducibles
    O = replab.rep.orbits(rep);
    n = O.nBlocks;
    sub = cell(1, n);
    for b = 1:n
        block = O.block(b);
        d = length(block);
        basis = sparse(block, 1:d, ones(1, d), rep.dimension, d);
        sub{b} = rep.subRep(basis);
    end
end
