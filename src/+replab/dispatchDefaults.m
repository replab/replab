function dispatchDefaults
% Register here the default functions used in dispatch for the core of RepLAB
    
% Equivariant construction
    replab.dispatch('register', 'replab.makeEquivariant', 'ForRepByImages', 10, ...
                    @(repR, repC) replab.equivariant.ForRepByImages(repR, repC));
    replab.dispatch('register', 'replab.makeEquivariant', 'ForFiniteGroup', 5, ...
                    @(repR, repC) replab.equivariant.ForFiniteGroup(repR, repC));
    % Default method, works for all compact groups
    replab.dispatch('register', 'replab.makeEquivariant', 'ForCompactGroup', 0, ...
                    @(repR, repC) replab.equivariant.ForCompactGroup(repR, repC));

    replab.dispatch('register', 'replab.irreducible.decomposition', 'UsingSplit', 0, ...
                    @(rep) replab.irreducible.decompositionUsingSplit(rep));
                
    replab.dispatch('register', 'replab.irreducible.split', 'ReduceBlocks', 500, ...
                    @(rep, samples, sub) replab.irreducible.splitPermutations(rep, samples, sub));
    replab.dispatch('register', 'replab.irreducible.split', 'UsingCommutant', 0, ...
                    @(rep, samples, sub) replab.irreducible.splitUsingCommutant(rep, samples, sub));
end
