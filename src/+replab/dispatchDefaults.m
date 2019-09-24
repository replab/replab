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
    
    replab.dispatch('register', 'replab.rep.split', 'ReduceBlocks', 500, @replab.rep.splitPermutations);
    replab.dispatch('register', 'replab.rep.split', 'UsingCommutant', 0, @replab.rep.splitUsingCommutant);
end
