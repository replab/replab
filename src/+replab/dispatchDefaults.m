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
    
    % If it is already irreducible
    replab.dispatch('register', 'replab.rep.decompose', 'Irreducible', 1000, @replab.rep.decomposeIrreducible);
    replab.dispatch('register', 'replab.rep.decompose', 'ReduceBlocks', 500, @replab.rep.decomposeReduceBlocks);
    replab.dispatch('register', 'replab.rep.decompose', 'ExtractAllOnes', 200, @replab.rep.decomposeExtractAllOnes);
    replab.dispatch('register', 'replab.rep.decompose', 'ExtractTrivial', 100, @replab.rep.decomposeExtractTrivial);
    replab.dispatch('register', 'replab.rep.decompose', 'UsingCommutant', 0, @replab.rep.decomposeUsingCommutant);
end
