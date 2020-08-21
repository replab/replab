function sub = splitBlocks(rep, trivialSamples, hermitianInvariantSamples)
% Identifies an explicit block form if present, and split accordingly
    if isa(rep, 'replab.rep.DirectSumRep') && rep.nFactors > 1
        start = 1;
        d = rep.dimension;
        n = rep.nFactors;
        sub = cell(1, n);
        for i = 1:rep.nFactors
            di = rep.factor(i).dimension;
            B_internal = sparse(start:start+di-1, 1:di, ones(1,di), d, di);
            E_internal = B_internal';
            sub{1, i} = rep.subRep(B_internal, E_internal);
            start = start + di;
        end
    elseif isa(rep.group, 'replab.rep.FiniteGroup') && ~rep.group.isTrivial
        mask = (rep.image_internal(rep.group.generator(1)) == 0);
        for i = 2:rep.group.nGenerators
            mask = mask | (rep.image_internal(rep.group.generator(i)) == 0);
        end
        partition = replab.UndirectedGraph.fromAdjacencyMatrix(mask).connectedComponents;
        if partition.nBlocks == 1
            sub = replab.DispatchNext('No blocks identified');
            return
        end
        n = partition.nBlocks;
        sub = cell(1, n);
        d = rep.dimension;
        for i = 1:n
            b = partition.block(i);
            di = length(b);
            B_internal = sparse(b, 1:di, ones(1, di), d, di);
            E_internal = B_internal';
            sub{1, i} = rep.subRep(B_internal, E_internal);
        end
    else
        sub = replab.DispatchNext('Not a covered case');
    end
end
