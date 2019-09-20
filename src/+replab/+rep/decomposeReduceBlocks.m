function sub = decomposeReduceBlocks(rep)
% If a representation has a block diagonal structure, decompose orbits separately
    if rep.isExtraTrue('reducedBlocks')
        error('replab:dispatch:tryNext', 'try next');
    end 
    % This is what we known about the children
    extra = struct('reducedBlocks', true);
    if rep.isExtraFalse('hasTrivialSubspace')
        extra.hasTrivialSubspace = false;
    end
    sub = {}; % start with empty subrepresentation array
    % Compute the block structure of the representation 
    O = replab.rep.orbits(rep);
    n = O.nBlocks;
    for b = 1:n
        block = O.block(b);
        d = length(block);
        % basis vectors are row vectors
        basis = sparse(1:d, block, ones(1, d), d, rep.dimension);
        % it's ok to pass a sparse matrix to subRepUnitary as they are hidden from the user
        sub{1, end+1} = rep.subRepUnitary(basis, extra).collapseParent;
    end
end
