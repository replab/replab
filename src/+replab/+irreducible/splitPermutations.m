function sub = splitPermutations(rep, samples, U)
% Splits a permutation representation
    replab.irreducible.tell('Attempting splitPermutations');
    if nargin > 2 && ~replab.iseye(U)
        replab.irreducible.tell('Not full rep');
        error('replab:dispatch:tryNext', 'try next');
    end
    if ~isa(rep.group, 'replab.FiniteGroup')
        replab.irreducible.tell('Not finite group');
        error('replab:dispatch:tryNext', 'try next');
    end
    nG = rep.group.nGenerators;
    d = rep.dimension;
    % Transforms all group generators into permutations, bail out if impossible
    G = zeros(nG, d);
    for i = 1:nG
        rhog = rep.image(rep.group.generator(i));
        try
            G(i,:) = replab.Permutations.fromMatrix(rhog);
        catch ME
            replab.irreducible.tell('Not permutation image');
            error('replab:dispatch:tryNext', 'try next');
        end
    end
    replab.irreducible.tell('Running splitPermutations');
    if rep.overR
        trivialIrrepInfo = replab.IrrepInfo('1', 'R', []);
    else
        trivialIrrepInfo = replab.IrrepInfo('1', [], []);
    end
    P = replab.Partition.permutationsOrbits(G);
    sub = {};
    % For each block, extract the all-ones trivial representation, decompose the orthogonal component
    % in that block separately
    for i = 1:P.nBlocks
        block = P.block(i);
        dB = length(block);
        % construct the trivial representation
        Vtrivial = sparse(ones(1, dB), block, ones(1, dB), 1, d);
        NBtrivial = replab.NiceBasis.fromIntegerBasis(Vtrivial);
        subTrivial = rep.subRepUnitary(NBtrivial.U, NBtrivial, trivialIrrepInfo);
        sub{1, end+1} = subTrivial;
        if dB > 1
            % construct the orthogonal complement basis in that block
            Unontrivial = sparse(dB-1, d);
            Unontrivial(:, block) = null(ones(1, dB))'; % TODO: algebraic nice balanced exact expression?
            replab.irreducible.tell('down');
            otherSub = replab.irreducible.split(rep, samples, Unontrivial);
            replab.irreducible.tell('up');
            sub = horzcat(sub, otherSub);
        end
    end
end
