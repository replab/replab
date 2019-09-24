function sub = splitPermutations(rep, samples, U)
% Splits a permutation representation
    disp('sp1')
    if nargin > 2
        error('replab:dispatch:tryNext', 'try next');
    end
    disp('sp2')
    if ~isa(rep.group, 'replab.FiniteGroup')
        error('replab:dispatch:tryNext', 'try next');
    end
    disp('sp3')
    nG = rep.group.nGenerators;
    d = rep.dimension;
    % Transforms all group generators into permutations, bail out if impossible
    G = zeros(nG, d);
    for i = 1:nG
        rhog = rep.image(rep.group.generator(i));
        try
            G(i,:) = replab.Permutations.fromMatrix(rhog);
        catch ME
            error('replab:dispatch:tryNext', 'try next');
        end
    end
    disp('sp4')
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
        % construct the orthogonal complement basis in that block
        Unontrivial = sparse(dB-1, d);
        Unontrivial(:, block) = null(ones(1, dB))'; % TODO: algebraic nice balanced exact expression?
        otherSub = replab.rep.split(rep, samples, Unontrivial);
        sub = horzcat(sub, otherSub);
    end
end
