function sub = splitPermutations(rep, samples, sub)
% Splits a permutation representation
    replab.irreducible.tell('Attempting splitPermutations');
    if ~replab.iseye(sub.H_internal)
        sub = replab.DispatchNext('Not full representation');
        return
    end
    if ~isa(rep.group, 'replab.FiniteGroup')
        sub = replab.DispatchNext('Not finite group');
        return
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
            sub = replab.DispatchNext('Images are not all permutations');
            return
        end
    end
    replab.irreducible.tell('Running splitPermutations');
    trivialIrrepInfo = replab.irreducible.TrivialInfo(rep.field);
    P = replab.Partition.permutationsOrbits(G);
    sub = {};
    % For each block, extract the all-ones trivial representation, decompose the orthogonal component
    % in that block separately
    for i = 1:P.nBlocks
        block = P.block(i);
        dB = length(block);
        % construct the trivial representation
        Vtrivial = sparse(ones(1, dB), block, ones(1, dB), 1, d);
        subTrivial = rep.subRepUnitaryByIntegerBasis(Vtrivial, trivialIrrepInfo);
        sub{1, end+1} = subTrivial;
        if dB > 1
            stdBasis = replab.rep.standardBasis(dB);
            % construct the orthogonal complement basis in that block
            Vnontrivial = sparse(dB-1, d);
            Vnontrivial(:, block) = stdBasis(2:end, :); % TODO: algebraic nice balanced exact expression?
            rest = rep.subRepUnitaryByIntegerBasis(Vnontrivial);
            replab.irreducible.tell('down');
            otherSub = replab.irreducible.split(rep, samples, rest);
            replab.irreducible.tell('up');
            sub = horzcat(sub, otherSub);
        end
    end
end
