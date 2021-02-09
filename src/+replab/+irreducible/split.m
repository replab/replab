function [trivial, nonTrivialIrreps] = split(rep, sample1, sample2, forceNonUnitaryAlgorithms)
% Splits a representation into irreducible subrepresentations
%
% First extracts the trivial subrepresentations before splitting the nontrivial invariant subspace.
%
% It uses the first sample to establish a basis, and the second sample to
%
% 1) validate irreducibility, and
% 2) find and fix the form of real representations.
%
% The second sample can be reused later to find equivalent subrepresentations and harmonize their bases.
%
% Args:
%   rep (`+replab.Rep`): Representation to decompose
%   sample1 (double(\*,\*)): First sample of ``rep.commutant``
%   sample2 (double(\*,\*)): Second sample of ``rep.commutant``
%   forceNonUnitaryAlgorithms (logical): Whether to force the use of algorithms for not necessarily unitary representations
%
% Returns
% -------
%   trivial: `+replab.Isotypic`
%     Trivial component
%   nontrivialIrreps: cell(1,\*) of `+replab.SubRep`
%     cell(1,\*) of `+replab.SubRep`: Irreducible subrepresentations of ``rep``
    fnua = forceNonUnitaryAlgorithms;
    partition = rep.invariantBlocks;
    D = rep.dimension;
    Itriv = zeros(D, 0);
    Ptriv = zeros(0, D);
    nontrivialIrreps = cell(1, 0);
    candidates = cell(1, 0);
    for i = 1:partition.nBlocks
        block = partition.block(i);
        d = length(block);
        [trivial, nontrivial] = replab.irreducible.blockTrivialSplit(rep, block, forceNonUnitaryAlgorithms);
        % concatenate the trivial subspaces
        Itriv = horzcat(Itriv, trivial.injection('double/sparse'));
        Ptriv = vertcat(Ptriv, trivial.projection('double/sparse'));
        candidates = horzcat(candidates, replab.irreducible.absoluteSplitInParent(nontrivial, sample1, fnua));
    end
    nonTrivialIrreps = cell(1, 0);
    while ~isempty(candidates)
        splitFurther = cell(1, 0);
        for i = 1:length(candidates)
            c = candidates{i};
            irreps = replab.irreducible.identifyIrrepsInParent(c, sample2);
            if isempty(irreps)
                splitFurther{1,end+1} = c;
            else
                nonTrivialIrreps = horzcat(nonTrivialIrreps, irreps);
            end
        end
        candidates = cell(1, 0);
        if ~isempty(splitFurther)
            % we need fresh samples
            X = rep.commutant.sample;
            for i = 1:length(splitFurther)
                sub = splitFurther{i};
                candidates = horzcat(candidates, replab.irreducible.absoluteSplitInParent(sub, X, fnua));
            end
        end
    end
end
