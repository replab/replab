function irreps = absoluteSplit_usingInvariantBlocks(rep, forceNonUnitaryAlgorithms)
% Splits a representation using existing block structure
%
% See `+replab.Rep.absoluteSplit`
%
% First extracts the trivial subrepresentations before splitting the nontrivial invariant subspace;
% this is done block by block.
%
% Args:
%   rep (`+replab.Rep`): Representation to decompose
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: Subrepresentations of ``rep``
    partition = rep.invariantBlocks;
    D = rep.dimension;
    Itriv = sparse(D, 0);
    Ptriv = sparse(0, D);
    nontrivialIrreps = cell(1, 0);
    for i = 1:partition.nBlocks
        block = partition.block(i);
        [trivial, nontrivial] = replab.rep.blockSplitTrivial(rep, block, forceNonUnitaryAlgorithms);
        Itriv = horzcat(Itriv, trivial.injection('double/sparse'));
        Ptriv = vertcat(Ptriv, trivial.projection('double/sparse'));
        nontrivialIrreps = horzcat(nontrivialIrreps, nontrivial.absoluteSplitInParent);
    end
    subT = rep.subRep(Itriv, 'projection', Ptriv);
    trivialIrreps = replab.Isotypic.fromTrivialSubRep(subT).irreps;
    irreps = horzcat(trivialIrreps, nontrivialIrreps);
end
