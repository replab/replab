function irreps = split(rep)
% Splits a representation into irreducible subrepresentations, plain version
%
% First extracts the trivial subrepresentations before splitting the nontrivial invariant subspace.
%
% Args:
%   rep (`+replab.Rep`): Representation to decompose
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: Irreducible subrepresentations of ``rep``
    trivial = rep.trivialComponent('double');
    if trivial.dimension == 0
        irreps = cell(1, 0);
        nontrivial = replab.SubRep.identical(rep);
    else
        nontrivial = rep.maschke(trivial);
        irreps = trivial.irreps;
    end
    nontrivialIrreps = nontrivial.splitInParent;
    for i = 1:length(nontrivialIrreps)
        nontrivialIrreps{i}.cache('trivialDimension', 0, '==');
    end
    irreps = horzcat(irreps, nontrivialIrreps);
end
