function irreps = decomposeIrreducible(rep)
% If a representation is already irreducible, no need to decompose it further
    if ~rep.isExtraTrue('isIrreducible')
        error('replab:dispatch:tryNext', 'try next');
    end
    irreps = {rep};
end
