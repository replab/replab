function I1 = conjugateIrreducible(I, r, W)
% Left conjugates the r-th irreducible representation by W
    m = I.multiplicities;
    d = I.dimensions;
    if ~isequal(I.divisionAlgebras{r}.shortName, 'R')
        warning('Unknown behavior when applied to representations not of real type.');
    end
    sizes = m.*d;
    % Size of previous components to skip
    shift = sum(sizes(1:r-1));
    U1 = I.U;
    for i = 1:I.multiplicities
        ind = shift+(1:d(r));
        U1(:, ind) = U1(:, ind) * W'; % TODO verify
        shift = shift + d(r);
    end
    I1 = replab.IrreducibleDecomposition(I.rep, U1, I.dimensions, I.multiplicities, I.divisionAlgebras);
end
