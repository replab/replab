function sub = splitAllOnes(rep, trivialSamples, hermitianInvariantSamples)
% Checks if the vector of all ones is invariant
    if ~isa(rep.group, 'replab.FiniteGroup') || ~isequal(rep.isUnitary, true) || rep.dimension < 2
        sub = replab.DispatchNext('Split of all ones works only on unitary finite group representations of dimension > 1');
        return
    end
    d = rep.dimension;
    v = ones(d, 1);
    for i = 1:rep.group.nGenerators
        v1 = rep.matrixRowAction(rep.group.generator(i), v);
        if ~isequal(v, v1)
            sub = replab.DispatchNext('Vector of all ones is not invariant');
            return
        end
    end
    sub = cell(1, 2);
    U1 = ones(d, 1)/sqrt(d);
    U2 = replab.rep.standardBasis(d);
    U2 = U2./sqrt(diag(U2*U2'));
    U2 = U2';
    U2 = U2(:,2:end);
    sub1 = replab.SubRep(rep, U1, U1');
    sub1.trivialDimension = 1;
    sub1.isIrreducible = true;
    sub1.frobeniusSchurIndicator = 1;
    sub2 = replab.SubRep(rep, U2, U2');
    sub = {sub1 sub2};
end
