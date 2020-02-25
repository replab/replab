function iso = buildIsotypic(rep, sub)
% Builds an isotypic canonical component from equivalent subrepresentations
%
% As part of the operation, we identify the division algebra type (if the representations are over R),
% put those division algebras in their canonical basis, and make the irreducible
% representations not only equivalent but identical.
%
% Args:
%   rep (`+replab.Rep`): Representation being decomposed, basis not overlapping
%   sub (cell(1,\*) of `+replab.SubRep`): Equivalent irreducible subrepresentations of ``rep``
%
% Returns:
%   `+replab.Isotypic`: The isotypic component
    assert(length(sub) >= 1, 'Isotypic component cannot be empty');
    m = length(sub);
    d = sub{1}.dimension;
    if isequal(rep.isUnitary, true)
        b = cellfun(@(s) isequal(s.B_internal, s.E_internal'), sub);
        if all(b)
            % all bases are unitary, parent is unitary, we use orthogonality
            iso = replab.Isotypic(rep, sub);
            return
        end
    end
    Bs = cell(1, m);
    for i = 1:m
        Bs{1,i} = sub{i}.B_internal;
    end
    Biso = [Bs{:}];
    subiso = replab.subRep(Biso);
    Fiso = subiso.F_internal;
    sub1 = cell(1, m);
    for i = 1:m
        range = (i-1)*d+(1:d);
        sub1{i} = rep.subRep(Biso(:,range), Fiso(range,:));
        sub1{i}.isIrreducible = sub{i}.isIrreducible;
        sub1{i}.trivialDimension = sub{i}.trivialDimension;
        sub1{i}.frobeniusSchurIndicator = sub{i}.frobeniusSchurIndicator;
        % this translates as this canonicity is a property of images, not of the
        % embedding map
        sub1{i}.isDivisionAlgebraCanonical = sub{i}.isDivisionAlgebraCanonical;
    end
    iso = replab.Isotypic(rep, sub1);
end
