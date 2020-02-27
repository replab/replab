function hi = harmonizeIsotypic(iso, context)
% Harmonizes an isotypic component
    n = iso.nIrreps;
    C = iso.parent.commutant.sampleInContext(context, 1);
    if iso.overR
        irr1 = replab.irreducible.canonicalDivisionAlgebra(iso.irrep(1));
        W = irr1.A_internal;
        Winv = irr1.Ainv_internal
    end
    irreps1 = cell(1, n);
    irreps1{1} = replab.rep.collapse(irr1);
    for i = 2:n
        A = iso.changeOfBasis(1, i, context);
        Ainv = iso.changeOfBasis(i, 1, context);

        irri = iso.irrep(i).
    end
    sub1 = sub{1};
    if rep.overR
        first1 = replab.irreducible.canonicalDivisionAlgebra(first, context);

        % identifies the type of real representations and
        % make the division algebra basis canonical
        realType = replab.irreducible.computeRealType(rep, samples, first);
        switch realType
          case 'R'
            ii = replab.irreducible.Info(realType, []);
            first = rep.subRepUnitary(first.U, first.niceBasis, ii);
          case 'C'
            W = replab.irreducible.enforceComplexEncoding(rep, samples, first);
            ii = replab.irreducible.Info(realType, true);
            first = rep.subRepUnitary(W*first.U, [], ii);
          case 'H'
            W = replab.irreducible.enforceQuaternionEncoding(rep, samples, first);
            ii = replab.irreducible.Info(realType, true);
            first = rep.subRepUnitary(W*first.U, [], ii);
        end
    end
    copies = cell(1, n);
    copies{1} = first;
    for j = 2:n
        other = sub{j};
        W = replab.irreducible.findCommonBasis(rep, samples, first, other);
        copies{j} = rep.subRepUnitary(W * other.U, [], first.irrepInfo);
    end
    iso = replab.Isotypic(rep, copies);

end
