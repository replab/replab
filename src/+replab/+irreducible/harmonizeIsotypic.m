function hi = harmonizeIsotypic(iso)
    %%% WORK IN PROGRESS

    first = sub{1};
    if rep.overR
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
