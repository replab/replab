function I = decompose(rep)
% Decomposes the given representation into irreducible subrepresentations
    assert(isa(rep, 'replab.Rep'));
    assert(rep.overR || rep.overC);
    sub1 = replab.rep1.orbitDecomposition(rep);
    trivial = cell(1, 0);
    nontrivial = cell(1, 0);
    for i = 1:length(sub1)
        s1 = sub1{i};
        [Ut Ur] = replab.rep1.extractTrivial(s1);
        for j = 1:size(Ut, 2)
            trivial{1, end+1} = s1.subRep(Ut(:, j)).in(rep);
        end
        s2 = s1.subRep(Ur);
        sub3 = replab.rep1.decomposeUsingCommutant(s2);
        for j = 1:length(sub3)
            s3 = sub3{j};
            nontrivial{1, end+1} = s3.in(rep);
        end
    end
    % regroup equivalent representations
    C = rep.commutant.sampleSelfAdjoint;
    nNT = length(nontrivial);
    mask = logical(zeros(nNT, nNT));
    tol = replab.Settings.doubleEigTol;
    for i = 1:nNT
        for j = 1:nNT
            subI = nontrivial{i};
            subJ = nontrivial{j};
            mask(i,j) = replab.isNonZeroMatrix(subI.U' * C * subJ.U, tol);
        end
    end
    cc = replab.Partition.connectedComponents(mask).blocks;
    nNT = length(cc);
    NT = cell(1, nNT);
    for i = 1:nNT
        iso = cc{i};
        first = nontrivial{iso(1)};
        if rep.overR
            realType = replab.rep1.realType(first);
            switch realType
              case 'C'
                W = replab.rep1.enforceComplexEncoding(first);
                first = rep.subRep(first.U * W);
              case 'H'
                W = replab.rep1.enforceQuaternionEncoding(first);
                first = rep.subRep(first.U * W);
            end
        else
            realType = [];
        end
        copies = cell(1, length(iso));
        copies{1} = first;
        for j = 2:length(iso)
            other = nontrivial{iso(j)};
            E = first.U' * C * other.U;
            E = E * sqrt(first.dimension/real(trace(E*E')));
            copies{j} = rep.subRep(other.U * E);
        end
        NT{i} = replab.Isotypic(rep, copies, realType);
    end
    if length(trivial) > 0
        if rep.overR
            trivialRealType = 'R';
        else
            trivialRealType = [];
        end
        trivialIsotypic = replab.Isotypic(rep, trivial, trivialRealType);
        components = horzcat({trivialIsotypic}, NT);
    else
        components = NT;
    end
    I = replab.Irreducible(rep, components);
end
