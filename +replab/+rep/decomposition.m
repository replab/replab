function I = decomposition(rep)
% Decomposes the given representation into (structured) irreducible subrepresentations
%
% Args:
%   rep (:class:`+replab.Rep`): Representation to decompose
%
% Returns:
%   :class:`+replab.Irreducible`: A structured decomposition into irreducible representations
    assert(isa(rep, 'replab.Rep'));
    assert(rep.overR || rep.overC);
    if isa(rep.group, 'replab.PermutationGroup')
        n = rep.group.domainSize;
        if rep.group.order == factorial(vpi(n))
            I = replab.rep.decomposeSymmetric(rep);
            if ~isequal(I, [])
                return
            end
        end
    end
    if rep.overR
        trivialRealDivisionAlgebra = replab.DivisionAlgebra.real;
    else
        trivialRealDivisionAlgebra = [];
    end
    trivial = cell(1, 0);
    nontrivial = cell(1, 0);
    sub1 = replab.rep.orbitDecomposition(rep);
    for i = 1:length(sub1)
        s1 = sub1{i};
        [Ut Ur] = replab.rep.extractTrivial(s1);
        for j = 1:size(Ut, 1)
            Uthis = Ut(j, :) * s1.U;
            trivial{1, end+1} = replab.Irrep(rep, Uthis, trivialRealDivisionAlgebra);
        end
        if size(Ur, 2) > 0
            s2 = s1.subRep(Ur);
            sub3 = replab.rep.decomposeUsingCommutant(s2);
            for j = 1:length(sub3)
                s3 = sub3{j};
                nontrivial{1, end+1} = s3.collapseParent.collapseParent;
            end
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
            mask(i,j) = replab.isNonZeroMatrix(subI.U * C * subJ.U', tol);
        end
    end
    cc = replab.Partition.connectedComponents(mask).blocks;
    nNT = length(cc);
    NT = cell(1, nNT);
    for i = 1:nNT
        iso = cc{i};
        first = nontrivial{iso(1)};
        if rep.overR
            realType = replab.rep.realType(first);
            switch realType
              case 'R'
                firstIrrep = replab.Irrep(rep, first.U0, ...
                                          replab.DivisionAlgebra.real);
              case 'C'
                W = replab.rep.enforceComplexEncoding(first);
                firstIrrep = replab.Irrep(rep, W * first.U, ...
                                          replab.DivisionAlgebra.complex);
              case 'H'
                W = replab.rep.enforceQuaternionEncoding(first);
                firstIrrep = replab.Irrep(rep, W * first.U, ...
                                          replab.DivisionAlgebra.quaternion);
            end
        else
            firstIrrep = replab.Irrep(rep, first.U0, []);
        end
        copies = cell(1, length(iso));
        copies{1} = firstIrrep;
        for j = 2:length(iso)
            other = nontrivial{iso(j)};
            E = replab.rep.findCommonBasis(rep, firstIrrep, other, C);
            copies{j} = replab.Irrep(rep, E * other.U, ...
                                     firstIrrep.realDivisionAlgebra);
        end
        NT{i} = replab.Isotypic(rep, copies);
    end
    % Sort by dimension first and then multiplicity
    ranks = cellfun(@(iso) iso.copyDimension*100000 + iso.multiplicity, NT);
    [~, I] = sort(ranks);
    NT = NT(I);
    if length(trivial) > 0
        trivialIsotypic = replab.Isotypic(rep, trivial);
        components = horzcat({trivialIsotypic}, NT);
    else
        components = NT;
    end
    I = replab.Irreducible(rep, components);
end
