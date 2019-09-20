function I = decomposition(rep)
% Decomposes the given representation into irreducible subrepresentations
%
% Args:
%   rep (replab.Rep): Representation to decompose
%
% Returns:
%   replab.Irreducible: Formal irreducible decomposition
    assert(isa(rep, 'replab.Rep'));
    if rep.overR
        trivialRealDivisionAlgebra = replab.DivisionAlgebra.real;
    else
        trivialRealDivisionAlgebra = [];
    end
    % use recursively replab.rep.decompose
    irreps = {};
    rest = {rep.subRepUnitary(speye(rep.dimension))};
    while ~isempty(rest)
        head = rest{1};
        rest = rest(2:end);
        if head.isExtraTrue('isIrreducible')
            
            irreps{1,end+1} = head;
        else
            subs = replab.rep.decompose(head);
            rest = horzcat(rest, subs);
        end
    end
    trivial = cell(1, 0);
    nontrivial = cell(1, 0);
    for i = 1:length(irreps)
        sub = irreps{i};
        if sub.isExtraTrue('hasTrivialSubspace')
            irrep = replab.Irrep(rep, sub.U0, trivialRealDivisionAlgebra);
            trivial{1,end+1} = irrep;
        elseif sub.isExtraFalse('hasTrivialSubspace')
            nontrivial{1,end+1} = sub;
        else
            error('Irreducible subrepresentations must have had detection of trivial equivalence');
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
