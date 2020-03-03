function I = decompositionUsingSplit(rep)
% Generic method for representation decomposition
%
% First it splits the representation into irreducibles, before recognizing which
% irreducible representations are part of the same isotypic component.
    context = replab.Context.make;
    subs = rep.splitIntoIrreducibles(context);
    context.close;
    % sort by trivial / non trivial
    trivial = {};
    nontrivial = {};
    for i = 1:length(subs)
        s = subs{i};
        if s.trivialDimension == 1
            trivial{1,end+1} = s;
        else
            nontrivial{1,end+1} = s;
        end
    end
    context = replab.Context.make;
    % regroup equivalent representations
    C = rep.commutant.sampleInContext(context, 1);
    nNT = length(nontrivial);
    mask = logical(zeros(nNT, nNT));
    tol = replab.Parameters.doubleEigTol;
    for i = 1:nNT
        subI = nontrivial{i};
        CI = subI.E_internal * C;
        for j = 1:nNT
            subJ = nontrivial{j};
            mask(i,j) = replab.isNonZeroMatrix(CI * subJ.B_internal, tol);
        end
    end
    cc = replab.Partition.connectedComponents(mask).blocks;
    % the blocks of the partition cc represent isotypic components
    nNT = length(cc);
    NT = cell(1, nNT);
    for i = 1:nNT
        iso = cc{i};
        subreps = nontrivial(cc{i});
        NT{i} = replab.irreducible.harmonizeIsotypic(replab.irreducible.buildIsotypic(rep, subreps), context);
    end
    % Sort by dimension first and then multiplicity
    dims = cellfun(@(iso) iso.irrepDimension, NT);
    muls = cellfun(@(iso) iso.multiplicity, NT);
    [~, I] = sortrows([dims(:) muls(:)]);
    NT = NT(I);
    if length(trivial) > 0
        trivialIsotypic = replab.irreducible.harmonizeIsotypic(replab.irreducible.buildIsotypic(rep, trivial), context);
        components = horzcat({trivialIsotypic}, NT);
    else
        components = NT;
    end
    I = replab.Irreducible(rep, components);
end
