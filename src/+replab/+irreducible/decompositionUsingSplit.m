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
            trivial{1,end+1} = sub;
        else
            nontrivial{1,end+1} = sub;
        end
    end
    context = replab.Context.make;
    % regroup equivalent representations
    C = rep.commutant.sampleInContext(context, 1);
    C = (C + C')/2;
    nNT = length(nontrivial);
    mask = logical(zeros(nNT, nNT));
    tol = replab.Parameters.doubleEigTol;
    for i = 1:nNT
        subI = nontrivial{i};
        C2 = subI.U * C;
        for j = 1:nNT
            subJ = nontrivial{j};
            mask(i,j) = replab.isNonZeroMatrix(C2 * subJ.U', tol);
        end
    end
    cc = replab.Partition.connectedComponents(mask).blocks;
    % the blocks of the partition cc represent isotypic components
    nNT = length(cc);
    NT = cell(1, nNT);
    for i = 1:nNT
        iso = cc{i};
        subreps = nontrivial(cc{i});
        NT{i} = replab.irreducible.buildIsotypic(rep, samples, subreps);
    end
    % Sort by dimension first and then multiplicity
    ranks = cellfun(@(iso) iso.irrepDimension*100000 + iso.multiplicity, NT);
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
