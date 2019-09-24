function I = decompositionUsingSplit(rep)
% Generic method for representation decomposition
%
% First it splits the representation into irreducibles, before recognizing which
% irreducible representations are part of the same isotypic component.
    replab.irreducible.tell('decompositionUsingSplit')
    samples = replab.irreducible.OnDemandSamples(rep);
    % obtain all irreps
    replab.irreducible.tell('down')
    mainSub = rep.subRepUnitaryByIntegerBasis(speye(rep.dimension));
    subs = replab.irreducible.split(rep, samples, mainSub);
    replab.irreducible.tell('up')    
    % sort by trivial / non trivial
    trivial = {};
    nontrivial = {};
    for i = 1:length(subs)
        sub = subs{i};
        if isequal(sub.irrepInfo.label, '1') % is it trivial
            trivial{1,end+1} = sub;
        else
            nontrivial{1,end+1} = sub;
        end
    end
    % regroup equivalent representations
    C = samples.commutantSample(2);
    C = (C + C')/2;
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
    % the blocks of the partition cc represent isotypic components
    nNT = length(cc);
    NT = cell(1, nNT);
    for i = 1:nNT
        iso = cc{i};
        subreps = nontrivial(cc{i});
        NT{i} = replab.irreducible.buildIsotypic(rep, samples, subreps);
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
