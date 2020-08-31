function permCT = SymmetricGroupCharacterTable(n)
    permGroup = replab.S(n);
    allData = replab.sym.CTData.instance(n);
    nData = allData(n);
    irreps = cellfun(@(part) strrep(replab.shortStr(replab.sym.findPartitions.conjugatePart(part)), ' ', ''),nData.partitionList,'UniformOutput', false);
    [classes partCell] = replab.sym.findConjClasses(permGroup);
    classNames = cellfun(@(c) strrep(replab.shortStr(c), ' ', ''), partCell, 'uniform', 0);
    chars = replab.cyclotomic.fromDoubles(permCTValues(n));
    permCT = replab.CharacterTable(permGroup, replab.ConjugacyClasses(permGroup, classes), chars, 'irrepNames', irreps, 'classNames', classNames);

    function CT = permCTValues(n)
        nParts = nData.nParts;
        CTunNormed = zeros(nParts);
        for row = 1:nParts
            partition = nData.partitionList{row};
            lin = replab.sym.SymPoly.powerToMonomial(partition,n);
            CTunNormed(:,row) = lin.coeffs';
        end
        CT = grahamSchmidt(CTunNormed,@(x,y) round(sum(x.*y.*nData.conjSizes/nData.fact)));
        CT = flip(round(CT));
        function mat = grahamSchmidt(mat,innerProd)
            d = nData.nParts;
            for i = d:-1:2
                for j = 1:(i-1)
                    mat(j,:) = mat(j,:)-innerProd(mat(j,:),mat(i,:))*mat(i,:);
                end
            end
        end
    end
end
