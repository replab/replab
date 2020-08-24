function [basis, basisMat] = findAllRepCoeffs2(rep)
    n = rep.group.domainSize;
    [basisCell1,irreps] =  replab.sym.repBlockBasis(rep);
    adjSwapIms = arrayfun(@(j) rep.image(swap(j)),1:(n-1),'UniformOutput',false);
    basis = arrayfun(@(i) irrepBasis(i),1:numel(irreps),'UniformOutput',false);
    basisMat = cellfun(@(matCell) [matCell{:}],basis,'UniformOutput',false);
    basisMat = [basisMat{:}];
    
    function irrepBasCell = irrepBasis(i)
        dim = replab.sym.findPartitions.dimension(irreps{i});
        if dim == 1
            irrepBasCell = basisCell1(i);
            return
        end
        [connectList,transList,offDiagList,diagList] = replab.sym.tableauxTree(irreps{i});
        irrepBasCell = cell(1,dim);
        irrepBasCell(1) = basisCell1(i);
        for index = 2:(dim)
            rec(index);
        end
       
        function shiftedBasis = rec(index)
            if ~isempty(irrepBasCell{index})
                shiftedBasis = irrepBasCell{index};
                return
            end
                parentIndex = connectList(index-1);
                basisOfParent = rec(parentIndex);
                shiftedBasis = (adjSwapIms{transList(index-1)}*basisOfParent-...
                    diagList(index-1)*basisOfParent)/offDiagList(index-1);
                irrepBasCell{index} = shiftedBasis;
            end
    end

    function perm = swap(k)
        perm = 1:n;
        perm([k k+1]) = perm([k+1 k]);
    end

end