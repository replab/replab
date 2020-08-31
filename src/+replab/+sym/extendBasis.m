function irrepBasCell = extendBasis(basisCell1,partition,adjSwapIms,isRat)
        dim = replab.sym.IntegerPartition.dimension(partition);
        if dim == 1
            irrepBasCell = basisCell1;
            return
        end
        [connectList,transList,offDiagList,diagList] = replab.sym.tableauxTree(partition);
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
                if ~isRat
                    shiftedBasis = (adjSwapIms{transList(index-1)}*basisOfParent-...
                        diagList(index-1)*basisOfParent)/offDiagList(index-1);
                else
                    shiftedBasis = (adjSwapIms{transList(index-1)}*basisOfParent-...
                        diagList(index-1)*basisOfParent);
                end   
                irrepBasCell{index} = shiftedBasis;
            end
    end