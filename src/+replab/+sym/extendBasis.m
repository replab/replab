function irrepBasCell = extendBasis(basisElem1,partition,adjSwapIms,isSymb)
% Given one multiplicity space in a symmetric group irrep, generate all
% others, using Eq 4-36 in Chen's Group Representation Theory for Physicists
 % Args:
        % basisElem1 (double(*\,*\): First multiplicity space
        %
        % partition (replab.sym.IntegerPartition): Representation to be decomposed. 
        %
        % adjSwapIms (cell(1,*\) of double(*\,n-1)): Images of adjacent
        % transpositions of the rep being decomposed
        %
        % symb (boolean): Do we want a symbolic result? Here this is
        % equivalent to getting the rational result.
        %
    % Returns:
        % irrepBasCell (cell(1,*\) of double(*\,*\)): All multiplicity spaces,
        % for this irrep
        %
        % irreps: (cell(1,*\) of replab.sym.IntegerPartition):
        % Partitions corresponding to the irreducible components of
        % representation.

        dim = replab.sym.IntegerPartition.dimension(partition);
        if dim == 1
            irrepBasCell = {basisElem1};
            return
        end
        [connectList,transList,offDiagList,diagList] = replab.sym.tableauxTree(partition);
        irrepBasCell = cell(1,dim);
        irrepBasCell(1) = basisElem1;
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
                if ~isSymb
                    shiftedBasis = (adjSwapIms{transList(index-1)}*basisOfParent-...
                        diagList(index-1)*basisOfParent)/offDiagList(index-1);
                else
                    shiftedBasis = (adjSwapIms{transList(index-1)}*basisOfParent-...
                        diagList(index-1)*basisOfParent);
                end   
                irrepBasCell{index} = shiftedBasis;
            end
    end