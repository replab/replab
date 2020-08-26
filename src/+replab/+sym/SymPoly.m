classdef SymPoly
    properties
        coeffs
        m% Sum of partitions/ degree of polynomials
    end
    methods(Static)
        function lin = fromCell(partStandForm)
        %m=degree
            len = numel(partStandForm);
            deg = sum(partStandForm{1});
            allData = replab.sym.CTData.instance;
            mData = allData(deg);
            arr = zeros(1,mData.nParts);
            for j =1:len
                ind = mData.partitions.index(partStandForm{j});
                arr(ind) = arr(ind)+1;
            end
            lin = replab.sym.SymPoly(arr,deg);

        end

    function lin = multRowWithAugmented(rowPart,part)
                len = numel(part);
                if len == 0
                    lin =  replab.sym.SymPoly.fromCell(rowPart);
                    return
                elseif len == 1
                    [rowPart,part] = deal(part,rowPart);
                end
                if (numel(rowPart)~=1)
                    lin =   replab.sym.SymPoly.fromCell({[]});
                    return
                end
                compositions = arrayfun(@(i) addElem(part,i,rowPart),1:len,'UniformOutput',false);
                compositions{len+1} = sort([rowPart part],'descend');
                lin =   replab.sym.SymPoly.fromCell(compositions);
                
                function arr = addElem(arr,ind,value)
                    arr(ind) = arr(ind) + value;
                    arr = sort(arr,'descend');
                end
            end

           function bilin = applyBilinear(lin1,lin2,func)
                allData = replab.sym.CTData.instance;
                coeffs3 = kron(lin1.coeffs,lin2.coeffs);
                data1 = allData(lin1.m);
                data2 = allData(lin2.m);
                data3 = allData(lin1.m+lin2.m);
                parts1 = data1.partitions;
                parts2 = data2.partitions;
                len1 = data1.nParts;
                len2 = data2.nParts;
                parts3 = zeros(len1*len2,data3.nParts);
                for j = 1:len1
                    for k = 1:len2
                        if coeffs3((j-1)*len2+k)
                            lin = func(parts1.list{j}.partition,parts2.list{k}.partition);
                            parts3((j-1)*len2+k,:) = lin.coeffs;
                        end
                    end
                end
                bilin =   replab.sym.SymPoly(coeffs3*parts3,lin1.m+lin2.m);
            end

           function augMonom = powerToAugmentedMonomial(part)
                len = numel(part);
                if len == 1
                    augMonom =    replab.sym.SymPoly.fromCell({part});
                    return
                end
                augMonom = replab.sym.SymPoly.applyBilinear( replab.sym.SymPoly.fromCell({part(1)}), ...
                replab.sym.SymPoly.fromCell({part(2)}),@replab.sym.SymPoly.multRowWithAugmented);
                for i = 3:len
                    lin2 =    replab.sym.SymPoly.fromCell({part(i)});
                    augMonom = replab.sym.SymPoly.applyBilinear(lin2,augMonom,@replab.sym.SymPoly.multRowWithAugmented);
                end
            end

           function monom = powerToMonomial(part,n)
                allData = replab.sym.CTData.instance;
                nData = allData(n);
                augMonom = replab.sym.SymPoly.powerToAugmentedMonomial(part);
                monom = augMonom.multDiag(nData.stabSizes);
           end
    end

    methods
        function self = SymPoly(arr,m)
            self.coeffs = arr;
            self.m = m;
        end

        function self = multDiag(self,diag)
            self.coeffs = self.coeffs.*diag;
        end


    end


end