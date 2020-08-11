function conjClasses = findConjClasses(group)
    pData = replab.sym.findPartitions(group.domainSize);
    conjClasses = cellFun(@(part) group.conjugacyClass(conjClassElem(part),pData.partCells,'UniformOutput', false);
    function arr = conjClassElem(part)
        count = 1;
        arr = 1:sum(part);
        for p = 1:length(part)
            arr(count:count+part(p)-1) = circshift(arr(count:count+part(p)-1),1);
            count = count + part(p);
        end
    end
end