function [conjClasses partCell] = findConjClasses(group)
    pData = replab.sym.findPartitions(group.domainSize);
    conjClasses = cellfun(@(part) group.conjugacyClass(conjClassElem(part)),pData.partCell,'UniformOutput', false);
    partCell = pData.partCell;
    function arr = conjClassElem(part)
        count = 1;
        arr = 1:sum(part);
        for p = 1:length(part)
            arr(count:count+part(p)-1) = circshift(arr(count:count+part(p)-1),1);
            count = count + part(p);
        end
    end
end