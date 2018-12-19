classdef Dec
    
    properties (Abstract, SetAccess = immutable)
        algebra;
        fromBlock;
        U;
        nComponents;
        compDims;
        repDims;
        repMuls;
    end
        
end
