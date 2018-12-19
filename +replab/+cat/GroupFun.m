classdef GroupFun < replab.cat.Group
    
    properties (SetAccess = protected)
        description;
        canEqv;
        canHash;
        canSample;
        eqvFun;
        hashFun;
        sampleFun;
        composeFun;
        inverseFun;
    end
    
    methods
        
        function self = GroupFun(description, eqvFun, hashFun, sampleFun, composeFun, identity, inverseFun)
            self.description = description;
            self.canEqv = ~isempty(eqvFun);
            self.canHash = ~isempty(hashFun);
            self.canSample = ~isempty(sampleFun);
            self.eqvFun = eqvFun;
            self.hashFun = hashFun;
            self.sampleFun = sampleFun;
            self.description = description;
            self.composeFun = composeFun;
            self.identity = identity;
            self.inverseFun = inverseFun;
        end
        
    end    
    
end
