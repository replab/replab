classdef DomainFun < replab.cat.Domain
    
    properties (SetAccess = protected)
        canEqv;
        canHash;
        canSample;
        description;
        eqvFun;
        hashFun;
        sampleFun;
    end
    
    methods
        
        function self = DomainFun(description, eqvFun, hashFun, sampleFun)
            self.description = description;
            self.eqvFun = eqvFun;
            self.hashFun = hashFun;
            self.sampleFun = sampleFun;
            self.canEqv = ~isempty(eqvFun);
            self.canHash = ~isempty(hashFun);
            self.canSample = ~isempty(sampleFun);
        end
        
        function b = eqv(self, x, y)
            b = self.eqvFun(x, y);
        end
        
        function h = hash(self, x)
            h = self.hashFun(x);
        end
        
        function x = sample(self)
            x = self.sampleFun();
        end

    end
    
end
