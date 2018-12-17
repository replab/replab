classdef DomainFun < replab.cat.Domain
    
    properties (SetAccess = protected)
        description;
        eqvFun;
        hashFun;
        sampleFun;
    end
    
    methods
        
        function self = DomainFun(description, eqvFun, hashFun, sampleFun)
            self.parentOption = [];
            self.description = description;
            self.eqvFun = eqvFun;
            self.sampleFun = sampleFun;
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
