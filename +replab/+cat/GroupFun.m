classdef GroupFun < replab.cat.DomainFun & replab.cat.Group
    
    properties (SetAccess = protected)
        composeFun;
        inverseFun;
    end
    
    methods
        
        function self = GroupFun(description, eqvFun, hashFun, sampleFun, composeFun, identity, inverseFun)
            self@replab.cat.DomainFun(description, eqvFun, hashFun, sampleFun);
            self.G = self;
            self.N = replab.cat.Domain.integerRange(-10, 10);
            self.description = description;
            self.sampleFun = sampleFun;
            self.eqvFun = eqvFun;
            self.composeFun = composeFun;
            self.identity = identity;
            self.inverseFun = inverseFun;
        end
        
        function b = eqv(self, x, y)
            b = self.eqvFun(x, y);
        end

        function xInv = inverse(self, x)
            xInv = self.inverseFun(x);
        end
        
        function z = compose(self, x, y)
            z = self.composeFun(x, y);
        end
        
    end    
    
end
