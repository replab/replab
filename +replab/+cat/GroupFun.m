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
            self.composeFun = composeFun;
            self.identity = identity;
            self.inverseFun = inverseFun;
        end
        
    end    
    
end
