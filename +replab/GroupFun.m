classdef GroupFun < replab.MonoidFun & replab.Group
    properties (SetAccess = protected)
        inverseFun;
    end
    methods
        function self = GroupFun(description, eqvFun, sampleFun, ... % Domain
                                 composeFun, ... % Semigroup
                                 identity, ... % Monoid
                                 inverseFun) % Group
            self = self@replab.MonoidFun(description, eqvFun, sampleFun, composeFun, identity);
            self.inverseFun = inverseFun;
        end
    end    
end
