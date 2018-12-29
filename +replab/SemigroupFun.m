classdef SemigroupFun < replab.DomainFun & replab.Semigroup
    properties (SetAccess = protected)
        composeFun;
    end
    methods
        function self = SemigroupFun(description, eqvFun, sampleFun, ... % Domain
                                     composeFun) % Semigroup
            self = self@replab.DomainFun(description, eqvFun, sampleFun);
            self.composeFun = composeFun;
        end
    end    
end
