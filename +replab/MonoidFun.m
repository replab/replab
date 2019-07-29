classdef MonoidFun < replab.DomainFun & replab.Monoid
    
    properties (SetAccess = protected)
        composeFun;
    end

    methods
        
        function self = MonoidFun(description, eqvFun, sampleFun, ... % Domain
                                     composeFun, identity) % Monoid
            self@replab.DomainFun(description, eqvFun, sampleFun);
            self.composeFun = composeFun;
            self.identity = identity;
        end
        
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.DomainFun(self), ...
                hiddenFields@replab.Monoid(self), ...
                {'composeFun'} ...
                );
        end
    end
end
