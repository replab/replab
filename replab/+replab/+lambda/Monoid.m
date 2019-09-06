classdef Monoid < replab.Monoid
    
    properties (SetAccess = protected)
        header;
        eqvFun;
        sampleFun;
        composeFun;
    end

    methods
        
        function self = Monoid(header, eqvFun, sampleFun, ... % Domain
                               composeFun, identity) % Monoid
            self.header = header;
            self.eqvFun = eqvFun;
            self.composeFun = composeFun;
            self.identity = identity;
        end
        
        function str = headerStr(self)
            str = self.header;
        end
        
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Monoid(self), ...
                {'header'} ...
                );
        end
        
    end
    
end
