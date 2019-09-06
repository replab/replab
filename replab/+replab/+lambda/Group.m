classdef Group < replab.Group
    
    properties (SetAccess = protected)
        header;
        eqvFun;
        sampleFun;
        composeFun;        
        inverseFun;
    end
    
    methods
        
        function self = Group(header, eqvFun, sampleFun, ... % Domain
                              composeFun, identity, ... % Monoid
                              inverseFun) % Group
            self.header = header;
            self.eqvFun = eqvFun;
            self.composeFun = composeFun;
            self.identity = identity;
            self.inverseFun = inverseFun;
        end
        
        function str = headerStr(self)
            str = self.header;
        end

        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Group(self), ...
                {'header'} ...
                );
        end
        
    end
    
end
