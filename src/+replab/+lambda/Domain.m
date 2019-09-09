classdef Domain < replab.Domain
    
    properties (SetAccess = protected)
        header;
        eqvFun;
        sampleFun;
    end
    
    methods
        
        function self = Domain(header, eqvFun, sampleFun)
            self.header = header;
            self.eqvFun = eqvFun;
            self.sampleFun = sampleFun;
        end
        
        function str = headerStr(self)
            str = self.header;
        end
        
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Domain(self), ...
                {'header'} ...
                );
        end
    
    end
    
end
