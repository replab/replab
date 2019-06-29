classdef StrFun < replab.Str
    properties (SetAccess = protected)
        shortStrFun;
        longStrFun;
    end
    
    methods
        
        function self = StrFun(shortStrFun, longStrFun)
            if nargin < 1
                shortStrFun = [];
            end
            if nargin < 2
                longStrFun = [];
            end
            self.shortStrFun = shortStrFun;
            self.longStrFun = longStrFun;
        end
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Str(self);
            names{end+1} = 'shortStrFun';
            names{end+1} = 'longStrFun';
        end
            
    end
    
end
