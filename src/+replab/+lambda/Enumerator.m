classdef Enumerator < replab.Enumerator
    
    properties (SetAccess = protected)
        atFun; % Handle that implements Enumerator.at
        findFun; % Handle that implements Enumerator.find
    end
    
    methods
        
        function self = Enumerator(size, atFun, findFun)
            if isa(size, 'vpi')
                self.size = size;
            else
                self.size = vpi(size);
            end
            self.atFun = atFun;
            self.findFun = findFun;            
        end
        
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Enumerator(self), ...
                {'atFun' 'sampleFun'} ...
                );
        end
        
        function obj = at(self, ind)
            obj = self.atFun(vpi(ind));
        end
        
        function ind = find(self, obj)
            ind = self.findFun(obj);
        end
        
    end
    
end
