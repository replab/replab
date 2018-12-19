classdef GroupMorphismFun < replab.cat.GroupMorphism
    
    properties
        imageFun; % morphism
    end
    
    methods
        
        function self = GroupMorphismFun(imageFun, S, T)
            self.imageFun = imageFun;
            self.S = S;
            self.T = T;
        end
            
    end
    
end

