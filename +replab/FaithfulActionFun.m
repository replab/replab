classdef FaithfulActionFun < replab.ActionFun & replab.FaithfulAction
    
    properties (SetAccess = protected)
        findMovedElementFun;
    end
    
    methods
        
        function self = FaithfulActionFun(description, G, P, leftActionFun, findMovedElementFun)
            self@replab.ActionFun(description, G, P, leftActionFun);
            self.findMovedElementFun = findMovedElementFun;
        end
        
        function names = hiddenFields(self)
            names1 = hiddenFields@replab.ActionFun(self);
            names2 = hiddenFields@replab.FaithfulAction(self);
            names = vertcat(names1(:), names2(:));
            names{end+1, 1} = 'findMovedElementFun';
            names = unique(names);
        end

    end
    
end
