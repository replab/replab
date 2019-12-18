classdef VarNode < replab.lobster.Node
    
    properties 
       name = ''; 
    end
    
    methods
        
        function self = VarNode(fragment)
           self@replab.lobster.Node(fragment); 
        end
        
        function process_fragment(self, fragment)
            self.name = fragment;
        end
        
        function str = render(self, context)
            var = replab.lobster.eval_with_context(self.name, context);
            if ischar(var) && isrow(var)
                str = var;
            else
                str = replab.shortStr(var);
            end
        end
        
    end
    
end
