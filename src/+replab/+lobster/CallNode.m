classdef CallNode < replab.lobster.Node
    
    properties 
       expression = ''; 
    end
    
    methods
        function self = CallNode(fragment)
            self@replab.lobster.Node(fragment);
        end
        
        function process_fragment(self, fragment)
            self.expression = strtrim(fragment);
        end
        
        function str = render(self, context)
            str = replab.lobster.eval_with_context(self.expression, context);
            
            if ~ischar(str)
                error('Lobster:CallError', ...
                    'The output of <%s> was not a string.', self.expression);
            end
        end
    end
    
end