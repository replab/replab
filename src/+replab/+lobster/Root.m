classdef Root < replab.lobster.Node
  
    methods
        function self = Root()
            self@replab.lobster.Node('');
        end
        
        function str = render(self, context)
            if ~exist('context', 'var')
                context = struct();
            end
            
            str = self.render_children(context);
        end        
        
        function process_fragment(self, fragment)
        end
        
    end
    
end
