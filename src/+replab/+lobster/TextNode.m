classdef TextNode < replab.lobster.Node
    
    properties
       text = ''; 
    end
    
    methods
        
        function self = TextNode(fragment)
            process_fragment(self, fragment);
        end
        
        function process_fragment(self, fragment)
            self.text = fragment;
        end
        
        function str = render(self, ~)
           str = self.text; 
       end
       
    end
    
end
