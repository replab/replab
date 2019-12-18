classdef ElseNode < replab.lobster.Node
   
    methods
        
        function self = ElseNode()
            self@replab.lobster.Node('');
        end
        
        function str = render(~, ~)
            str = '';
        end
        
        function process_fragment(self, fragment)
        end
        
    end
    
end
