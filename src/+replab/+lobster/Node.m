classdef Node < handle 
       
    properties
       creates_scope = false
       children = cell(0)
    end
    
    methods 
        
        function self = Node(fragment)
            if ~exist('fragment', 'var')
                fragment = '';
            end
            
            self.process_fragment(fragment);
        end
        
        function process_fragment(self, fragment)
            error('Abstract');
        end
        
        function enter_scope(self)
            error('Abstract');
        end
        
        function exit_scope(self)
            error('Abstract');
        end
        
        function str = render(self, context)
            error('Abstract');
        end
        
        function str = render_children(self, context, children)
            if ~exist('children', 'var')
                children = self.children;
            end
            
            rendered_children = cellfun(@(x) x.render(context), ...
                children, 'Uniform', false);
            
            str = strjoin(rendered_children, '');
        end
            
        function add_child(self, child)
            if ~isa(child, 'replab.lobster.Node')
                error('Lobster:GenericError', ...
                    'Attempted to add a <%s> to the children of the root node', ...
                    class(child));
            end
            
           self.children{end+1} = child; 
        end
    end
    
    
end