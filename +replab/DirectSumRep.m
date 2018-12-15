classdef DirectSumRep < replab.Rep
    
    properties (SetAccess = protected)
        children;
    end
    
    methods
        
        function self = DirectSumRep(group, children, field)
            self.group = group;
            self.d = sum(cellfun(@(x) x.d, children));
            self.children = children;
            self.field = field;
        end
        
        function prettyDisp(self, spaces)
            disp(sprintf('%s Direct sum of total dimension %d', spaces, self.d));
            for i = 1:length(self.children)
                self.children{i}.prettyDisp([spaces '  ']);
            end
        end
       
        function M = image(self, permutation)
            blocks = cell(1, length(self.children));
            for i = 1:length(self.children)
                blocks{i} = self.children{i}.image(permutation);
            end
            M = blkdiag(blocks{:});
        end
                
    end
    
end
