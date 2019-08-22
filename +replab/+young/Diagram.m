classdef Diagram < replab.Str
% Young diagram
%
% A Young diagram is described by a partition $L = (l_1, ..., l_k)$
% of a nonnegative integer $n = l_1 + ... + l_k$, such that
% $l_{i+1} \le l_i$.
    
    properties
        partition % Partition of n
    end
    
    methods
        
        function self = Diagram(partition)
            self.partition = partition;
        end
        
        function str = shortStr(self, maxColumns)
            s1 = cellfun(@num2str, self.partition, 'uniform', 0);
            s2 = strjoin(s1, '|');
            str = ['Young diagram corresponding to partition: ' s2];
        end
        
        function str = longStr(self, maxRows, maxColumns)
            p = self.partition;
            n = length(p);
            if n <= maxRows
                rows = cell(n, 1);
                for i = 1:n
                    s1 = cell(1, p(i));
                    s1(:) = {' '};
                    s2 = strjoin(s1, '|');
                    rows{i} = ['[' s2 ']'];
                end
                str = rows;
            else
                str = {self.shortStr(maxColumns)};
            end
        end
        
        function str = asIdentifier(self)
            numbers = arrayfun(@num2str, self.partition, 'uniform', 0);
            str = ['Y' strjoin(numbers, '_')];
        end
        
    end
    
end
