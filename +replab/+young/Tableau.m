classdef Tableau < replab.Str
    
    properties
        rows;
    end
    
    methods

        function self = Tableau(rows)
            self.rows = rows;
        end
        
        function str = headerStr(self, maxColumns)
            l = cellfun(@length, self.rows);
            str = sprintf('Young tableau of %d elements, %d lines', sum(l), length(l));
        end
        
        function str = shortStr(self, maxColumns)
            n = length(self.rows);
            rs = cell(1, n);
            for i = 1:n
                r1 = cellfun(@num2str, self.rows{i}, 'uniform', 0);
                rs{i} = strjoin(r1, ' ');
            end
            r = strjoin(rs, ';');
            str = ['Young diagram with rows: ' r];
            if length(str) > maxColumns
                str = self.headerStr(self, maxColumns);
            end
        end
        
        function str = longStr(self, maxRows, maxColumns)
            n = length(self.rows);
            if n <= maxRows
                l = cellfun(@length, self.rows);
                n = length(l);
                table = cell(n, 2*max(l)+1);
                spec = repmat('l', 1, 2*max(l)+1);
                for i = 1:n
                    table{i,1} = '[';
                    row = self.rows{i};
                    for j = 1:length(row)
                        table{i,2*j} = num2str(row(j));
                        if j < length(row)
                            table{i,2*j+1} = '|';
                        else
                            table{i,2*j+1} = ']';
                        end
                    end
                end
                str = replab.str.align(table, spec);
            else
                str = {self.shortStr(self, maxColumns)};
            end
        end
        
        function D = diagram(self)
            partition = cellfun(@(x) length(x), self.rows);
            D = replab.young.Diagram(partition);
        end
        
    end
    
end
