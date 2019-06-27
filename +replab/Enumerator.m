classdef Enumerator < replab.Str
    
    properties (SetAccess = protected)
        size; % Number of elements contained in this enumerator (vpi)
    end
    
    methods (Access = protected)
        
        function s = describe(self, i)
        % Returns a one line string that describes the i-th element
        % used by Enumerator.str
            sizeNChars = length(strtrim(num2str(self.size)));
            msg = ['at(%' num2str(sizeNChars) 's) = %s'];
            el = replab.strOf(self.at(i));
            s = sprintf(msg, strtrim(num2str(i)), el);
        end
        
    end
    
    methods % Abstract
        
        function obj = at(self, ind)
        % Returns the element at the "ind" position
        % ind is an integer encoded as a double, string or vpi object
            obj = self.atFun(vpi(ind));
        end
        
        function ind = find(self, obj)
        % Returns the position of the given element as a vpi object
        % or [] if the element cannot be found
            ind = self.findFun(obj);
        end
        
    end
    
    methods
        
        function self = Enumerator(size)
            self.size = size;
        end
        
        function s = shortStr(self, maxColumns)
            s = sprintf('Enumerator of %s elements', replab.shortStr(self.size, maxColumns));
            if length(s) > maxColumns
                s = 'Enumerator';
            end
        end
        
        function [lines overLimit] = longStr(self, maxRows, maxColumns)
            if self.size > maxRows - 1
                n = maxRows - 2;
                start = ceil(n/2);
                finish = n - start;
                omit = self.size - start - finish;
                omitting = 1;
            else
                start = self.size;
                omit = [];
                omitting = 0;
                finish = 0;
            end
            table = cell(start + omitting + finish, 3);
            ind = 1;
            for i = 1:start
                % TODO: we do not register the overLimit from inner printing
                table{ind,1} = replab.shortStr(i, maxColumns);
                table{ind,2} = ' = ';
                table{ind,3} = replab.shortStr(self.at(i), maxColumns);
                ind = ind + 1;
            end
            if omitting
                table{ind,1} = ['.. ' replab.shortStr(omit, maxColumns)];
                table{ind,2} = '';
                table{ind,3} = 'elements omitted';
                ind = ind + 1;
            end
            for i = 1:finish
                index = self.size - finish + i;
                table{ind,1} = replab.shortStr(index, maxColumns);
                table{ind,2} = ' = ';
                table{ind,3} = replab.shortStr(self.at(index), maxColumns);
                ind = ind + 1;
            end
            lines = vertcat({self.shortStr(maxColumns)}, replab.str.align(table, 'rcl'));
            [lines overLimit] = replab.str.longFit(lines, false, maxRows, maxColumns);
        end
        
        function obj = sample(self)
        % Returns an element uniformly sampled from this Enumerator
            obj = self.at(randint(self.size));
        end
        
        function C = toCell(self)
            if self.size == 0
                C = cell(1, 0);
            else
                n = self.size;
                msg = 'Enumerator of size %s too big to enumerate in a matrix';
                assert(n < intmax('int32'), msg, num2str(self.size));
                n = double(n);
                C = cell(1, n);
                for i = 1:n
                    C{i} = self.at(i);
                end
            end
        end
        
    end
    
end
