classdef Enumerator < replab.Str
    
    properties (SetAccess = protected)
        size; % Number of elements contained in this enumerator (vpi)
    end
    
    methods (Access = protected)
        
        function s = describe(self, i)
        % Returns a one line string that describes the i-th element
        % used by Enumerator.str
            sizeNChars = length(strtrim(self.stringOnOneLine(num2str(self.size))));
            msg = ['at(%' num2str(sizeNChars) 's) = %s'];
            el = replab.strOf(self.at(i));
            s = sprintf(msg, strtrim(self.stringOnOneLine(num2str(i))), el);
        end
        
        function s = stringOnOneLine(self, s)
        % Helper function that returns a string possibly spanning several
        % lines onto just one line
            if (size(s,1) >= 2)
                s = reshape(s', 1, numel(s));
            end
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
        
        function s = str(self)
            s = sprintf('Enumerator of %s elements', strtrim(self.stringOnOneLine(num2str(self.size))));
            if self.size < 10
                for i = 1:double(self.size)
                    s = [s char(10) self.describe(i)];
                end
            else
                for i = 1:3
                    s = [s char(10) self.describe(i)];
                end        
                s = [s char(10) '..' strtrim(self.stringOnOneLine(num2str(self.size - 6))) ' elements omitted..'];
                for i = 2:-1:0
                    s = [s char(10) self.describe(self.size - i)];
                end
            end
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
