classdef Enumerator < replab.Str
    
    properties
        D; % Domain in which the elements are contained
        size; % Number of elements contained in this enumerator
        sizeNChars; % Number of chars required to represent any index as a string
        atFun; % Handle that implements Enumerator.at
        findFun; % Handle that implements Enumerator.find
    end
    
    methods (Access = protected)
        
        function s = describe(self, i)
        % Returns a one line string that describes the i-th element
        % used by Enumerator.str
            msg = ['at(%' num2str(self.sizeNChars) 's) = %s'];
            el = replab.strOf(self.at(i));
            s = sprintf(msg, strtrim(num2str(i)), el);
        end
        
    end
    
    methods
        
        function self = Enumerator(D, size, atFun, findFun)
            self.D = D;
            self.size = size;
            self.sizeNChars = length(strtrim(num2str(self.size)));
            self.atFun = atFun;
            self.findFun = findFun;
        end
        
        function s = str(self)
            s = sprintf('Enumerator of %s elements', strtrim(num2str(self.size)));
            if self.size < 10
                for i = 1:double(self.size)
                    s = [s char(10) self.describe(i)];
                end
            else
                for i = 1:3
                    s = [s char(10) self.describe(i)];
                end        
                s = [s char(10) '..' strtrim(num2str(self.size - 6)) ' elements omitted..'];
                for i = 2:-1:0
                    s = [s char(10) self.describe(self.size - i)];
                end
            end
        end            
        
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
        
        function obj = sample(self)
        % Returns an element uniformly sampled from this Enumerator
            obj = self.at(randint(self.size));
        end
        
        function M = toMatrix(self)
            if self.size == 0
                M = [];
            else
                msg = 'Enumerator of size %s too big to enumerate in a matrix';
                assert(self.size < intmax('int32'), msg, num2str(self.size));
                size = double(self.size);
                x1 = self.at(1);
                msg = 'Elements must be vectors';
                assert(isvector(x1) && isnumeric(x1), msg);
                d = length(x1);
                M = zeros(size, d, class(x1));
                M(1, :) = x1;
                for i = 2:size
                    M(i, :) = self.at(i);
                end
            end
        end
        
    end
    
end