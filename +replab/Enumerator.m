classdef Enumerator < replab.Str
    
    properties
        D;
        size;
        sizeLength;
        atFun;
        findFun;
    end
    
    methods (Access = protected)
        
        function s = describe(self, ind)
            msg = ['at(%' num2str(self.sizeLength) 's) = %s'];
            el = replab.strOf(self.at(ind), true);
            s = sprintf(msg, strtrim(num2str(ind)), el);
        end
        
    end
    
    methods
        
        function self = Enumerator(D, size, atFun, findFun)
            self.D = D;
            self.size = size;
            self.sizeLength = length(strtrim(num2str(self.size)));
            self.atFun = atFun;
            self.findFun = findFun;
        end
        
        function s = str(self)
            s = sprintf('Enumerator of %s elements', strtrim(num2str(self.size)));
            if self.size < 10
                for i = 1:double(self.size)
                    s = [s newline self.describe(i)];
                end
            else
                for i = 1:3
                    s = [s newline self.describe(i)];
                end        
                s = [s newline '..' strtrim(num2str(self.size - 6)) ' elements omitted..'];
                for i = 2:-1:0
                    s = [s newline self.describe(self.size - i)];
                end
            end
        end            
        
        function obj = at(self, ind)
            obj = self.atFun(ind);
        end
        
        function ind = find(self, obj)
            ind = self.findFun(obj);
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