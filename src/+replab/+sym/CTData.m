classdef CTData
    
    properties
        n
        fact %factorial
        partitions %replab.sym.IntegerPartitions
        stabSizes %double(1,\*): stabilizer sizes list
        conjSizes %double(1,\*): conjugacy class sizes
        innerProdDenom % %double(1,\*): Inner product is <y,x> = y'*diag(1./ innerProdDenom)*x
        nParts
    end
    
    methods
        function self = CTData(n)
        self.n=n;
        %Factorial
        self.fact = factorial(n);
        self.partitions = replab.sym.IntegerPartitions(n);
        self.nParts = self.partitions.nParts;
        self.stabSizes = cellfun(@(part) prod(factorial(nonzeros(part.powers))),self.partitions.list);
        self.innerProdDenom = round(self.stabSizes.*cellfun(@(part) prod(part.partition),self.partitions.list));
        self.conjSizes = round(self.fact./self.innerProdDenom);
        end
    end
    
    methods(Static)
        function nData = instance(n)
            persistent instance_
            persistent m
            if isempty(m)
                m = 1;
            end
            if nargin == 0
                newInstance = false;
            elseif nargin == 1
                newInstance = (m < n)|| n==1;
            end
            if newInstance
                if isempty(instance_)
                    instance_ = replab.sym.CTData.empty;
                end
                for i = m:n
                    instance_(i) = replab.sym.CTData(i);
                end
                m = n;
            end
            nData = instance_;
        end
     end
        

end







