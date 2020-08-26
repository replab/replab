classdef CTData
    
    properties
        n
        fact %factorial
        partitionHash %canonical ordering
        partitionList %list of partitions
        stabSizes %stabilizer sizes list
        conjSizes % conjugacy class sizes
        cycSizes %'Group notation' Cycle sizes list
        nCycles %'Group notation' Number sizes list
        nParts
    end
    
    methods
        function self = CTData(n)
        self.n=n;
        %Factorial
        self.fact = factorial(n);
        genPartitions(n);
                function genPartitions(n) 
                    parts = replab.sym.findPartitions(n);
                    [self.partitionList,self.partitionHash,self.nParts] = deal(parts.partCell,parts.partitionHash,parts.nParts);
                    self.cycSizes = parts.cycleSizes;
                    self.nCycles = parts.nCycles;
                    self.stabSizes = zeros(1,self.nParts);
                    self.conjSizes = zeros(1,self.nParts);
                    for i = 1:self.nParts
                        part = self.partitionList{i};
                        cell(1,self.nParts);
                        self.stabSizes(i) = prod(factorial(self.nCycles{i}));
                        self.conjSizes(i)= self.fact/prod(part)/self.stabSizes(i);
                    end
                end
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






