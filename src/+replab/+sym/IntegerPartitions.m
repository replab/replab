classdef IntegerPartitions
    properties
        list % cell of replab.sym.IntegerPartition
        set % replab.perm.Set
        nParts% double: number of partitions
        n % double
    end
    
    methods (Static,Access = protected)
        function padded = padded(partition,n)
            padded = [partition zeros(1,n-numel(partition))];
        end
    end
    methods
        function self = IntegerPartitions(n)
            self.n = n;
            parts = replab.sym.IntegerPartition.enumerate(n);
            self.nParts = numel(parts);
            self.list = cellfun(@(part) replab.sym.IntegerPartition(part,n),parts,'UniformOutput',0);
            self.set = replab.perm.Set(n);
            for i =1:self.nParts
                self.set.insert(self.padded(parts{i},n)');
            end
        end
        
        function ind = index(self,partition)
            ind = self.set.find(self.padded(partition,self.n)');
        end
    end
end