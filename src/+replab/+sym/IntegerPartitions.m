classdef IntegerPartitions

    properties
        n % (integer): Integer to partition
        nParts % (integer): Number of partitions
        list % (cell(1,\*) of `.IntegerPartition`): List of integer partitions of `.n`
        set % (`+replab.+perm.Set`): Hash table of partitions for retrieval
    end

    properties % from CTData
        fact % (double): Factorial of `.n`
        stabSizes % (double(1,\*)): Stabilizer sizes list TODO
        conjSizes % (double(1,\*)): Conjugacy class sizes
        innerProdDenom % %double(1,\*): Inner product is ``<y,x> = y'*diag(1./ innerProdDenom)*x``
    end

    methods (Static)

        function ip = make(n)
            persistent cache
            if isempty(cache)
                cache = cell(1, 0);
            end
            if n+1 > length(cache) || isempty(cache{n+1})
                cache{1,n+1} = replab.sym.IntegerPartitions(n);
            end
            ip = cache{n+1};
        end

    end

    methods

        function self = IntegerPartitions(n)
            self.n = n;
            parts = replab.sym.IntegerPartition.enumerate(n);
            self.nParts = numel(parts);
            self.list = cellfun(@(part) replab.sym.IntegerPartition(part,n),parts,'UniformOutput',0);
            self.set = replab.perm.Set(n);
            for i = 1:self.nParts
                self.set.insert(self.padded(parts{i},n)');
            end
            self.fact = factorial(n);
            self.stabSizes = cellfun(@(part) prod(factorial(nonzeros(part.powers))), self.list);
            self.innerProdDenom = round(self.stabSizes.*cellfun(@(part) prod(part.partition), self.list));
            self.conjSizes = round(self.fact./self.innerProdDenom);
        end

        function ind = index(self, partition)
        % Finds an integer partition in the set
            ind = self.set.find(self.padded(partition,self.n)');
        end

        function conjClasses = conjugacyClasses(self)
            conjClasses = cellfun(@(intPart) intPart.conjugacyClass, self.list, 'UniformOutput', 0);
        end

        function names = conjugacyClassNames(self)
            names = cellfun(@(intPart) strrep(replab.shortStr(intPart.partition), ' ', ''), self.list, 'uniform', 0);
        end

    end

    methods (Static, Access = protected)

        function padded = padded(partition, n)
        % Pads an integer partition with zeros at the end
        %
        % Args:
        %   partition (integer(1,\*)): Integer partition
        %   n (integer): Integer being partitioned
        %
        % Returns:
        %   integer(1,\*): A padded integer row vector of length ``n`` with zeroes at the end
            padded = [partition zeros(1,n-numel(partition))];
        end

    end

end