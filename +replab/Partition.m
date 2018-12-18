classdef Partition
   
    properties
        n;      % domain size
        blockIndex;  % index of the block for each element
        start;  % starting index for each block
        next;   % next index in the same block, or 0 if at the end
    end
    
    methods (Access = protected)
        
        function self = Partition(n, blockIndex, start, next)
            self.n = n;
            self.blockIndex = blockIndex;
            self.start = start;
            self.next = next;
        end
        
    end
    
    methods
        
        function n = nBlocks(self)
            n = length(self.start);
        end
        
        function B = block(self, i)
            el = self.start(i);
            B = [];
            while el > 0
                B = [B el];
                el = self.next(el);
            end
        end
        
    end

    methods (Static)
        
        function P = permutationsOrbits(permutationsAsMatrix)
        % permutations are a nG x domainSize double matrix
            n = size(permutationsAsMatrix, 2);
            nG = size(permutationsAsMatrix, 1);
            rest = 1:n;
            blockIndex = zeros(1, n);
            start = [];
            next = zeros(1, n);
            block = 1;
            for i = 1:n
                if blockIndex(i) == 0
                    % new block discovered, starting with i
                    start = [start i];
                    added = i;
                    blockIndex(i) = block;
                    test = i;
                    while ~isempty(test)
                        t = test(1);
                        test = test(2:end);
                        for j = 1:nG
                            ti = permutationsAsMatrix(j, t);
                            if blockIndex(ti) == 0
                                added = [added ti];
                                blockIndex(ti) = block;
                                test = [test ti];
                            end
                        end
                    end
                    added = sort(added);
                    for i = 1:length(added) - 1
                        next(added(i)) = added(i + 1);
                    end
                    block = block + 1;
                end
            end
            P = replab.Partition(n, blockIndex, start, next);
        end
                
    end

end
