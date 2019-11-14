classdef GenericShape < replab.tensor.Shape

    properties (SetAccess = immutable)
        fullIndices % integer matrix: Matrix of all (unsymmetrized) subindices
        partition % `+replab.Partition`: Partition of the unsymmetrized indices under action of the group
    end
    
    methods
       
        function self = GenericShape(dimensions, group, isOrderColumnMajor)
            self = self@replab.tensor.Shape(dimensions, group, isOrderColumnMajor);
            if isOrderColumnMajor
                fullIndices = replab.tensor.GenericShape.colMajorEnumerate(dimensions);
            else
                fullIndices = replab.tensor.GenericShape.rowMajorEnumerate(dimensions);
            end
            n = length(dimensions);
            nf = size(fullIndices, 1);
            ng = group.nGenerators;
            J = zeros(ng, nf);
            sub = zeros(1, n);
            for i = 1:group.nGenerators
                g = group.generator(i);
                for j = 1:size(fullIndices, 1)
                    sub(g) = fullIndices(j, :);
                    [~, res] = ismember(sub, fullIndices, 'rows');
                    assert(length(res) == 1);
                    J(i, j) = res;
                end
            end
            self.fullIndices = fullIndices;
            self.partition = replab.Partition.permutationsOrbits(J);
        end
        
        function n = nComponents(self)
            n = self.partition.nBlocks;
        end
        
        function sub = indToSub(self, ind)
            sub = self.fullIndices(self.partition.start(ind), :);
        end
        
        function ind = subToInd(self, sub)
            [~, fullInd] = ismember(sub, self.fullIndices);
            ind = self.partition.blockIndex(fullInd);
        end
        
        function orbit = subOrbit(self, sub)
            [~, fullInd] = ismember(sub, self.fullIndices);
            block = self.partition.block(self.partition.blockIndex(fullInd));
            orbit = self.fullIndices(block, :);
        end
        
    end
    
    methods (Static)
        
        function sub = rowMajorEnumerate(dimensions)
            switch length(dimensions)
              case 0
                  sub = zeros(1, 0);
              case 1
                sub = (1:dimensions)';
              otherwise
                d = dimensions(1);
                rest = replab.tensor.GenericShape.rowMajorEnumerate(dimensions(2:end));
                drest = size(rest, 1);
                sub = [kron((1:d)', ones(size(rest, 1), 1)) repmat(rest, d, 1)];
            end
        end
        
        function sub = colMajorEnumerate(dimensions)
            switch length(dimensions)
              case 0
                  sub = zeros(1, 0);
              case 1
                sub = (1:dimensions)';
              otherwise
                d = dimensions(end);
                rest = replab.tensor.GenericShape.colMajorEnumerate(dimensions(1:end-1));
                drest = size(rest, 1);
                sub = [repmat(rest, d, 1) kron((1:d)', ones(size(rest, 1), 1))];
            end
        end
        
    end
    
end
