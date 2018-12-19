classdef BlockDiagonalMatrices < replab.cat.Domain

    properties (SetAccess = protected)
        blocks;
        canEqv = false;
        canHash = false;
        canSample = true;
    end
    
    methods
        
        function s = str(self)
            s = sprintf('Block diagonal matrix with %d blocks:\n', self.nBlocks);
            for i = 1:self.nBlocks
                s = sprintf('%s  %s\n', s, replab.strOf(self.blocks{i}));
            end
        end
        
        function self = BlockDiagonalMatrices(blocks)
            for i = 1:length(blocks)
                assert(isa(blocks{i}, 'replab.Matrices'));
                assert(blocks{i}.nRows == blocks{i}.nCols);
            end
            self.blocks = blocks;
        end
        
        function n = nBlocks(self)
            n = length(self.blocks);
        end
        
        function y = inverse(self, x)
            n = self.nBlocks;
            y = cell(1, n);
            for i = 1:n
                y{i} = inv(x{i});
            end
        end
        
        function z = compose(self, x, y)
            n = self.nBlocks;
            z = cell(1, n);
            for i = 1:n
                z{i} = x{i} * y{i};
            end
        end
        
        function x = sample(self)
            n = self.nBlocks;
            x = cell(1, n);
            for i = 1:n
                x{i} = self.blocks{i}.sample;
            end
        end
        
    end
    
end
