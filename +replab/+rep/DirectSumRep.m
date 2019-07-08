classdef DirectSumRep < replab.Rep
    
    properties
        blocks;
    end
    
    methods
        
        function self = DirectSumRep(blocks)
            assert(length(blocks) >= 1);
            d = 0;
            for i = 1:length(blocks)
                assert(isa(blocks{i}, 'replab.Rep'));
                d = d + blocks{i}.dimension;
            end
            self.dimension = d;
            for i = 2:length(blocks)
                assert(isequal(blocks{1}.group, blocks{i}.group));
                assert(isequal(blocks{1}.field, blocks{i}.field));
            end
            self.blocks = blocks;
            self.group = blocks{1}.group;
            self.field = blocks{1}.field;
        end
        
        function n = nBlocks(self)
            n = length(self.blocks);
        end
        
        function block = block(self, i)
            block = self.blocks{i};
        end
                
        % Str
                
        function names = hiddenFields(self)
            names = hiddenFields@replab.Rep(self);
            names{1, end+1} = 'blocks';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            for i = 1:self.nBlocks
                names{1, end+1} = sprintf('block(%d)', i);
                values{1, end+1} = self.block(i);
            end
        end
        
        % Rep

        function rho = image(self, g)
            rhos = cellfun(@(rep) rep.image(g), self.blocks, 'uniform', 0);
            rho = blkdiag(rho{:});
        end
        
    end
end
