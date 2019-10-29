classdef DirectSumRep < replab.Rep
% A direct sum of representations, such that images are diagonal by blocks
    
    properties
        blocks % row cell array of replab.Rep: Contained subrepresentations
    end
    
    methods
        
        function self = DirectSumRep(blocks)
        % Constructs a direct sum from a cell array of representations
        %
        % All the subrepresentations should be defined on the same group, and on the same field.
        %
        % Args:
        %   blocks (row cell array of replab.Rep): Subrepresentations
            assert(length(blocks) >= 1);
            d = 0;
            for i = 1:length(blocks)
                assert(isa(blocks{i}, 'replab.Rep'));
                d = d + blocks{i}.dimension;
            end
            self.dimension = d;
            blocksAreUnitary = cellfun(@(x) x.isUnitary, blocks, 'uniform', 0);
            self.isUnitary = replab.trileanAnd(blocksAreUnitary{:});
            for i = 2:length(blocks)
                assert(blocks{1}.group == blocks{i}.group);
                assert(isequal(blocks{1}.field, blocks{i}.field));
            end
            self.blocks = blocks;
            self.group = blocks{1}.group;
            self.field = blocks{1}.field;
        end
        
        function n = nBlocks(self)
        % Returns the number of blocks in the direct sum
        %
        % Returns:
        %   integer: Number of subrepresentations composing this representation
            n = length(self.blocks);
        end
        
        function block = block(self, i)
        % Returns a block in the direct sum
        %
        % Args:
        %   i (integer): Index of block
        %
        % Returns:
        %   replab.Rep: Representation corresponding to the i-th block
            block = self.blocks{i};
        end

        %% Str methods
                
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
        
        %% Rep methods

        function rho = image(self, g)
            rhos = cellfun(@(rep) rep.image(g), self.blocks, 'uniform', 0);
            rho = blkdiag(rhos{:});
        end
        
        function rho = inverseImage(self, g)
            rhos = cellfun(@(rep) rep.inverseImage(g), self.blocks, 'uniform', 0);
            rho = blkdiag(rhos{:});
        end

    end

end
