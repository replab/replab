classdef DirectSumRep < replab.Rep
% A direct sum of representations, such that images are diagonal by blocks

    properties
        blocks % (cell(1,*) of `+replab.Rep`): Contained subrepresentations
    end

    methods

        function self = DirectSumRep(group, field, blocks)
        % Constructs a direct sum from a cell array of representations
        %
        % All the subrepresentations should be defined on the same group, and on the same field.
        %
        % Args:
        %   group (`+replab.CompactGroup`): Common group
        %   field ({'R', 'C'}): Real or complex field
        %   blocks (cell(1,*) of `+replab.Rep`): Subrepresentations
            replab.rep.assertCompatibleFactors(group, field, blocks);
            % own properties
            self.blocks = blocks;
            % replab.Rep immutable
            self.group = group;
            self.field = field;
            self.dimension = sum(cellfun(@(b) b.dimension, blocks));
            % replab.Rep mutable
            blocksAreUnitary = cellfun(@(x) x.isUnitary, blocks, 'uniform', 0);
            self.isUnitary = replab.trileanAnd(blocksAreUnitary{:});
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

        function rho = image_internal(self, g)
            rhos = cellfun(@(rep) rep.image_internal(g), self.blocks, 'uniform', 0);
            rho = blkdiag(rhos{:});
        end

        function rho = inverseImage(self, g)
            rhos = cellfun(@(rep) rep.inverseImage_internal(g), self.blocks, 'uniform', 0);
            rho = blkdiag(rhos{:});
        end

    end

end
