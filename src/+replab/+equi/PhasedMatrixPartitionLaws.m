classdef PhasedMatrixPartitionLaws < replab.Laws
% Law checks for a PhasedMatrixPartition

    properties (SetAccess = protected)
        P % (`.PhasedMatrixPartition`): Phased matrix partition tested
    end

    methods

        function self = PhasedMatrixPartitionLaws(P)
            self.P = P;
        end

        function law_phases_are_consistent_(self)
            self.assert(all(all(self.P.phase >= 0)));
            self.assert(all(all(self.P.phase < self.P.phaseOrder)));
        end

        function law_no_empty_block_(self)
            L = cellfun(@(blk) size(blk, 2), self.P.blocks);
            self.assert(all(L > 0));
        end

        function law_zero_block_is_consistent_(self)
            ind1 = self.P.zeroBlock';
            [I, J] = find(self.P.block == 0);
            ind2 = [I(:) J(:)];
            self.assert(isequal(sortrows(ind1), sortrows(ind2)));
        end

        function law_blocks_are_consistent_(self)
            for i = 1:self.P.nBlocks
                [I, J] = find(self.P.block == i);
                ind1 = self.P.blocks{i}';
                ind2 = [I(:) J(:)];
                self.assert(isequal(sortrows(ind1), sortrows(ind2)));
            end
        end

    end

end
