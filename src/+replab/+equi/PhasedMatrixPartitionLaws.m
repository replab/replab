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

        function law_no_empty_subset_(self)
            L = cellfun(@(blk) size(ss, 2), self.P.subsets);
            self.assert(all(L > 0));
        end

        function law_zero_subset_is_consistent_(self)
            ind1 = self.P.zeroSubset';
            [I, J] = find(self.P.subsetIndex == 0);
            ind2 = [I(:) J(:)];
            self.assert(isequal(sortrows(ind1), sortrows(ind2)));
        end

        function law_subsets_are_consistent_(self)
            for i = 1:self.P.nSubsets
                [I, J] = find(self.P.subsetIndex == i);
                ind1 = self.P.subsets{i}';
                ind2 = [I(:) J(:)];
                self.assert(isequal(sortrows(ind1), sortrows(ind2)));
            end
        end

    end

end
