classdef DirectProductGroup_compact < replab.DirectProductGroup & replab.CompactGroup

    methods

        function self = DirectProductGroup_compact(factors)
        % Constructs a direct product of groups
        %
        % Args:
        %   factors (cell(1,\*) of `.CompactGroup`): Factor groups
            assert(all(cellfun(@(x) isa(x, 'replab.CompactGroup'), factors)));
            self.factors = factors;
            self.identity = cellfun(@(f) f.identity, factors, 'uniform', 0);
        end

    end

    methods % Implementations

        function b = hasReconstruction(self)
            b = ~isempty(self.torusBlocks);
        end

        function [mu, R] = reconstruction(self)
            blocks = self.torusBlocks;
            N = sum(cellfun(@length, blocks));
            n = self.nFactors;
            f = @(t) arrayfun(@(i) self.factor(i).reconstruction.imageElement(t(blocks{i})), 1:n, 'uniform', 0);
            T = replab.T(N);
            sets = cell(1, 0);
            for i = 1:n
                Fi = self.factor(i);
                [~, Ri] = Fi.reconstruction;
                e = self.injection(i);
                sets = horzcat(sets, cellfun(@(S) cellfun(@(s) e.imageElement(s), S, 'uniform', 0), Ri.sets, 'uniform', 0));
            end
            R = replab.SetProduct(self, sets, true);
            mu = T.morphismByFunction(self, f);
        end

    end

end
