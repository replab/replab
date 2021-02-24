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

end
