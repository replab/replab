classdef DirectProductGroup_compact < replab.DirectProductGroup & replab.CompactGroup

    methods

        function self = DirectProductGroup_compact(factors)
        % Constructs a direct product of groups
        %
        % Args:
        %   factors (cell(1,\*) of `.CompactGroup`): Factor groups
            assert(all(cellfun(@(x) isa(x, 'replab.CompactGroup'), factors)));
        end

    end

end
