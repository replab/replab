classdef DirectProductOfCompactGroups < replab.DirectProductGroup & replab.CompactGroup

    methods

        function self = DirectProductOfCompactGroups(factors)
            assert(all(cellfun(@(x) isa(x, 'replab.CompactGroup'), factors)));
            self@replab.DirectProductGroup(factors);
        end

    end

end
