classdef SemidirectProductOfCompactGroups < replab.SemidirectProductGroup & replab.CompactGroup

    methods

        function self = SemidirectProductOfCompactGroups(phi)
            assert(isa(phi.G, 'replab.CompactGroup'));
            assert(isa(phi.P, 'replab.CompactGroup'));
            self@replab.SemidirectProductGroup(phi);
        end

    end

end
