classdef PermutationBSGSGroup < replab.BSGSGroup & replab.PermutationGroup
    methods
        function self = PermutationBSGSGroup(parent, generators, orderOpt)
            if nargin < 3
                orderOpt = [];
            end
            self = self@replab.PermutationGroup(parent.domainSize);
            self = self@replab.BSGSGroup(parent, parent.naturalAction, generators, orderOpt);
        end
    end
end
