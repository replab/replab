classdef PermutationBSGSGroup < replab.PermutationGroup & replab.BSGSGroup
    methods
        function self = PermutationBSGSGroup(parent, generators, orderOpt)
            if nargin < 3
                orderOpt = [];
            end
            self = self@replab.PermutationGroup(parent.domainSize);
            self.action = parent.naturalAction;
            self.generators = generators;
            self.order_ = orderOpt;
        end
    end
end
