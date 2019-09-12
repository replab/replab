classdef SignedPermutationSubgroup < replab.SignedPermutationGroup
% Represents a subgroup of a signed permutation group
    
    methods

        function self = SignedPermutationSubgroup(parent, generators, order)
            % stuff derived from parent
            self.parent = parent;
            self.domainSize = parent.domainSize;
            self.identity = parent.identity;
            % own stuff
            if nargin > 2 && ~isempty(order)
                self.order_ = order;
            end
            self.generators = generators;
        end
        
    end
        
end
