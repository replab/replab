classdef PermutationSubgroup < replab.PermutationGroup
% Represents a subgroup of a permutation group
    
    methods
        
        function self = PermutationSubgroup(parent, generators, order)
            % stuff derived from parent
            self.parent = parent;
            self.domainSize = parent.domainSize;
            self.identity = parent.identity;
            self.niceMonomorphism = parent.niceMonomorphism;
            % own stuff
            if nargin > 2 && ~isempty(order)
                self.order_ = order;
            end
            self.generators = generators;
        end
        
    end
        
end
