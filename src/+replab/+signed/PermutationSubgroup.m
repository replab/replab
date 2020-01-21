classdef PermutationSubgroup < replab.signed.PermutationGroup
% Represents a subgroup of a signed permutation group

    methods

        function self = PermutationSubgroup(parent, generators, order)
            % stuff derived from parent
            self.parent = parent;
            self.domainSize = parent.domainSize;
            self.identity = parent.identity;
            % own stuff
            if nargin > 2 && ~isempty(order)
                self.order_ = order;
            end
            self.generators = generators;
            for i = 1:length(generators)
                assert(~parent.isIdentity(generators{i}), 'Generator cannot be identity');
            end
        end

    end

end
