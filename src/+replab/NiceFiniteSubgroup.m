classdef NiceFiniteSubgroup < replab.NiceFiniteGroup
% A generic implementation of a subgroup of a nice finite group

    methods

        function self = NiceFiniteSubgroup(parent, generators, order)
            self.parent = parent;
            self.identity = parent.identity;
            % own stuff
            if nargin > 2 && ~isempty(order)
                self.order_ = order;
            end
            for i = 1:length(generators)
                assert(~parent.isIdentity(generators{i}), 'Generator cannot be identity');
            end
            self.generators = generators;
        end

    end

end
