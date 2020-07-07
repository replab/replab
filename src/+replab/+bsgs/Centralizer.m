classdef Centralizer
% Computes the centralizer of an element

    properties
        group
        other
    end

    methods

        function self = Centralizer(group, other)
            self.group = group;
            self.other = other;
        end

        function s = subgroup(self)
            c = replab.bsgs.subgroupSearch(self.group.chain, @(g) self.prop(g));
            s = replab.PermutationGroup.fromChain(c, self.group.parent);
        end

        function b = prop(self, g)
            for i = 1:self.other.nGenerators
                gen = self.other.generator(i);
                if ~isequal(g(gen), gen(g))
                    b = false;
                    return
                end
            end
            b = true;
        end

    end

end
