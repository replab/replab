classdef Centralizer
% Computes the centralizer of a subgroup

    properties
        group % `+replab.PermutationGroup`: Group
        other % `+replab.PermutationGroup`: Subgroup to centralize
    end

    methods

        function self = Centralizer(group, other)
            self.group = group;
            self.other = other;
        end


        function s = subgroup(self)
            orbits = self.other.orbits.blocks;
            % remove singletons
            orbits = orbits(cellfun(@(o) length(o) > 1, orbits));
            [~, I] = sort(-cellfun(@length, orbits));
            orbits = orbits(I); % largest orbits first
            % now we have the blocks of size > 1
            base = [orbits{:}];
            L = length(base);
            chain = self.group.chain.mutableCopy;
            chain.baseChange(base);
            chain.makeImmutable;
            l = 1;
            tests = cell(1, L);
            for i = 1:length(orbits)
                orbit = orbits{i};
                % first element of an orbit stores the current element tested
                tests{l} = @(g, data) deal(true, g);
                l = l + 1;
                for j = 2:length(orbit)
                    b = orbit(j);
                    assert(b == base(l));
                    tests{l} = @(g, g0) deal(g0(b) == g(b), g0);
                    l = l + 1;
                end
            end
            c = replab.bsgs.subgroupSearch(chain, @(g) self.prop(g), tests);
            s = replab.PermutationGroup.fromChain(c, self.group.type);
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
