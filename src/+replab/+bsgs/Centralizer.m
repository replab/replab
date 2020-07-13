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
            n = self.group.domainSize;
            orbits = self.other.orbits.blocks;
            % compute the size of the orbit of each point
            orbitSizes = zeros(1, n);
            for i = 1:length(orbits)
                orbitSizes(orbits{i}) = length(orbits{i});
            end
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
                b = orbit(1);
                % first element of an orbit stores the current element tested
                % tests{l} = @(g, data) deal(true, g);
                tests{l} = @(g, data) deal(orbitSizes(b) == orbitSizes(g(b)), g(orbit));
                l = l + 1;
                for j = 2:length(orbit)
                    b = orbit(j);
                    assert(b == base(l));
                    %tests{l} = @(g, g0) deal(g0(b) == g(b), g0);
                    tests{l} = @(g, O) deal(true, O);
                    l = l + 1;
                end
            end
            tests = {};
            startData = [];
            initSubgroup = [];
            c = replab.bsgs.subgroupSearch(chain, @(g) self.prop(g), tests, startData, initSubgroup, true);
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

    methods (Static)

        function [chain base orbits allPoints orbitNumber orbitIndex] = orbits(originalChain, partition, largestFirst)
        % Arranges the base of the given chain so that the given orbits are contiguous
        %
        % Redundant base points are eliminated; however, they are still accounted for in the returned values.
        %
        % To the base point ``base(l)`` we associated the points ``allPoints{l}``, for which ``allPoints{l}(1) == base(l)``,
        % and ``allPoints{l}(2:end)`` correspond to the redundant base points after ``base(l)`` which have been eliminated.
        %
        % To ``allPoints{l}`` we associate ``orbitNumber{l}`` and ``orbitIndex{l}`` so that
        % ``orbits{orbitNumber{l}{i}}(orbitIndex{l}{i}) == allPoints{l}(i)``.

        % Args:
        %   chain (`+replab.bsgs.Chain`): Chain
        %   partition (`+replab.Partition`): Partition of the domain of the chain in orbits
        %   largestFirst (logical, optional): Whether to put the largest orbits first, default true
        %
        % Returns
        % -------
        %
        % chain:
        %   `+replab.bsgs.Chain`: Chain with reordered basis
        % base:
        %   `integer(1,\*)`: Prescribed basis, is a partial basis for ``chain``
        % orbits:
        %   `cell(1,\*) of integer(1,\*)`: Orbits present in the base
        % allPoints:
        %   `cell(1,\*) of integer(1,\*)`: Base points including redundant one, see description above
        % orbitNumber:
        %   `cell(1,\*) of integer(1,\*)`: Orbit numbers corresponding to ``allPoints``
        % orbitIndex:
        %   `cell(1,\*) of integer(1,\*)`: Index in the corresponding orbits for ``allPoints``
            if nargin < 3 || isempty(largestFirst)
                largestFirst = true;
            end
            blocks = partition.blocks;
            singletons = cellfun(@(block) length(block) == 1, blocks);
            blocks = blocks(~singletons);
            criterion = cellfun(@length, blocks);
            if largestFirst
                criterion = -criterion;
            end
            [~, ind] = sort(-cellfun(@length, criterion));
            orbits = blocks(I); % largest orbits first
            chain1 = chain.mutableCopy;
            base = [];
            allPoints = {};
            orbitNumber = {};
            orbitIndex = {};
            for i = 1:length(blocks)
                l = length(base) + 1;
                block = blocks{i};
                remred = true;
                chain.baseChange([base block], remred);

            end


        end

    end

end
