classdef FactorizationEnumeration < replab.mrp.Factorization
% Factorizes permutations using breath-first orbit computation

    properties (SetAccess = protected)
        elements % (`+replab.perm.Set`): Set of group elements
        words % (cell(1,\*) of integer(1,\*)): Words corresponding to the group elements, described by their letters
    end

    methods

        function self = FactorizationEnumeration(group, elements, words)
            self.group = group;
            self.elements = elements;
            self.words = words;
        end

    end

    methods % Implementations

        function letters = factorize(self, g)
            ind = self.elements.find(g');
            assert(~isempty(ind), 'The permutation %s is not a member of the group.', replab.shortStr(g));
            letters = self.words{ind};
        end

    end

    methods (Static)

        function m = make(group, generators, useInverses)
        % Constructs a group factorization object by enumerating all elements of a permutation group
        %
        % Args:
        %   group (`+replab.PermutationGroup`): Group to decompose elements of
        %   generators (cell(1,\*) of elements of ``group``): Group generators (default: ``group.generators``)
        %   useInverses (logical, optional): Whether to use inverses in the decomposition (default: true)
        %
        % Returns:
        %   `+replab.+mrp.Factorization`: The factorization object that can compute preimages in words of the generators
            if nargin < 2 || isempty(generators)
                generators = group.generators;
            end
            if nargin < 3 || isempty(useInverses)
                useInverses = true;
            end
            n = group.domainSize;
            o = double(group.order);
            nG = length(generators);
            gens = generators;
            if useInverses
                invGens = cellfun(@(g) group.inverse(g), gens, 'uniform', 0);
            end
            elements = replab.perm.Set(n);
            elements.insert(group.chain.allElements);
            ind = elements.find((1:n)');
            ok = false(1, o);
            ok(ind) = true;
            words = cell(1, o);
            words{ind} = [];
            computed = ind;
            l = 1;
            while l <= length(computed)
                p = elements.at(computed(l))';
                w = words{computed(l)};
                if ~isempty(w)
                    i = abs(w(1));
                    range = [i 1:i-1 i+1:nG];
                else
                    range = 1:nG;
                end
                for i = range
                    g = gens{i};
                    p1 = g(p);
                    w1 = [i w];
                    ind = elements.find(p1');
                    if ~ok(ind)
                        ok(ind) = true;
                        words{ind} = w1;
                        computed = [computed ind];
                    end
                end
                if useInverses
                    for i = range
                        g = invGens{i};
                        p1 = g(p);
                        w1 = [-i w];
                        ind = elements.find(p1');
                        if ~ok(ind)
                        ok(ind) = true;
                        words{ind} = w1;
                        computed = [computed ind];
                        end
                    end
                end
                l = l + 1;
            end
            m = replab.mrp.FactorizationEnumeration(group, elements, words);
        end

    end

end
