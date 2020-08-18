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

        function letters = preimageElement(self, g)
            ind = self.elements.find(g');
            assert(~isempty(ind), 'The permutation %s is not a member of the group.', replab.shortStr(g));
            letters = self.words{ind};
        end

    end

    methods (Static)

        function m = make(group)
            n = group.domainSize;
            o = double(group.order);
            nG = group.nGenerators;
            gens = group.generators;
            invGens = cellfun(@(g) group.inverse(g), gens, 'uniform', 0);
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
                l = l + 1;
            end
            m = replab.mrp.FactorizationEnumeration(group, elements, words);
        end

    end

end
