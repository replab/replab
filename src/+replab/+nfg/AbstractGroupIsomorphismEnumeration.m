classdef AbstractGroupIsomorphismEnumeration < replab.nfg.AbstractGroupIsomorphism
% Describes an isomorphism from an abstract group to its realization (permutation group) using brute force enumeration

    properties
        elements
        words
    end

    methods

        function self = AbstractGroupIsomorphismEnumeration(source, elements, words)
        % Constructs the nice isomorphism from an abstract finite group to its permutation group
        %
        % Args:
        %   source (`+replab.AbstractGroup`): Source of the isomorphism
        %   elements (`+replab.perm.Set`): Set of group elements
        %   words (cell(1,\*) of integer(1,\*)): Words corresponding to the group elements, described by their letters
            self.source = source;
            self.target = source.permutationGroup;
            self.elements = elements;
            self.words = words;
        end

        function iso1 = withUpdatedSource(self, source1)
            iso1 = replab.nfg.AbstractGroupIsomorphismEnumeration(source1, self.elements, self.words);
        end

    end

    methods (Static)

        function m = make(source)
            target = source.permutationGroup;
            n = target.domainSize;
            o = double(target.order);
            nG = target.nGenerators;
            gens = target.generators;
            invGens = cellfun(@(g) target.inverse(g), gens, 'uniform', 0);
            elements = replab.perm.Set(n);
            elements.insert(target.chain.allElements);
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
            m = replab.nfg.AbstractGroupIsomorphismEnumeration(source, elements, words);
        end

    end

    methods

        function T = imageGroup(self, S)
            T = S.niceGroup;
        end

        function t = imageElement(self, s)
            t = self.source.niceImage(s);
        end

        function s = preimageElement(self, t)
            s = self.source.fromLetters(self.words{self.elements.find(t')});
        end

        function S = preimageGroup(self, T)
            gens = cellfun(@(t) self.preimageElement(t), T.generators, 'uniform', 0);
            S = self.source.niceSubgroup(gens, T.order, T);
        end

    end

end
