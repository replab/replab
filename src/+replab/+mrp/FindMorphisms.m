classdef FindMorphisms
% Search for morphisms between groups up to conjugation
%
% Implements EPIMORPHISMS from
% D. Holt et al., Handbook of Computational Group Theory, Chapman & Hall/CRC, 2004, pp. 325-332


    properties (SetAccess = protected)
        r % (integer): Number of generators in `.F`
        F % (`+replab.AbstractGroup`): Abstract group to find homomorphisms from
        G % (`+replab.FiniteGroup`): Finite group to find homomorphisms to
        A % (`+replab.FiniteGroup`): Supergroup of `.G`
        I % (cell(1,r) of cell(1,\*) of elements of `.G`): Candidates
        CI % (cell(1,r) of cell(1,\*) of `+replab.FiniteGroup`): Candidates centralizers in `.A`
        relatorLetters % (cell(1,\*) of integer(1,\*)): Relators of F in letter form
        relatorSubsets % (integer(1,nFGens)): How many relators contain only the first ``k`` generators, see line 1 of EPIMORPHISMS
        filter % ({'morphisms', 'epimorphisms', 'isomorphisms'): Whether to look for morphisms, epimorphisms, isomorphisms
        single % (logical): Whether to return only the first result
    end

    methods

        function self = FindMorphisms(F, G, A, filter, single)
        % Constructor
            assert(G == A); % for now, contact D. Holt to know how to compute representative of conjugacy classes in subgroups
            r = F.nGenerators;
            h = cellfun(@(c) c.representative, G.conjugacyClasses.classes, 'uniform', 0);
            Ch = cellfun(@(c) c.representativeCentralizer, G.conjugacyClasses.classes, 'uniform', 0);
            heo = cellfun(@(g) G.elementOrder(g), h);
            I = cell(1, r);
            CI = cell(1, r);
            relatorLetters = cellfun(@(rel) F.factorizeLetters(rel), F.relators, 'uniform', 0);
            relmax = cellfun(@(rl) max(abs(rl)), relatorLetters);
            for i = 1:r
                relatorSubsets(i) = find([relmax > i, true], 1) - 1;
            end
            for i = 1:r
                t = 0;
                for j = 1:length(relatorLetters)
                    rel = relatorLetters{j};
                    if all(abs(rel) == i)
                        t = gcd(t, sum(sign(rel)));
                    end
                end
                Ii = cell(1, 0);
                CIi = cell(1, 0);
                for j = 1:length(h)
                    if isequal(filter, 'isomorphisms')
                        condition = (t == 0) || (t == heo(j));
                    else
                        condition = mod(t, heo(j)) == 0;
                    end
                    if condition
                        Ii{1,end+1} = h{j};
                        CIi{1,end+1} = Ch{j};
                    end
                end
                I{i} = Ii;
                CI{i} = CIi;
            end
            self.r = r;
            self.F = F;
            self.G = G;
            self.A = A;
            self.I = I;
            self.CI = CI;
            self.relatorLetters = relatorLetters;
            self.relatorSubsets = relatorSubsets;
            self.filter = filter;
            self.single = single;
        end

        function res = searchAll(self)
            assert(~self.single);
            red = self.searchUpToConjugation;
            res = cell(1, 0);
            for i = 1:length(red)
                f = red{i};
                G = f.image;
                innerElements = G/G.center;
                innerElements = innerElements.transversal;
                for i = 1:length(innerElements)
                    res{1,end+1} = f.andThen(G.conjugatingAutomorphism(innerElements{i}));
                end
            end
        end

        function res = searchUpToConjugation(self)
            I1 = self.I{1};
            CI1 = self.CI{1};
            res = cell(1, 0);
            % At first level, we do not need to conjugate the candidates, because we are looking for morphisms up to conjugation
            for i = 1:length(I1)
                g1 = I1{i};
                Cg1 = CI1{i};
                res1 = self.rec(2, {Cg1}, {g1}, {self.A.identity}, {g1});
                res = horzcat(res, res1);
                if self.single && ~isempty(res)
                    return
                end
            end
        end

        function res = testRelators(self, k, im)
        % Tests the relators for the given partial generator array
            res = false;
            for i = 1:self.relatorSubsets(k)
                if ~self.G.isIdentity(self.G.composeLetters(im, self.relatorLetters{i}))
                    return
                end
            end
            res = true;
        end

        function res = rec(self, k, Lambda, g, alpha, im)
        % Searchs for epimorphisms
        %
        % Note that our definition of conjugacy is such that conjugacy is a left action; thus we need to reverse some definitions.
        %
        % Args:
        %   k (integer, > 1): Level to investigate
        %   Lambda (cell(1,k-1) of `+replab.FiniteGroup`): Intersection of centralizers as per description in p. 329, Holt
        %   g (cell(1,k-1) of elements of `.G`): Partial list of generator candidates
        %   alpha (cell(1,k-1) of elements of `.A`): Conjugating elements
        %   im (cell(1,k-1) of elements of `.G`): Partial list of generators
            res = cell(1, 0);
            if k > self.r
                if ~isequal(self.filter, 'morphisms')
                    IM = self.G.subgroup(im);
                    if IM.order ~= self.G.order
                        return % empty res
                    end
                end
                res = {self.F.morphismByImages(self.G, 'preimages', self.F.generators, 'images', im, 'nChecks', 0)};
                return
            end
            Ik = self.I{k};
            CIk = self.CI{k};
            for i = 1:length(Ik)
                gk = Ik{i};
                Cgk = CIk{i};
                delta = self.A.doubleCosets(Lambda{k-1}, Cgk).transversal;
                for j = 1:length(delta)
                    deltaj = delta{j};
                    g{k} = gk;
                    alpha{k} = deltaj;
                    im{k} = self.A.leftConjugate(deltaj, gk);
                    if self.testRelators(k, im)
                        Lkprev = Lambda{k-1};
                        Lambda{k} = Lkprev.centralizer(im{k});
                        res1 = self.rec(k + 1, Lambda, g, alpha, im);
                        res = horzcat(res, res1);
                        if self.single && ~isempty(res)
                            return
                        end
                    end
                end
            end
        end

    end

end
