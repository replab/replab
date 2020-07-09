classdef Standard < replab.Atlas

    methods

        function self = Standard
            self@replab.Atlas(1000);
        end

        function E = dihedral(self, n)
            assert(n > 2);
            name = sprintf('Dihedral group of order %d', 2*n);
            % Permutation realization
            X = [n:-1:1];
            A = [2:n 1];
            prmGroup = replab.PermutationGroup.of(X, A);
            % Presentation from the groupprops wiki
            % < x, a | a^n = x^2 = 1, x a x^-1 = a^-1 >
            [F x a] = replab.FreeGroup.of('x', 'a');
            relators = {a^n, x^2, x*a/x*a};
            fpGroup = F / relators;
            fpGroup.setPermutationImages(prmGroup.generators);
            E = replab.AtlasEntry(self, name, fpGroup, prmGroup);
        end

        function R = recognizeDihedral(self, G)
        % Recognizes if the given group is the dihedral group and provides the generators according to the standard presentation
        %
        % The standard presentation is ``<x, a| a^n = x^2 = id, x a x^-1 = a>``
            assert(G.order > 2);
            R = [];
            x = [];
            a = [];
            if G.order > 2^53-1
                return
            end
            if mod(G.order, 2) == 1
                return
            end
            n = double(G.order/2);
            if mod(n, 2) == 1
                Zn = G.derivedSubgroup;
                if ~(Zn.isCyclic && Zn.order == n)
                    return
                end
            else
                G1 = G.derivedSubgroup;
                if ~(G1.isCyclic && G1.order == n/2)
                    return
                end
                T = G.rightCosetsOf(G1).transversal;
                good = false;
                for i = 1:4
                    Zn = G1.closure(T{i});
                    if Zn.isCyclic && Zn.order == n
                        good = true;
                        break
                    end
                end
                if ~good
                    return
                end
            end
            x = G.sample;
            while Zn.contains(x)
                x = G.sample;
            end
            if ~(G.elementOrder(x) == 2 && all(cellfun(@(a) G.isIdentity(G.composeAll({x a x a})), Zn.generators)))
                return
            end
            Rcyclic = self.recognizeCyclic(Zn);
            a = Rcyclic.standardGenerators{1};
            entry = self.dihedral(n);
            R = replab.AtlasResult(G, entry, {x a});
        end

        function E = symmetric(self, n)
            assert(n > 2);
            name = sprintf('Symmetric group S(%d) of degree %d', n, n);
            % Permutation realization
            S = [2:n 1];
            T = [2 1 3:n];
            prmGroup = replab.PermutationGroup.of(S, T);
            % this is the presentation from page 2100 of
            % https://www.ams.org/journals/tran/2003-355-05/S0002-9947-03-03040-X/S0002-9947-03-03040-X.pdf
            [F s t] = replab.FreeGroup.of('s', 't');
            relators = {s^n, t^2, (s*t)^(n-1)};
            for j = 2:floor(n/2)
                comm = inv(t)*inv(s^j)*t*s^j;
                relators{1,end+1} = comm^2;
            end
            fpGroup = F / relators;
            fpGroup.setPermutationImages(prmGroup.generators);
            outer = {replab.Morphism.identity(prmGroup)};
            if n == 6
                imgS = [6 1 5 4 3 2];
                imgT = [2 1 4 3 6 5];
                outer{1,2} = prmGroup.morphismByImages(prmGroup, {imgS, imgT});
            end
            E = replab.AtlasEntry(self, name, fpGroup, prmGroup, outer);
        end

        function R = recognizeSymmetric(self, G)
            R = [];
            [n r] = replab.atlas.unfactorial(G.order);
            if r ~= 0
                return
            end
            n = double(n);
            C = G.conjugacyClasses;
            entry = self.symmetric(n);
            for i = 1:length(C)
                S = C{i};
                s = S.representative;
                if G.elementOrder(s) == n
                    for j = 1:length(C)
                        T = C{j};
                        if G.elementOrder(T.representative) == 2
                            U = T.elements;
                            for k = 1:length(U)
                                t = U{k};
                                if entry.fpGroup.imagesDefineMorphism(G, {s t})
                                    if G.subgroup({s, t}).order == G.order
                                        R = replab.AtlasResult(G, entry, {s t});
                                        return
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        function E = cyclic(self, n)
            assert(n >= 2);
            name = sprintf('Cyclic group C(%d) of order %d', n, n);
            % Permutation realization
            X = [2:n 1];
            prmGroup = replab.PermutationGroup.of(X);
            % standard presentation
            [F x] = replab.FreeGroup.of('x');
            % < x | x^n = 1 >
            relators = {x^n};
            fpGroup = F / relators;
            fpGroup.setPermutationImages(prmGroup.generators);
            E = replab.AtlasEntry(self, name, fpGroup, prmGroup);
        end

        function R = recognizeCyclic(self, G)
        % Recognizes if the given group is the cyclic group and provides the group generator
        %
        % The standard presentation is ``<x| x^n = id>``
            R = [];
            if G.isTrivial
                return
            end
            if G.order > 2^53-1
                return
            end
            n = double(G.order);
            x = [];
            if ~G.isCyclic
                return
            end
            x = G.sample;
            while G.elementOrder(x) ~= G.order
                x = G.sample;
            end
            entry = self.cyclic(n);
            R = replab.AtlasResult(G, entry, {x});
        end

        function E = alternating(G, n)
            assert(n > 2);
            name = sprintf('Alternating group A(%d) of degree %d', n, n);
            isEven = mod(n, 2) == 0;
            % Permutation realization
            T = [2 3 1 4:n];
            if isEven
                S = [2 1 4:n 3];
            else
                sS= [1 2 4:n 3];
            end
            prmGroup = replab.PermutationGroup.of(S, T);
            % this is the presentation from page 2100 of
            % https://www.ams.org/journals/tran/2003-355-05/S0002-9947-03-03040-X/S0002-9947-03-03040-X.pdf
            [F s t] = replab.FreeGroup.of('s', 't');
            if isEven
                relators = {s^(n-2), t^3, (s*t)^(n-1), (inv(t)*inv(s)*t*s)^2};
            else
                relators = {s^(n-2), t^3, (s*t)^n};
                for k = 1:floor((n-3)/2)
                    relators{1,end+1} = (t*s^(-k)*t*s^k)^2;
                end
            end
            fpGroup = F / relators;
            fpGroup.setPermutationImages(prmGroup.generators);
            E = replab.AtlasEntry(self, name, fpGroup, prmGroup);
        end

        function R = recognize(self, G)
            R = self.recognizeCyclic(G);
            if ~isempty(R)
                return
            end
            R = self.recognizeDihedral(G);
            if ~isempty(R)
                return
            end
            R = self.recognizeSymmetric(G);
            if ~isempty(R)
                return
            end
            R = [];
        end

    end

end
