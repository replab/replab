classdef Standard < replab.Atlas
% A collection of descriptions of common groups, with presentations used in the mathematical literature
%
% Example:
%   >>> D4 = replab.DihedralGroup(4);
%   >>> D4.recognize.atlasGroup.name
%       'Dihedral group of order 8'
%   >>> S5 = replab.SymmetricGroup(5);
%   >>> S5.recognize.atlasGroup.name
%       'Symmetric group S(5) of degree 5'

    methods (Static)

        function A = instance
            persistent A_
            if isempty(A_)
                A_ = replab.atlas.Standard;
            end
            A = A_;
        end

    end

    methods

        function self = Standard
            self@replab.Atlas(1000);
        end

        function A = trivial(self)
        % Constructs the atlas entry corresponding to the trivial group
            name = 'Trivial group';
            prmGroup = replab.S(1); % this is a legit permutation representation
            generators = cell(1, 0);
            relators = cell(1, 0);
            % Presentation from the groupprops wiki
            % < x | x = 1 >
            A = replab.AbstractGroup(generators, prmGroup, relators, name);
        end

        function R = recognizeTrivial(self, G)
        % Recognizes if the given group is the trivial group and provides the generators according to the standard presentation
        %
        % The standard presentation is ``<x| x = id>``
            R = [];
            if ~G.isTrivial
                return
            end
            entry = self.trivial;
            R = replab.AtlasResult(G, entry, cell(1, 0));
        end

        function A = dihedral(self, n)
        % Constructs the atlas entry corresponding to the dihedral group of order 2*n
            assert(n > 2);
            name = sprintf('Dihedral group of order %d', 2*n);
            % Permutation realization
            X = [n:-1:1];
            A = [2:n 1];
            prmGroup = replab.PermutationGroup.of(X, A);
            % Presentation from the groupprops wiki
            % < x, a | a^n = x^2 = 1, x a x^-1 = a^-1 >
            relators = {['a^' num2str(n)] 'x^2' 'x a x^-1 a'};
            A = replab.AbstractGroup({'x' 'a'}, prmGroup, relators, name);
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

        function A = klein(self)
        % Constructs the atlas entry corresponding to the klein four-group
            name = sprintf('Klein four-group of order %d', 4);
            % Permutation realization
            X = [2,1,4,3];
            A = [3,4,1,2];
            prmGroup = replab.PermutationGroup.of(X, A);
            % Presentation from the groupprops wiki
            % < x, a | a^2 = x^2 = 1, x a x^-1 = a^-1 >
            relators = {'a^2' 'x^2' 'x a x^-1 a'};
            A = replab.AbstractGroup({'x' 'a'}, prmGroup, relators, name);
        end

        function R = recognizeKlein(self, G)
        % Recognizes if the given group is the Klein four-group and provides the generators according to the standard presentation
        %
        % The standard presentation is ``<x, a| a^2 = x^2 = id, x a x^-1 = a>``
            R = [];
            if G.order ~= 4
                return
            end
            if G.isCyclic
                return
            end
            x = G.generator(1);
            a = [];
            for i = 2:G.nGenerators
                g = G.generator(i);
                if ~G.eqv(x, g)
                    a = g;
                    break
                end
            end
            assert(~isempty(a));
            entry = self.klein;
            R = replab.AtlasResult(G, entry, {x a});
        end

        function A = symmetric(self, n)
        % Constructs the atlas entry corresponding to the symmetric group of degree n
            assert(n > 2);
            name = sprintf('Symmetric group S(%d) of degree %d', n, n);
            % Permutation realization
            S = [2:n 1];
            T = [2 1 3:n];
            prmGroup = replab.PermutationGroup.of(S, T);
            % this is the presentation from page 2100 of
            % https://www.ams.org/journals/tran/2003-355-05/S0002-9947-03-03040-X/S0002-9947-03-03040-X.pdf
            relators = {['s^' num2str(n)], 't^2', ['(s*t)^' num2str(n-1)]};
            for j = 2:floor(n/2)
                relators{1,end+1} = sprintf('(t^-1 s^-%d t s^%d)^2', j, j);
            end
            A = replab.AbstractGroup({'s' 't'}, prmGroup, relators, name);
            %outer = {replab.FiniteIsomorphism.identity(prmGroup)};
            %if n == 6
            %    imgS = [6 1 5 4 3 2];
            %    imgT = [2 1 4 3 6 5];
            %    outer{1,2} = prmGroup.morphismByImages(prmGroup, {imgS, imgT});
            %end
        end

        function R = recognizeSymmetric(self, G)
        % Recognizes if the given group is the symmetric group and provides the generators according to the standard presentationx
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
                                if entry.imagesDefineMorphism(G, {s t})
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

        function A = cyclic(self, n)
        % Constructs the cyclic group of order n
            assert(n >= 2);
            name = sprintf('Cyclic group C(%d) of order %d', n, n);
            % Permutation realization
            X = [2:n 1];
            prmGroup = replab.PermutationGroup.of(X);
            % standard presentation
            % < x | x^n = 1 >
            A = replab.AbstractGroup({'x'}, prmGroup, {['x^' num2str(n)]}, name);
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

        function A = alternating(self, n)
        % Constructs the alternating group of degree n
            assert(n >= 4);
            name = sprintf('Alternating group A(%d) of degree %d', n, n);
            isEven = mod(n, 2) == 0;
            % Permutation realization
            T = [2 3 1 4:n];
            if isEven
                S = [2 1 4:n 3];
            else
                S = [1 2 4:n 3];
            end
            prmGroup = replab.PermutationGroup.of(S, T);
            % this is the presentation from page 2100 of
            % https://www.ams.org/journals/tran/2003-355-05/S0002-9947-03-03040-X/S0002-9947-03-03040-X.pdf
            if isEven
                relators = {['s^' num2str(n-2)], 't^3', ['(s t)^' num2str(n-1)], '(t^-1 s^-1 t s)^2'};
            else
                relators = {['s^' num2str(n-2)], 't^3', ['(s t)^' num2str(n)]};
                for k = 1:floor((n-3)/2)
                    relators{1,end+1} = sprintf('(t s^-%d t s^%d)^2', k, k);
                end
            end
            A = replab.AbstractGroup({'s' 't'}, prmGroup, relators, name);
        end

        function R = recognizeAlternating(self, G)
        % Recognizes the alternating group and returns the generators corresponding to the standard presentation
            R = [];
            [n r] = replab.atlas.unfactorial(G.order*2);
            if r ~= 0
                return
            end
            n = double(n);
            C = G.conjugacyClasses;
            entry = self.alternating(n);
            for i = 1:length(C)
                S = C{i};
                s = S.representative;
                if G.elementOrder(s) == n-2
                    for j = 1:length(C)
                        T = C{j};
                        if G.elementOrder(T.representative) == 3
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

        function R = recognize(self, G)
            R = self.recognizeTrivial(G);
            if ~isempty(R)
                return
            end
            R = self.recognizeCyclic(G);
            if ~isempty(R)
                return
            end
            R = self.recognizeKlein(G);
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
            R = self.recognizeAlternating(G);
            if ~isempty(R)
                return
            end
            R = [];
        end

    end

end
