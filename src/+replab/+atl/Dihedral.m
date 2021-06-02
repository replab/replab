classdef Dihedral
% The family of dihedral groups

    methods (Static)

        function ct = characterTable(G, classes, stop, field)
        % Generates the real or complex character table for the dihedral group Dn
        %
        % From Fässler, A., Stiefel, E., & Wong, B. D. (1992). Group theoretical methods and their applications.
        % Boston: Birkhäuser, 23-24.
        %
        % Args:
        %   G (`+replab.PermutationGroup`): Permutation realization of the dihedral group
        %   classes (`+replab.ConjugacyClasses`): Conjugacy classes
        %   stop (integer): What is it
        %
        % Returns:
        %   ct (`+replab.CharacterTable`)
            nClasses = classes.nClasses;
            n = G.domainSize;
            ord = 2*n;

            % Irreps are generated first in 1D and then in 2D
            irreps = cell(1, nClasses);
            w = replab.cyclotomic.E(n);
            irreps{1} = G.repByImages('C', 1, 'images', {1, 1});
            irreps{2} = G.repByImages('C', 1, 'images', {1, -1});
            n1D = 2;
            if even(n)
                irreps{3} = G.repByImages('C', 1, 'images', {-1, 1});
                irreps{4} = G.repByImages('C', 1, 'images', {-1, -1});
                stop = stop - 1;
                n1D = 4;
            end
            for j = 1:stop
                g1 = replab.cyclotomic.zeros(2, 2);
                if field == 'R'
                    c = (w^j + w^(-j))/2;  % cos(2*pi*j/n)
                    s = (w^j - w^(-j))/2i; % sin(2*pi*j/n)
                    g1 = [c -s
                          s c];
                    g2 = replab.cyclotomic([1 0; 0 -1]);
                else
                    g1(1, 1) = w^j;
                    g1(2, 2) = w^(-j);
                    g2 = replab.cyclotomic([0 1; 1 0]);
                end
                irreps{n1D + j} = G.repByImages('C', 2, 'images', {g1, g2});
            end
            % Characters can be assigned to rotations then reflections
            chars = replab.cyclotomic.zeros(nClasses, nClasses);
            chars(1, :) = replab.cyclotomic(1);
            chars(1:n1D, 1) = replab.cyclotomic(1);
            if nClasses > n1D
                chars(n1D+1:nClasses, 1) = replab.cyclotomic(2);
            end
            if even(n)
                for k = 1:nClasses - 3
                    chars(2, k+1) = replab.cyclotomic(1);
                    chars(3:4, k+1) = replab.cyclotomic((-1)^k);
                    for j = 1:stop
                        pow = mod(j*k, ord);
                        chars(n1D+j, k+1) = w^pow + w^(-pow);
                    end
                end
                for k = 0:1
                    chars(2, nClasses - 1 + k) = replab.cyclotomic(-1);
                    chars(3, nClasses - 1 + k) = replab.cyclotomic((-1)^k);
                    chars(4, nClasses - 1 + k) = replab.cyclotomic((-1)^(k+1));
                    for j = 1:stop
                        chars(n1D+j, nClasses - 1 + k) = replab.cyclotomic(0);
                    end
                end
            else
                for k = 1:nClasses - 2
                    chars(2, k+1) = replab.cyclotomic(1);
                    for j = 1:stop
                        pow = mod(j*k, ord);
                        chars(n1D+j, k+1) = w^pow + w^(-pow);
                    end
                end
                chars(2, nClasses) = replab.cyclotomic(-1);
                for j = 1:stop
                    chars(n1D+j, nClasses) = replab.cyclotomic(0);
                end
            end
            if field == 'R'
                ct = replab.RealCharacterTable(G, classes, chars, 'irreps', irreps);
            else
                ct = replab.ComplexCharacterTable(G, classes, chars, 'irreps', irreps);
            end
        end

        function R = recognize(G)
        % Recognizes if the given group is the dihedral group and provides the generators according to the standard presentation
        %
        % The standard presentation is ``<a, x| a^n = x^2 = id, x a x^-1 = a>``
            assert(G.order > 2);
            R = [];
            x = [];
            a = [];
            if G.order <= 4 || G.order > 2^53-1
                % too big, or is the Klein four-group
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
            Rcyclic = replab.atl.Cyclic.recognize(Zn);
            a = Rcyclic.imageElement(Rcyclic.source.generator(1));
            entry = replab.atl.Dihedral.make(n);
            R = entry.isomorphismByImages(G, 'images', {a, x});
        end


        function [C, stop] = classReps(n)
        % Returns conjugacy classes representatives
        %
        % Classes will be listed with first rotations and then reflections
            if even(n)
                nClasses = 4 + n/2 - 1;
            else
                nClasses = 2 + (n - 1)/2;
            end
            C = cell(1, nClasses);
            d = [2:n, 1];
            s = fliplr(1:n);
            if even(n)
                stop = n/2;
            else
                stop = (n-1)/2;
            end
            rep = 1:n;
            for i = 1:stop + 1
                C{i} = rep;
                rep = rep(d);
            end
            C{i + 1} = s;
            if even(n)
                C{i + 2} = s(d);
            end
        end

        function G = make(n)
        % Constructs the atlas entry corresponding to the dihedral group of order 2*n
            assert(n > 2);
            name = sprintf('Dihedral group of order %d', 2*n);
            % Permutation realization
            A = [2:n 1];
            X = [n:-1:1];
            % Presentation from the groupprops wiki
            % < a, x | a^n = x^2 = 1, x a x^-1 = a^-1 >
            relators = {sprintf('a^%d', n) 'x^2' 'x a x^-1 a'};
            G = replab.AbstractGroup({'a', 'x'}, relators, 'permutationGenerators', {A, X}, 'order', vpi(2*n), 'inAtlas', true);
            [C, stop] = replab.atl.Dihedral.classReps(n);
            classList = cellfun(@(r) G.permutationGroup.conjugacyClass(r), C, 'uniform', 0);
            classes = replab.ConjugacyClasses(G.permutationGroup, classList);
            G.cache('conjugacyClasses', classes.imap(G.niceMorphism.inverse), 'error');
            ctR = replab.atl.Dihedral.characterTable(G.permutationGroup, classes, stop, 'R');
            ctC = replab.atl.Dihedral.characterTable(G.permutationGroup, classes, stop, 'C');
            G.cache('realCharacterTable', ctR.imap(G.niceMorphism.inverse), 'error');
            G.cache('complexCharacterTable', ctC.imap(G.niceMorphism.inverse), 'error');
        end

    end

end