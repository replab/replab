classdef Cyclic

    methods (Static)

        function R = recognize(G)
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
            entry = replab.atl.Cyclic.make(n);
            R = entry.isomorphismByImages(G, 'images', {x});
        end

        function G = make(n)
        % Constructs the cyclic group of order n and its character table
        %
        % From Fässler, A., Stiefel, E., & Wong, B. D. (1992). Group theoretical methods and their applications.
        % Boston: Birkhäuser, 22.
        %
        % Parts written by Joscelyn
            assert(n >= 2);
            name = sprintf('Cyclic group C(%d) of order %d', n, n);
            % Permutation realization
            X = [2:n 1];
            % standard presentation
            % < x | x^n = 1 >
            G = replab.AbstractGroup({'x'}, {['x^' num2str(n)]}, 'permutationGenerators', {X}, 'order', vpi(n), 'name', name, 'inAtlas', true);
            classArray = arrayfun(@(r) G.conjugacyClass(sprintf('x^%d', r)), 0:n-1, 'uniform', 0);
            classes = replab.ConjugacyClasses(G, classArray);

            % Generate irreps with images as increasing powers of E(n)
            w = replab.cyclotomic.E(n);
            irreps = arrayfun(@(x) G.repByImages('C', 1, 'images', {w^x}), 0:n-1, 'uniform', 0);

            % Generate characters with increasing powers of conjugacy classes and irreps
            chars = replab.cyclotomic.zeros(n, n);
            for i = 1:n
                for j = 1:n
                    chars(i, j) = w^mod((i-1) * (j-1), n);
                end
            end
            ct = replab.ComplexCharacterTable(G, classes, chars, 'irreps', irreps);
            G.cache('complexCharacterTable', ct, 'error');
            G.cache('conjugacyClasses', classes, 'error');
        end

    end

end