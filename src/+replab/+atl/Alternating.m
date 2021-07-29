classdef Alternating

    methods (Static)

        function R = recognize(G)
        % Recognizes the alternating group and returns the generators corresponding to the standard presentation
            R = [];
            [n r] = replab.util.unfactorial(G.order*2);
            if r ~= 0
                return
            end
            n = double(n);
            C = G.conjugacyClasses.classes;
            entry = replab.atl.Alternating.make(n);
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
                                if entry.isMorphismByImages(G, 'images', {s t})
                                    if G.subgroup({s, t}).order == G.order
                                        R = entry.isomorphismByImages(G, 'images', {s t});
                                        return
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        function G = make(n)
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
            G = replab.AbstractGroup({'s' 't'}, relators, 'permutationGenerators', {S, T}, 'order', replab.util.factorial(n)/2, 'name' , name, 'inAtlas', true);
        end

    end

end
