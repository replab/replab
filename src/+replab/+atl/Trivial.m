classdef Trivial

    methods (Static)

        function G = make
            name = 'Trivial group';
            generators = cell(1, 0);
            relators = cell(1, 0);
            A = replab.AbstractGroup(cell(1, 0), cell(1, 0), 'permutationGenerators', cell(1, 0), 'order', vpi(1), 'name', 'TrivialGroup', 'inAtlas', true);
            classes = replab.ConjugacyClasses(A, {A.conjugacyClass(A.identity)});
            irrepsC = {A.trivialRep('C', 1)};
            irrepsR = {A.trivialRep('R', 1)};
            values = replab.cyclotomic(1);
            A.cache('conjugacyClasses', classes, 'error');
            A.cache('realCharacterTable', replab.RealCharacterTable(A, classes, values, 'irreps', irrepsR);
            A.cache('complexCharacterTable', replab.RealCharacterTable(A, classes, values, 'irreps', irrepsC);
        end

        function R = recognize(G)
        % Recognizes if the given group is the symmetric group and provides the generators according to the standard presentation
            if G.isTrivial
                entry = replab.atl.Trivial.make;
                R = entry.isomorphismByImages(G);
            else
                R = [];
            end
        end

    end

end
