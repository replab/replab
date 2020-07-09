classdef AtlasResult < replab.Str
% Identifies a user-defined group as a standard group present in an atlas

    properties
        group % (`+replab.NiceFiniteGroup`): User-defined group
        entry % (`+replab.AtlasEntry`): Entry of the group in an atlas
        standardGenerators % (cell(1,\*) of elements of group): Generators for `group` that respect the standard presentation in `entry`
    end

    methods

        function self = AtlasResult(group, entry, standardGenerators)
            self.group = group;
            self.entry = entry;
            self.standardGenerators = standardGenerators;
        end

        function G = allStandardGenerators(self)
        % Returns all variants of standard generators as rows in a cell array
            prmGroup = self.entry.prmGroup;
            fpGroup = self.entry.fpGroup;
            center = prmGroup.center;
            assert(~isempty(center));
            assert(~isempty(self.outerAutoRepresentatives));
            cosets = center \ prmGroup;
            T = cosets.transversal;
            O = self.entry.outerAutoRepresentatives;
            A = {};
            % We consider the group of all automorphisms of G, Aut(G)
            % This group contains the normal subgroup of inner automorphisms of G, Inn(G)
            % In turn, Inn(G) is G / Center(G)
            %
            % Thus, with T the transversal of the cosets G / Center(G) in G given
            % as group elements acting by left conjugation,
            %
            % and O a transversal of the cosets Aut(G) / Inn(G) in Aut(G) given
            % as morphisms
            %
            % we can enumerate all automorphism by the composition of g -> o(t g t^-1)
            %
            %
            % and Out(G) = Inn(G) \ Auto(G)
            fpGens = fpGroup.generators;
            prmImages = prmGroup.generators;
            grpImages = cellfun(@(g) self.isomorphism(g), fpGens, 'uniform', 0);
            h = prmGroup.morphismByImages(self.group, grpImages);
            G = cell(length(T)*length(O), length(fpGens));
            ind = 1;
            for i = 1:length(T)
                t = prmGroup.leftConjugateMorphism(T{i});
                for j = 1:length(O)
                    o = O{j};
                    G(ind,:) = cellfun(@(g) h(o(t(g))), prmGroup.generators, 'uniform', 0);
                end
            end
        end

    end

end
