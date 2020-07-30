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

        function s = headerStr(self)
            s = ['AtlasResult(' self.entry.name ')'];
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Str(self);
            names{1,end+1} = 'presentationString';
            values{1,end+1} = self.presentationString;
        end

        function f = presentationString(self)
            f = self.entry.abstractGroup.presentationString;
        end

        function G = allStandardGenerators(self)
        % Returns all variants of standard generators as rows in a cell array
            abGroup = self.entry.abstractGroup;
            prmGroup = abGroup.permutationGroup;
            center = prmGroup.center;
            assert(~isempty(center));
            assert(~isempty(self.entry.outerRepresentatives));
            cosets = center \ prmGroup;
            T = cosets.transversal;
            O = self.entry.outerRepresentatives;
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
            % and Out(G) = Inn(G) \ Auto(G)
            abGens = abGroup.generators;
            prmImages = prmGroup.generators;
            grpImages = cellfun(@(g) self.niceMorphism.imageElement(g), abGens, 'uniform', 0);
            h = prmGroup.morphismByImages(self.group, grpImages);
            G = cell(length(T)*length(O), length(abGens));
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
