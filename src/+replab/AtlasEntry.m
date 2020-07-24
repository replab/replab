classdef AtlasEntry < replab.Str
% Describes the properties of a finite group
%
% It contains a human readable name of the group, a presentation and a realization as a permutation group.
%
% It can optionally contain information to enumerate all automorphism of the group, as to find all possible
% sets of standard generators when recognizing a user-defined group.

    properties
        name % (charstring): Group name
        abstractGroup % (`+replab.AbstractGroup`): Group presentation
        prmGroup % (`+replab.PermutationGroup`): Realization as permutation group

        outerRepresentatives % (cell(1,\*) of `+replab.Morphism`): Right coset representatives of the outer automorphism group
                             %
                             %                                     Given as morphisms prmGroup -> prmGroup
    end

    methods

        function self = AtlasEntry(name, abstractGroup, prmGroup, outerRepresentatives)
            self.name = name;
            self.abstractGroup = abstractGroup;
            self.prmGroup = prmGroup;
            if nargin == 4
                self.outerRepresentatives = outerRepresentatives;
            else
                self.outerRepresentatives = [];
            end
        end

    end

end
