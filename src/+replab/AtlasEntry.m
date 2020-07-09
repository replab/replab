classdef AtlasEntry < replab.Str
% Describes the properties of a finite group
%
% It contains a human readable name of the group, a presentation and a realization as a permutation group.
%
% It can optionally contain information to enumerate all automorphism of the group, as to find all possible
% sets of standard generators when recognizing a user-defined group.

    properties
        atlas % (`+replab.Atlas`): Atlas in which this group is present
        name % (charstring): Group name
        fpGroup % (`+replab.FiniteFPGroup`): Group presentation
        prmGroup % (`+replab.PermutationGroup`): Realization as permutation group

        outerRepresentatives % (cell(1,\*) of `+replab.Morphism`): Right coset representatives of the outer automorphism group
                             %
                             %                                     Given as morphisms prmGroup -> prmGroup
        % TODO: conjugacy classes
        % TODO: character table
    end

    methods

        function self = AtlasEntry(atlas, name, fpGroup, prmGroup, outerRepresentatives)
            self.atlas = atlas;
            self.name = name;
            self.fpGroup = fpGroup;
            self.prmGroup = prmGroup;
            if nargin == 5
                self.outerRepresentatives = outerRepresentatives;
            else
                self.outerRepresentatives = [];
            end
        end

    end

end
