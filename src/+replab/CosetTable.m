classdef CosetTable < replab.Str
% Describes a coset table
%
% Note that the action of words on cosets is a left action, contrary to the algorithms described in Holt.
% Thus, we implement the algorithms of Holt in a different class, which matches the pseudocode closely.
% This class, however, translates between those conventions.

    properties (Access = protected)
        generatorNames % (cell(1,\*) of charstring): Names of the generators
        internal % (`+replab.+fp.CosetTable`): Internal coset table
    end

    methods (Static)


    end

    methods

        function self = CosetTable(generatorNames, internal)
            self.generatorNames = generatorNames;
            self.internal = internal;
        end

        function t = table(self)
        % Returns the coset table as a prettyprintable table
        %
        % Returns:
        %   `+replab.+str.Table`: String table
            nG = length(self.generatorNames);
            C = self.internal.C;
            C = C(:,[nG+1:2*nG 1:nG]); % swap generators and generator inverses, so that we show a left action
            t = replab.str.Table(C);
            invNames = cellfun(@(x) ['inv(', x, ')'], self.generatorNames, 'uniform', 0);
            t.addColumnNames([self.generatorNames, invNames]);
            t.addRowNames(num2cell(1:size(self.internal.C, 1)));
            t.setRowSep(1, '-');
            t.setColSep(1, ' | ');
        end

    end

    methods (Static)

        function ct = fromPresentation(generatorNames, relatorWords, subgroupGeneratorWords)
        % Enumerate cosets
        %
        % Args:
        %   generatorNames (cell(1,\*) of charstring): Group generator names
        %   relatorWords (cell(1,\*) of charstring): Group relators to parse
        %   subgroupGeneratorWords (cell(1, \*) of charstring): Generators for a subgroup, expressed as words in the group generators
            nGenerators = length(generatorNames);
            subgroupGenerators = cellfun(@(w) replab.fp.parseLetters(w, generatorNames), subgroupGeneratorWords, 'uniform', 0);
            relators = cellfun(@(w) replab.fp.parseLetters(w, generatorNames), relatorWords, 'uniform', 0);
            internal = replab.fp.CosetTable(nGenerators, relators, subgroupGenerators, 2^51 - 1); % TODO: estimate memory use
            internal.cosetEnumerationR;
            ct = replab.CosetTable(generatorNames, internal);
        end

    end

end
