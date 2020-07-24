classdef CharacterTable < replab.Obj
% Describes the character table of a group
%
% Example:
%   >>> s3ct = replab.CharacterTable.forPermutationGroup(replab.S(3));
%   >>> disp(s3ct.table)
%            [1, 2, 3]  [1, 3, 2]  [2, 3, 1]
%       X.1      1          1          1
%       X.2      2          0         -1
%       X.3      1         -1          1

    properties (SetAccess = protected)
        group % (`+replab.Group`): group represented by character table
        irreps % (cell(1, nclasses) of `+replab.Rep`): irreducible representations
        classes % (cell(1, nclasses) of `+replab.ConjugacyClass`): conjugacy classes
        chars % (double(nclasses, nclasses)): matrix of characters
        table % (`+replab.str.Table`): table object for display purposes
    end

    methods

        function self = CharacterTable(group, irreps, classes, chars)
        %  Assigns properties and generates Table with given properties
            self.group = group;
            self.irreps = irreps;
            self.classes = classes;
            self.chars = chars;

            % make Table object of character table
            colnames = cellfun(@(v) v.representative, self.classes, 'UniformOutput', false);
            nirreps = length(self.irreps);
            rownames = cellfun(@(n)['X.', num2str(n)], num2cell(1:nirreps), 'UniformOutput', false);
            self.table = replab.str.Table(self.chars);
            self.table.addRowNames(rownames);
            self.table.addColumnNames(colnames);

        end

        function s = headerStr(self)
            group_str = replab.headerStr(self.group);
            if length(group_str) > 1
                group_str = [lower(group_str(1)), group_str(2:end)];
            end
            s = sprintf(['Character table for ', group_str]);
        end

        function useBorders(self, logical)
        %  Turn on and off borders on the table
        %
        % Args:
        %   logical (logical): whether borders are on (true) or off (false)
            if logical
                self.table.setColSep(0, '| ')
                self.table.setColSep(self.table.nColumns, ' |')
                self.table.setColSep(1:self.table.nColumns-1, ' | ')
                self.table.setRowSep(0:self.table.nRows, '-')
            else
                self.table.setColSep(0:self.table.nColumns, '  ')
                self.table.setRowSep(0:self.table.nRows, '')
            end
        end

        function setIrrepLabels(self, labels)
        % Set the labels for the irreducible representation in the table
        %
        % Args:
        %   labels ({cell(1,nirreps), 'default'): cell array of character string labels
        %                                         or 'default' to use default labels
            if iscell(labels)
                if length(labels) == length(self.table.getRowNames)
                    self.table.addRowNames(labels)
                end
            elseif isequal(labels, 'default')
                nirreps = length(self.irreps);
                rownames = cellfun(@(n)['X.', num2str(n)], num2cell(1:nirreps), 'UniformOutput', false);
                self.table.addRowNames(rownames)
            end
        end

        function setClassLabels(self, labels)
        % Set the labels for the conjugacy classes in the table
        %
        % Args:
        %   labels (cell(1,nclasses)): cell array of character string labels
            if iscell(labels)
                if length(labels) == length(self.table.getColumnNames)
                    self.table.addColumnNames(labels)
                end
            elseif isequal(labels, 'default')
                colnames = cellfun(@(v) v.representative, self.classes, 'UniformOutput', false);
                self.table.addColumnNames(colnames)
            end
        end

        function table = pointGroupTable(self)
        % Returns a table with labels in chemistry notation
            table = self.table;
            irrepLabels = table.getRowNames;
            for i = 1:length(irrepLabels)
                sym = self.mulliken(self.irreps{i});
                irrepLabels{i} = sym;
            end
            table.addRowNames(irrepLabels)
        end

        function sizes = classSizes(self)
            sizes = cellfun(@(c) c.size, self.classes);
        end

    end

    methods (Static)

        function sym = mulliken(rep)
        % Returns the Mulliken symbol of a representation
        %
        % Requires knowledge of the principal axis to fully determine
            if rep.dimension == 1
                sym = 'A/B';
            else
                sym = 'E';
            end
        end

        function table = forPermutationGroup(group)
        % Generates character table for a permutation group
        %
        % This currently assumes that the characters are integers, and the algorithm is slow.
        %
        % Args:
        %    group (`+replab.PermutationGroup`): Permutation group
        %
        % Returns:
        %    table (`+replab.CharacterTable`): Character table of group
            assert(isa(group, 'replab.PermutationGroup'));
            ord = double(group.order);
            decomp = group.naturalRep.decomposition.nice;
            irreps = decomp.components;
            classes = group.conjugacyClasses;
            k = length(classes);
            for i = 1:k
                % verify that g^n is conjugate to g for all g \in G and n \in Z with n and (g) coprime
                % https://math.stackexchange.com/questions/2792741/classification-of-groups-with-integer-valued-characters
                cl = classes{i};
                g = cl.representative;
                eo = group.elementOrder(g);
                for j = 2:eo-1
                    if gcd(j, eo) == 1
                        if ~cl.contains(group.composeN(g, j))
                            error('Group does not have integer characters');
                        end
                    end
                end
            end
            ccreps = cell(1, length(classes));
            cclens = cell(1, length(classes));
            for i = 1:k
                ccreps{i} = classes{i}.representative;
                cclens{i} = double(classes{i}.size);
            end
            nirreps = length(irreps);
            chars = cell(nirreps, k);
            for i = 1:nirreps
                irrep = irreps{i};
                chars(i,:) = cellfun(@(x) trace(irrep.image(x)), ccreps, 'UniformOutput', false);
            end
            chars = cell2mat(chars);
            start_irreps = 1;
            while nirreps ~= k
                for i = start_irreps:nirreps
                    for j = start_irreps:nirreps
                        new_rep = irreps{i}.kron(irreps{j});
                        ss = 0;
                        for n = 1:k
                            char = trace(new_rep.image(ccreps{n}));
                            ss = ss + cclens{n}*(char)^2;
                        end
                        if ss ~= ord
                            new_irreps = new_rep.decomposition.nice;
                            start_irreps = nirreps + 1;
                            for m = 1:new_irreps.nComponents
                                new_irrep = new_irreps.component(m);
                                new_char = cellfun(@(x) trace(new_irrep.image(x)), ...
                                                    ccreps, 'UniformOutput', false);
                                new_char = cell2mat(new_char);
                                if ~any(all(abs(chars - new_char) < 1/2, 2))
                                    irreps{nirreps + 1} = new_irrep;
                                    chars(nirreps + 1, :) = new_char;
                                    nirreps = nirreps + 1;
                                end
                            end
                        end
                    end
                end
            end
            chars = round(chars);
            table = replab.CharacterTable(group, irreps, classes, chars);
        end

    end
end