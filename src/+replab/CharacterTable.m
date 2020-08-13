classdef CharacterTable < replab.Obj
% Describes the character table of a group
%
% $$$ % Example:
% $$$ %   >>> s3ct = replab.CharacterTable.forPermutationGroup(replab.S(3));
% $$$ %   >>> disp(s3ct.table)
% $$$ %            [1, 2, 3]  [1, 3, 2]  [2, 3, 1]
% $$$ %       X.1      1          1          1
% $$$ %       X.2      2          0         -1
% $$$ %       X.3      1         -1          1
%
% The character values are stored as elements of the cyclotomic field, using the `.cyclotomic` class which requires
% external libraries and a Java Virtual Machine available.
%
% Instances of `.CharacterTable` are immutable; however, they can be enriched by calling the appropriate methods to
% add contextual information.

    properties (SetAccess = protected)
        group % (`+replab.FiniteGroup`): Group represented by character table
        classes % (cell(1, nClasses) of `+replab.ConjugacyClass`): Conjugacy classes of `.group`
        irreps % (cell(1, nClasses) of ``[]`` or `.RepByImages`): Explicit matrix representations (can contain empty values)
        characters % (`.cyclotomic` (nClasses, nClasses)): Character values
    end

% $$$     methods (Static)
% $$$
% $$$         function ct = make(group, conjugacyClasses, conjugacyClassNames, characterExpressions, irrepNames, irrepExpressions)
% $$$             if isempty(conjugacyClassNames)
% $$$                 conjugacyClassNames = cellfun(@(c) replab.shortStr(c.representative), conjugacyClasses, 'uniform', 0);
% $$$             end
% $$$             characterValues = cellfun(@(str) replab.cyclo.Parser.parse(str), characterExpressions);
% $$$             if isempty(irrepNames)
% $$$                 ind = 1;
% $$$                 for i = 1:size(characterValues, 1)
% $$$                     if all(characterValues(i, :) == 1)
% $$$                         irrepNames{i} = 't';
% $$$                     else
% $$$                         irrepNames{i} = sprintf('r_%d', ind);
% $$$                         ind = ind + 1;
% $$$                     end
% $$$                 end
% $$$             end
% $$$             ct = replab.CharacterTable(group, conjugacyClasses, conjugacyClassNames, irrepNames, characterExpressions, characterValues, irrepExpressions);
% $$$         end
% $$$
% $$$     end

    methods

        function self = CharacterTable(group, classes, characters)
        % Constructs a character table
            self.group = group;
            self.classes = classes;
            self.characters = characters;
        end

        function s = headerStr(self)
            group_str = replab.headerStr(self.group);
            if length(group_str) > 1
                group_str = [lower(group_str(1)), group_str(2:end)];
            end
            s = sprintf(['Character table for ', group_str]);
        end

        function v = findInTable(self, c, hshs, vars, vals)
            h = c.hash;
            inds = find(hshs == h);
            v = [];
            for i = inds
                if vals(i) == c
                    v = vars{i};
                    return
                end
            end
        end

        function [cells vars vals] = tableValues(self)
            chars = self.characters;
            cells = cell(size(chars));
            hshs = zeros(1, 0);
            vars = cell(1, 0);
            vals = replab.cyclotomic.zeros(1, 0);
            for i = 1:size(chars, 1)
                for j = 1:size(chars, 2)
                    c = chars(i, j);
                    d = double(c);
                    if c == 1i
                        cells{i,j} = 'i ';
                    elseif c == -1i
                        cells{i,j} = '-i ';
                    elseif isreal(d) && floor(d) == d
                        cells{i,j} = [strtrim(num2str(c)) ' '];
                    else
                        v = self.findInTable(c, hshs, vars, vals);
                        if ~isempty(v)
                            cells{i,j} = [v ' '];
                            continue
                        end
                        v = self.findInTable(conj(c), hshs, vars, vals);
                        if ~isempty(v)
                            cells{i,j} = [v '*'];
                            continue
                        end
                        v = self.findInTable(-c, hshs, vars, vals);
                        if ~isempty(v)
                            cells{i,j} = ['-' v ' '];
                            continue
                        end
                        v = self.findInTable(-conj(c), hshs, vars, vals);
                        if ~isempty(v)
                            cells{i,j} = ['-' v '*'];
                            continue
                        end
                        ind = length(vars) + 1;
                        hshs(1,ind) = c.hash;
                        vars{1,ind} = char('A'+ind-1);
                        vals = [vals c];
                        cells{i,j} = vars{1,ind};
                    end
                end
            end
        end

    end
% $$$
% $$$
% $$$         function useBorders(self, logical)
% $$$         %  Turn on and off borders on the table
% $$$         %
% $$$         % Args:
% $$$         %   logical (logical): whether borders are on (true) or off (false)
% $$$             if logical
% $$$                 self.table.setColSep(0, '| ')
% $$$                 self.table.setColSep(self.table.nColumns, ' |')
% $$$                 self.table.setColSep(1:self.table.nColumns-1, ' | ')
% $$$                 self.table.setRowSep(0:self.table.nRows, '-')
% $$$             else
% $$$                 self.table.setColSep(0:self.table.nColumns, '  ')
% $$$                 self.table.setRowSep(0:self.table.nRows, '')
% $$$             end
% $$$         end
% $$$
% $$$         function setIrrepLabels(self, labels)
% $$$         % Set the labels for the irreducible representation in the table
% $$$         %
% $$$         % Args:
% $$$         %   labels ({cell(1,nirreps), 'default'): cell array of character string labels
% $$$         %                                         or ``'default'`` to use default labels
% $$$             if iscell(labels)
% $$$                 if length(labels) == length(self.table.getRowNames)
% $$$                     self.table.addRowNames(labels)
% $$$                 end
% $$$             elseif isequal(labels, 'default')
% $$$                 nirreps = length(self.irreps);
% $$$                 rownames = cellfun(@(n)['X.', num2str(n)], num2cell(1:nirreps), 'UniformOutput', false);
% $$$                 self.table.addRowNames(rownames)
% $$$             end
% $$$         end
% $$$
% $$$         function setClassLabels(self, labels)
% $$$         % Set the labels for the conjugacy classes in the table
% $$$         %
% $$$         % Args:
% $$$         %   labels (cell(1,nclasses)): cell array of character string labels
% $$$             if iscell(labels)
% $$$                 if length(labels) == length(self.table.getColumnNames)
% $$$                     self.table.addColumnNames(labels)
% $$$                 end
% $$$             elseif isequal(labels, 'default')
% $$$                 colnames = cellfun(@(v) v.representative, self.classes, 'UniformOutput', false);
% $$$                 self.table.addColumnNames(colnames)
% $$$             end
% $$$         end
% $$$
% $$$         function table = pointGroupTable(self)
% $$$         % Get table with labels in chemistry notation
% $$$         %
% $$$         % Returns:
% $$$         %   table (`replab.str.Table`): character table using chemistry notation
% $$$             table = self.table;
% $$$             irrepLabels = table.getRowNames;
% $$$             for i = 1:length(irrepLabels)
% $$$                 sym = self.mulliken(self.irreps{i});
% $$$                 irrepLabels{i} = sym;
% $$$             end
% $$$             table.addRowNames(irrepLabels)
% $$$         end
% $$$
% $$$         function sizes = classSizes(self)
% $$$         % Gets the number of elements in each conjugacy class
% $$$         %
% $$$         % Returns:
% $$$         %   sizes (integer(1,nclasses)): vector of conjugacy class sizes
% $$$             sizes = cellfun(@(x) double(x.nElements), self.classes);
% $$$         end
% $$$
% $$$         function mults = multiplicities(self, rep)
% $$$         % Calculate the multiplicities of irreducible representations in rep
% $$$         %
% $$$         % Args:
% $$$         %   rep (`replab.Rep`): representation of self.group
% $$$         %
% $$$         % Returns:
% $$$         %   mults (integer(1,nirreps)): vector of multiplicities of self.irreps in the representation rep
% $$$             nirreps = length(self.irreps);
% $$$             mults = zeros(1, nirreps);
% $$$             ord = double(self.group.order);
% $$$             sizes = self.classSizes;
% $$$             repchars = self.charactersOfRepresentation(rep);
% $$$             for i = 1:nirreps
% $$$                 mults(i) = sum(sizes.*repchars.*conj(self.chars(i,:))) / ord;
% $$$             end
% $$$         end
% $$$
% $$$         function mults = tensorProdMultiplicities(self, irreps)
% $$$         % Find the multiplicities of irreducible representations in a tensor product of the irreducible representations
% $$$         %
% $$$         % Args:
% $$$         %   irreps (integer(1,\*)): vector of the locations of the irreducible
% $$$         %                           representations of which we take the tensor product
% $$$         %
% $$$         % Returns:
% $$$         %   mults (integer(1,nirreps)): vector of the multiplicities of the
% $$$         %                               irreps in the tensor representation
% $$$         %
% $$$         % Convention: to take tensor product of n copies of the same irrep, use
% $$$         %             n copies of irrep location in irreps
% $$$         %
% $$$         % Example:
% $$$         %   >>> S4 = replab.S(4);
% $$$         %   >>> s4ct = replab.CharacterTable.forPermutationGroup(S4);
% $$$         %   >>> s4ct.tensorProdMultiplicities([2,2,3])
% $$$         %       1     2     2     2     1
% $$$             tensorchars = ones(1, length(self.classes));
% $$$             for i = 1:length(irreps)
% $$$                 tensorchars = tensorchars .* self.chars(irreps(i), :);
% $$$             end
% $$$             nclass = length(self.classes);
% $$$             mults = zeros(1, nclass);
% $$$             ord = double(self.group.order);
% $$$             sizes = self.classSizes;
% $$$             for i = 1:nclass
% $$$                 mults(i) = sum(sizes.*tensorchars.*conj(self.chars(i,:))) / ord;
% $$$             end
% $$$         end



% $$$     methods (Static)
% $$$
% $$$         function sym = mulliken(rep)
% $$$         % Returns the Mulliken symbol of a representation
% $$$         %
% $$$         % Requires knowledge of the principal axis to fully determine
% $$$         %
% $$$         % Args:
% $$$         %   rep (`replab.Rep`): irreducible representation of a crystallographic group
% $$$         %
% $$$         % Returns:
% $$$         %   sym (charstring): Mulliken symbol of the irrecudible representation
% $$$             if rep.dimension == 1
% $$$                 sym = 'A/B';
% $$$             else
% $$$                 sym = 'E';
% $$$             end
% $$$         end
% $$$
% $$$         function chars = charactersOfRepresentation(rep)
% $$$         % Returns the characters of a representation
% $$$         %
% $$$         % Args:
% $$$         %   rep (`replab.Rep`): representation whose characters we want
% $$$         %
% $$$         % Returns:
% $$$         %   chars (double(1,nclasses)): vector with the character of each conjugacy class in rep
% $$$             classReps = cellfun(@(x) x.representative, rep.group.conjugacyClasses, ...
% $$$                                 'UniformOutput', false);
% $$$             images = cellfun(@(x) rep.image(x), classReps, 'UniformOutput', false);
% $$$             chars = cellfun(@trace, images);
% $$$         end
% $$$
% $$$         function table = forPermutationGroup(group)
% $$$         % Generates character table for a permutation group
% $$$         %
% $$$         % This currently assumes that the characters are integers, and the algorithm is slow.
% $$$         %
% $$$         % Args:
% $$$         %    group (`+replab.PermutationGroup`): Permutation group
% $$$         %
% $$$         % Returns:
% $$$         %    table (`+replab.CharacterTable`): Character table of group
% $$$             assert(isa(group, 'replab.PermutationGroup'));
% $$$             ord = double(group.order);
% $$$             decomp = group.naturalRep.decomposition.nice;
% $$$             irreps = decomp.components;
% $$$             classes = group.conjugacyClasses;
% $$$             k = length(classes);
% $$$             for i = 1:k
% $$$                 % verify that g^n is conjugate to g for all g \in G and n \in Z with n and (g) coprime
% $$$                 % https://math.stackexchange.com/questions/2792741/classification-of-groups-with-integer-valued-characters
% $$$                 cl = classes{i};
% $$$                 g = cl.representative;
% $$$                 eo = group.elementOrder(g);
% $$$                 for j = 2:eo-1
% $$$                     if gcd(j, eo) == 1
% $$$                         if ~cl.contains(group.composeN(g, j))
% $$$                             error('Group does not have integer characters');
% $$$                         end
% $$$                     end
% $$$                 end
% $$$             end
% $$$             ccreps = cell(1, length(classes));
% $$$             cclens = cell(1, length(classes));
% $$$             for i = 1:k
% $$$                 ccreps{i} = classes{i}.representative;
% $$$                 cclens{i} = double(classes{i}.nElements);
% $$$             end
% $$$             nirreps = length(irreps);
% $$$             chars = cell(nirreps, k);
% $$$             for i = 1:nirreps
% $$$                 irrep = irreps{i};
% $$$                 chars(i,:) = cellfun(@(x) trace(irrep.image(x)), ccreps, 'UniformOutput', false);
% $$$             end
% $$$             chars = cell2mat(chars);
% $$$             start_irreps = 1;
% $$$             while nirreps ~= k
% $$$                 for i = start_irreps:nirreps
% $$$                     for j = start_irreps:nirreps
% $$$                         new_rep = irreps{i}.kron(irreps{j});
% $$$                         ss = 0;
% $$$                         for n = 1:k
% $$$                             char = trace(new_rep.image(ccreps{n}));
% $$$                             ss = ss + cclens{n}*(char)^2;
% $$$                         end
% $$$                         if ss ~= ord
% $$$                             new_irreps = new_rep.decomposition.nice;
% $$$                             start_irreps = nirreps + 1;
% $$$                             for m = 1:new_irreps.nComponents
% $$$                                 new_irrep = new_irreps.component(m);
% $$$                                 new_char = cellfun(@(x) trace(new_irrep.image(x)), ...
% $$$                                                     ccreps, 'UniformOutput', false);
% $$$                                 new_char = cell2mat(new_char);
% $$$                                 if ~any(all(abs(chars - new_char) < 1/2, 2))
% $$$                                     irreps{nirreps + 1} = new_irrep;
% $$$                                     chars(nirreps + 1, :) = new_char;
% $$$                                     nirreps = nirreps + 1;
% $$$                                 end
% $$$                             end
% $$$                         end
% $$$                     end
% $$$                 end
% $$$             end
% $$$             chars = round(chars);
% $$$             table = replab.CharacterTable(group, irreps, classes, chars);
% $$$         end
% $$$
% $$$     end
% $$$
end