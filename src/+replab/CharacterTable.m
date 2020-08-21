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
% Instances of `.CharacterTable` are immutable.

    properties (SetAccess = protected)
        group % (`+replab.FiniteGroup`): Group represented by character table
        classes % (`.ConjugacyClasses`): Conjugacy classes of `.group`
        classNames % (cell(1,nClasses) of charstring): Names of conjugacy classes
        irrepNames % (cell(1,nIrreps) of charstring): Names of the irreducible representations/characters
        characters % (`.cyclotomic` (nClasses, nClasses)): Character values
        irreps % (cell(1, nClasses) of ``[]`` or `.RepByImages`): Explicit matrix representations (can contain empty values)
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

        function self = CharacterTable(group, classes, characters, varargin)
        % Constructs a character table
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Group represented by character table
        %   classes (`.ConjugacyClasses`): Conjugacy classes of `.group`
        %   characters (`.cyclotomic` (nClasses, nClasses)): Character values
        %
        % Keyword Args:
        %   irreps (cell(1,\*) of ``[]`` or `+replab.RepByImages`): Explicit matrix representations (can contain empty values)
        %   classNames (cell(1,\*) of charstring): Names of conjugacy classes
        %   irrepNames (cell(1,\*) of charstring): Names of irreducible representations
            nIrreps = size(characters, 1);
            nClasses = classes.nClasses;
            assert(size(characters, 2) == nClasses);
            assert(nIrreps == nClasses);
            args = struct('irrepNames', {replab.CharacterTable.defaultIrrepNames(nIrreps)}, ...
                          'classNames', {replab.CharacterTable.defaultClassNames(classes.classElementOrders)}, ...
                          'irreps', {cell(1, nIrreps)});
            args = replab.util.populateStruct(args, varargin);
            self.group = group;
            self.classes = classes;
            self.characters = characters;
            self.irrepNames = args.irrepNames;
            self.classNames = args.classNames;
            self.irreps = args.irreps;
        end

    end

    methods (Static)

        function names = defaultIrrepNames(n)
        % Constructs default names for irreducible representations
        %
        % Args:
        %   n (integer): Number of irreducible representations
        %
        % Returns:
        %   cell(1,\*) of charstring: Irreducible representation names
            names = arrayfun(@(i) sprintf('X.%d', i), 1:n, 'uniform', 0);
        end

        function names = defaultClassNames(elementOrders)
        % Constructs default names for conjugacy classes
        %
        % Args:
        %   elementOrders (integer(1,\*)): Orders of the conjugacy classes elements
        %
        % Returns:
        %   cell(1,\*) of charstring: Conjugacy class names
            orders = [];
            numbers = [];
            nClasses = length(elementOrders);
            names = cell(1, nClasses);
            for i = 1:nClasses
                o = elementOrders(i);
                j = find(orders == o);
                if isempty(j)
                    j = length(orders) + 1;
                    orders(j) = o;
                    numbers(j) = 1;
                else
                    numbers(j) = numbers(j) + 1;
                end
                names{i} = sprintf('%d%s', o, char('a' + numbers(j) - 1));
            end
        end

    end

    methods % Implementations

        function s = headerStr(self)
            group_str = replab.headerStr(self.group);
            if length(group_str) > 1
                group_str = [lower(group_str(1)), group_str(2:end)];
            end
            s = sprintf(['Character table for ', group_str]);
        end

        function lines = longStr(self, maxRows, maxColumns)
            switch replab.globals.formatCharacterTable
              case 'gap'
                lines = self.gapLongStr(maxRows, maxColumns);
              case 'plain'
                lines = self.plainLongStr(maxRows, maxColumns);
              otherwise
                error('Unknown format %s', replab.globals.formatCharacterTable);
            end
        end

    end

    methods (Access = protected)

        function strings = centralizerSizeTable(self)
            primes = unique(double(factor(self.group.order)));
            nC = self.classes.nClasses;
            strings = cell(length(primes), nC+1);
            for i = 1:length(primes)
                strings{i,1} = sprintf('%d', primes(i));
            end
            for i = 1:nC
                c = self.classes.classes{i};
                f = double(factor(c.representativeCentralizer.order));
                for j = 1:length(primes)
                    s = sum(f == primes(j));
                    if s == 0
                        strings{j,i+1} = '.';
                    else
                        strings{j,i+1} = sprintf('%d', sum(f == primes(j)));
                    end
                end
            end
        end

        function lines = gapLongStr(self, maxRows, maxColumns)
            ct = replab.str.CyclotomicTable(self.characters);
            pp = self.classes.powerMapDefaultPrimes;
            m = self.classes.powerMaps(pp);
            nC = self.classes.nClasses;
            powerMaps = cell(length(pp)+1, nC + 1);
            powerMaps{1,1} = '';
            powerMaps(1,2:end) = self.classNames;
            for i = 1:length(pp)
                powerMaps{i+1,1} = sprintf('%dP', pp(i));
                for j = 1:nC
                    powerMaps{i+1,j+1} = self.classNames{m(i,j)};
                end
            end
            chars = horzcat(self.irrepNames(:), ct.strings);
            cst = self.centralizerSizeTable;
            sep = repmat({''}, 1, nC+1);
            t = replab.str.Table(vertcat(cst, sep, powerMaps, sep, chars), 'colAlign', repmat('r', 1, nC+1));
            lines1 = strsplit(t.format(maxRows, maxColumns), '\n');
            lines2 = arrayfun(@(i) sprintf(' %s = %s', ct.variables{i}, num2str(ct.values(i))), 1:length(ct.variables), 'uniform', 0);
            lines = vertcat(lines1(:), {''}, lines2(:));
        end

        function lines = plainLongStr(self, maxRows, maxColumns)
            ct = replab.str.CyclotomicTable(self.characters);
            nC = self.classes.nClasses;
            header = cell(2, nC + 1);
            header{1,1} = 'Class';
            header(1,2:end) = self.classNames;
            header{2,1} = 'Size';
            header(2,2:end) = cellfun(@(s) strtrim(num2str(s)), self.classes.classSizes, 'uniform', 0);
            chars = horzcat(self.irrepNames(:), ct.strings);
            sep = repmat({''}, 1, nC+1);
            t = replab.str.Table(vertcat(header, sep, chars), 'colAlign', repmat('r', 1, nC+1));
            lines1 = strsplit(t.format(maxRows, maxColumns), '\n');
            lines2 = arrayfun(@(i) sprintf(' %s = %s', ct.variables{i}, num2str(ct.values(i))), 1:length(ct.variables), 'uniform', 0);
            lines = vertcat(lines1(:), {''}, lines2(:));
        end

    end

    methods

% $$$         function imap(self, f, imageGroup, preserveLexOrder)
% $$$         % Maps the conjugacy classes under an isomorphism
% $$$         %
% $$$         % Args:
% $$$         %   f (`.FiniteIsomorphism`): Isomorphism with ``self.group.isSubgroupOf(f.source)``
% $$$         %   imageGroup (`.FiniteGroup`, optional): Image of `.group` under ``f``, default ``[]`` (recompute)
% $$$         %   preserveLexOrder (logical, optional): Whether the isomorphism preserves the lexicographic order of group elements, default false
% $$$         %
% $$$         % Returns:
% $$$         %   `.CharacterTable`: The character table of the subgroup in the image of the isomorphism
% $$$             if nargin < 3 || isempty(imageGroup)
% $$$                 imageGroup = f.imageGroup(self.group);
% $$$             end
% $$$             if nargin < 4 || isempty(preserveLexOrder)
% $$$                 preserveLexOrder = false;
% $$$             end
% $$$             classes1 = cellfun(@(c) c.imap(f, imageGroup, preserveLexOrder), self.classes, 'uniform', 0);
% $$$             c1 = replab.ConjugacyClasses(imageGroup, classes1);
% $$$
% $$$         end


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