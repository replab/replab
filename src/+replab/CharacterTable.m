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
        %                                         or ``'default'`` to use default labels
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
        % Get table with labels in chemistry notation
        %
        % Returns:
        %   table (`replab.str.Table`): character table using chemistry notation
            table = self.table;
            irrepLabels = table.getRowNames;
            for i = 1:length(irrepLabels)
                sym = self.mulliken(self.irreps{i});
                irrepLabels{i} = sym;
            end
            table.addRowNames(irrepLabels)
        end

        function sizes = classSizes(self)
        % Gets the number of elements in each conjugacy class
        %
        % Returns:
        %   sizes (integer(1,nclasses)): vector of conjugacy class sizes
            sizes = cellfun(@(x) double(x.size), self.classes);
        end

        function mults = multiplicities(self, rep)
        % Calculate the multiplicities of irreducible representations in rep
        %
        % Args:
        %   rep (`replab.Rep`): representation of self.group
        %
        % Returns:
        %   mults (integer(1,nirreps)): vector of multiplicities of self.irreps in the representation rep
            nirreps = length(self.irreps);
            mults = zeros(1, nirreps);
            ord = double(self.group.order);
            sizes = self.classSizes;
            repchars = self.charactersOfRepresentation(rep);
            mults = (sizes.*repchars)*conj(self.chars(i,:))' / ord;
        end

        function mults = tensorProdMultiplicities(self, irreps)
        % Find the multiplicities of irreducible representations in a tensor product of the irreducible representations
        %
        % Args:
        %   irreps (integer(1,\*)): vector of the locations of the irreducible
        %                           representations of which we take the tensor product
        %
        % Returns:
        %   mults (integer(1,nirreps)): vector of the multiplicities of the
        %                               irreps in the tensor representation
        %
        % Convention: to take tensor product of n copies of the same irrep, use
        %             n copies of irrep location in irreps
        %
        % Example:
        %   >>> S4 = replab.S(4);
        %   >>> s4ct = replab.CharacterTable.forPermutationGroup(S4);
        %   >>> s4ct.tensorProdMultiplicities([2,2,3])
        %       1     2     2     2     1
            tensorchars = ones(1, length(self.classes));
            for i = 1:length(irreps)
                tensorchars = tensorchars .* self.chars(irreps(i), :);
            end
            nclass = length(self.classes);
            mults = zeros(1, nclass);
            ord = double(self.group.order);
            sizes = self.classSizes;
            mults(i) = sizes.*tensorchars*conj(self.chars')/ ord;
        end


    end

    methods (Static)

        function sym = mulliken(rep)
        % Returns the Mulliken symbol of a representation
        %
        % Requires knowledge of the principal axis to fully determine
        %
        % Args:
        %   rep (`replab.Rep`): irreducible representation of a crystallographic group
        %
        % Returns:
        %   sym (charstring): Mulliken symbol of the irrecudible representation
            if rep.dimension == 1
                sym = 'A/B';
            else
                sym = 'E';
            end
        end

        function chars = charactersOfRepresentation(rep)
        % Returns the characters of a representation
        %
        % Args:
        %   rep (`replab.Rep`): representation whose characters we want
        %
        % Returns:
        %   chars (double(1,nclasses)): vector with the character of each conjugacy class in rep
            classReps = cellfun(@(x) x.representative, rep.group.conjugacyClasses, ...
                                'UniformOutput', false);
            images = cellfun(@(x) rep.image(x), classReps, 'UniformOutput', false);
            chars = cellfun(@trace, images);
        end

        function permCT = permutationCharTable(permGroup)
            n = permGroup.domainSize;
            allData = replab.sym.CTData.instance(n);
            nData = allData(n);
            irreps = cellFun(@(part) num2str(part),nData.partitionList,'UniformOutput', false);
            classes = replab.sym.findConjClasses(permGroup);
            chars = permCTValues(n);
            permCT = CharacterTable(permGroup, irreps, classes, chars);
            
            function CT = permCTValues(n)
                nParts = nData.nParts;
                CTunNormed = zeros(nParts);
                for row = 1:nParts
                   partition = nData.partitionList{row};
                   lin = replab.sym.SymPoly.powerToMonomial(partition,n);
                   CTunNormed(:,row) = lin.coeffs';
                end
                CT = grahamSchmidt(CTunNormed,@(x,y) round(sum(x.*y.*nData.conjSizes/nData.fact)));
                CT = flip(round(CT));
                function mat = grahamSchmidt(mat,innerProd)
                   d = nData.nParts;
                   for i = d:-1:2
                       for j = 1:(i-1)
                           mat(j,:) = mat(j,:)-innerProd(mat(j,:),mat(i,:))*mat(i,:);
                       end
                   end
           end
            end
        end



    end
end