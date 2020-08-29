classdef CharacterTable < replab.Obj
% Describes the character table of a group
%
% Example:
%   >>> replab.CharacterTable.dihedral(3)
%       Class  1a   3a   2a
%        Size   1    2    3
%
%         X.1  1    1    1
%         X.2  1    1   -1
%         X.3  2   -1    0
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

    methods (Static)

        function ct = cyclic(n)
            ct = replab.ct.CyclicCharacterTable(n);
        end

        function ct = dihedral(n)
            ct = replab.ct.DihedralCharacterTable(n);
        end

        function ct = A5
            ct = replab.ct.A5CharacterTable;
        end

        function ct = A4
            ct = replab.ct.A4CharacterTable;
        end

        function ct = S4
            ct = replab.ct.S4CharacterTable;
        end

        function ct = forPermutationGroup(G)
            ct = replab.ct.PermutationCharacterTable(G);
        end

    end

    methods (Access = protected)

        function K = computeKronecker(self)
            n = self.classes.nClasses;
            K = zeros(n,n,n);
            for j = 1:n
                cj = self.character(j);
                for k = j:n
                    ck = self.character(k);
                    K(:,j,k) = self.multiplicities(cj*ck);
                    K(:,k,j) = K(:,j,k);
                end
            end
        end

    end

    methods

        function self = CharacterTable(group, classes, characters, varargin)
        % Constructs a character table
        %
        % Args:
        %   group (`.FiniteGroup`): Group represented by character table
        %   classes (`.ConjugacyClasses`): Conjugacy classes of `.group`
        %   characters (`.cyclotomic` (nClasses, nClasses)): Character values
        %
        % Keyword Args:
        %   irreps (cell(1,\*) of ``[]`` or `+replab.RepByImages`, optional): Explicit matrix representations (can contain empty values)
        %   classNames (cell(1,\*) of charstring, optional): Names of conjugacy classes
        %   irrepNames (cell(1,\*) of charstring, optional): Names of irreducible representations
        %   kronecker (integer(\*,\*,\*), optional): Kronecker coefficients
            nIrreps = size(characters, 1);
            nClasses = classes.nClasses;
            assert(size(characters, 2) == nClasses);
            assert(nIrreps == nClasses);
            args = struct('irrepNames', {replab.CharacterTable.defaultIrrepNames(nIrreps)}, ...
                          'classNames', {replab.CharacterTable.defaultClassNames(classes.classElementOrders)}, ...
                          'irreps', {cell(1, nIrreps)}, 'kronecker', []);
            args = replab.util.populateStruct(args, varargin);
            if ~isempty(args.kronecker)
                self.cache('kronecker', args.kronecker);
            end
            self.group = group;
            self.classes = classes;
            self.characters = characters;
            self.irrepNames = args.irrepNames;
            self.classNames = args.classNames;
            self.irreps = args.irreps;
        end

        function K = kronecker(self)
        % Returns the Kronecker coefficients corresponding to this character table
        %
        % This returns an integer matrix $K$ such that $K(i,j,k)$ is the multiplicity of the $i$-th irrep in
        % the product of the $j$-th and $k$-th irrep.
        %
        % Returns:
        %   integer(\*,\*,\*): Kronecker coefficients
            K = self.cached('kronecker', @() self.computeKronecker);
        end

        function c = character(self, ind)
        % Returns an irreducible character in this character table
            c = replab.Character(self.classes, self.characters(ind, :));
        end

        function ct = directProduct(self, ct2)
        % Returns the direct product of character tables
        %
        % Args:
        %   ct2 (`+replab.CharacterTable`): character table with which to take direct product
        %
        % Returns:
        %   ct (`+replab.CharacterTable`): The direct product of the character tables
            new_group = self.group.directProduct(ct2.group);
            % New characters are kronecker product of character matrices
            A = self.characters;
            B = ct2.characters;
            new_chars = replab.cyclotomic.zeros(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2));
            for i = 0:size(A, 1)-1
                for j = 0:size(A, 2)-1
                    new_chars(i * size(B,1) + 1:(i+1) * size(B,1), j * size(B,2) + 1:(j+1) * size(B,2)) = A(i+1, j+1) * B;
                end
            end
            % New conjugacy classes are cartesian product of input
            c1 = cellfun(@(c) c.representative, self.classes.classes, 'UniformOutput', false);
            c2 = cellfun(@(c) c.representative, ct2.classes.classes, 'UniformOutput', false);
            new_classes = cell(1, length(c1) * length(c2));
            for i = 0:length(c1) - 1
                new_classes(i*length(c2)+1:(i+1)*length(c2)) = cellfun(@(x) {c1{i+1}, x}, c2, 'UniformOutput', false);
            end
            classarray = cellfun(@(r) new_group.conjugacyClass(r), new_classes, 'UniformOutput', false);
            new_classes = replab.ConjugacyClasses(new_group, classarray);
            % New irreps are direct products of input irreps
            new_irreps = cell(1, length(self.irreps) * length(ct2.irreps));
            for i = 1:length(self.irreps)
                for j = 1:length(ct2.irreps)
                    B = ct2.irreps{j}.image(ct2.group.identity);
                    new_images1 = cell(1, self.group.nGenerators);
                    for k = 1:self.group.nGenerators
                        A = self.irreps{i}.images_internal{k};
                        % Kronecker product of cyclotomics
                        new_image = replab.cyclotomic.zeros(size(A,1) * size(B, 1), size(A, 2) * size(B, 2));
                        for m = 0:size(A, 1)-1
                            for n = 0:size(A, 2)-1
                                new_image(m*size(B,1)+1:(m+1)*size(B,1), n*size(B,2)+1:(n+1)*size(B,2)) = A(m+1, n+1) * B;
                            end
                        end
                        new_images1{k} = new_image;
                    end
                    new_images2 = cell(1, ct2.group.nGenerators);
                    A = self.irreps{i}.image(self.group.identity);
                    for k = 1:ct2.group.nGenerators
                        B = ct2.irreps{j}.images_internal{k};
                        % Kronecker product of cyclotomics
                        new_image = replab.cyclotomic.zeros(size(A,1) * size(B, 1), size(A, 2) * size(B, 2));
                        for m = 0:size(A, 1)-1
                            for n = 0:size(A, 2)-1
                                new_image(m*size(B,1)+1:(m+1)*size(B,1), n*size(B,2)+1:(n+1)*size(B,2)) = A(m+1, n+1) * B;
                            end
                        end
                        new_images2{k} = new_image;
                    end
                    new_dim = double(self.irreps{i}.dimension) * double(ct2.irreps{j}.dimension);
                    new_irrep = new_group.repByImages('C', new_dim, [new_images1, new_images2]);
                    new_irreps{j + (i-1)*length(ct2.irreps)} = new_irrep;
                end
            end
            ct = replab.CharacterTable(new_group, new_classes, new_chars, 'irreps', new_irreps);
        end

        function n = nClasses(self)
            n = self.classes.nClasses;
        end

        function n = nIrreps(self)
            n = self.nClasses;
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

        % Str

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

        % Obj

        function l = laws(self)
            l = replab.laws.CharacterTableLaws(self);
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

        function mults = multiplicities(self, arg)
        % Calculate the multiplicities of the irreducible characters in this table in a given representation or character
        %
        % The ordering of coefficients corresponds to the order of irreducible representations in this character table.
        %
        % Args:
        %   arg (`.Character` or `.Rep`): Character or representation of `.group`
        %
        % Returns:
        %    (integer(1,\*)): Multiplicities of irreducible representations
            n = self.nIrreps;
            mults = zeros(1, n);
            if isa(arg, 'replab.Character')
                for i = 1:n
                    mults(i) = self.character(i).dot(arg);
                end
            elseif isa(arg, 'replab.Rep')
                for i = 1:n
                    mults(i) = self.character(i).dotRep(rep);
                end
            else
                error('Invalid argument');
            end
            mults = round(mults);
        end

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
% $$$
end
