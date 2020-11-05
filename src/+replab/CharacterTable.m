classdef CharacterTable < replab.Obj
% Describes the character table of a group
%
% Example:
%   >>> G = replab.PermutationGroup.dihedral(3);
%   >>> G.characterTable
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
        characters % (`.cyclotomic` (nIrreps,nClasses)): Character values
    end

    properties (Access = protected)
        irreps_ % (cell(1, nClasses) of ``[]`` or `.RepByImages`): Explicit matrix representations (can contain empty values)
    end

    methods (Static)

        function ct = forPermutationGroup(G)
            ct = replab.sym.PermutationCharacterTable(G);
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
            self.irreps_ = args.irreps;
        end

        function dec = decomposition(self, rep)
        % Computes the decomposition of a representation into irreducibles using the characters
        %
        % Args:
        %   rep (`.Rep`): Representation to decompose, must be complex
        %
        % Returns:
        %   `.Irreducible`: Irreducible decomposition
            assert(rep.overC, 'Can only decompose complex representations');
            nI = self.nIrreps;
            I = cell(1, nI);
            for i = 1:self.nIrreps
                irrep = self.irrep(i);
                B = replab.irreducible.findIsotypicBasis(rep, irrep);
                sub = rep.subRep(B);
                E = sub.E_internal;
                d = irrep.dimension;
                m = size(B, 2) / d;
                irreds = cell(1, m);
                for j = 1:m
                    comp = rep.subRep(B(:, (j-1)*d + (1:d)));
                    comp.isIrreducible = true;
                    irreds{j} = comp;
                end
                I{i} = replab.HarmonizedIsotypic(rep, irreds, E);
            end
            dec = replab.Irreducible(rep, I);
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

        function ind = trivialCharacterIndex(self)
        % Returns the row index of the trivial character
            ind = 1:self.nIrreps;
            for i = 1:self.nClasses
                ind = ind(self.characters(ind, i) == 1);
            end
            assert(length(ind) == 1);
        end

        function ind = identityConjugacyClassIndex(self)
        % Returns the column index of the identity conjugacy class
            ind = self.classes.classIndex(self.group.identity);
        end

        function ind = linearCharacterIndices(self)
        % Returns the indices of the linear characters
        %
        % The linear characters correspond to representation of dimension ``1``.
        %
        % Returns:
        %   integer(1,\*): Indices of the linear characters
            ind = find(self.characters(:, self.identityConjugacyClassIndex) == 1);
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
            ct = replab.ct.directProduct(self, ct2);
        end

        function n = nClasses(self)
            n = self.classes.nClasses;
        end

        function n = nIrreps(self)
            n = self.nClasses;
        end

        function r = irreps(self)
        % Returns a cell vector of the irreducible representations corresponding to the characters in this table
        %
        % Returns:
        %   cell(1,nIrreps) of `.RepByImages`: Irreducible representations
            r = arrayfun(@(i) self.irrep(i), 1:self.nIrreps, 'uniform', 0);
        end

        function r = irrep(self, ind)
        % Returns the irreducible representation that corresponds to the character of given index
        %
        % Args:
        %   ind (integer): Index of the character
        %
        % Returns:
        %   `.RepByImages`: An irreducible representation with coefficients in the cyclotomic field
            r = self.irreps_{ind};
            if isempty(r)
                error('I do not know how to compute the %d-th irrep of this group', ind);
            end
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

    methods (Access = protected)

        function mul12 = multiplicityProduct(self, mul1, mul2)
            n = self.classes.nClasses;
            K = self.kronecker;
            mul12 = reshape(K, n, n*n) * reshape(mul1(:) * mul2(:)', n*n, 1);
            mul12 = mul12(:)';
        end

    end

    methods

        function mults = multiplicities(self, arg)
        % Calculate the multiplicities of the irreducible characters in this table in a given representation or character
        %
        % The ordering of coefficients corresponds to the order of irreducible representations in this character table.
        %
        % Note that this method is optimized when the representation is a tensor product.
        %
        % Example:
        %   >>> G = replab.PermutationGroup.dihedral(3);
        %   >>> ct = G.characterTable;
        %   >>> rep2 = ct.irreps{2};
        %   >>> rep3 = ct.irreps{3};
        %   >>> rep = kron(rep2, rep3);
        %   >>> isequal(ct.multiplicities(rep), [0 0 1])
        %       1
        %
        % Example:
        %   >>> G = replab.S(5);
        %   >>> ct = G.characterTable;
        %   >>> S5 = ct.group;
        %   >>> isequal(ct.multiplicities(S5.naturalRep), [1 0 1 0 0 0 0])
        %       1
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
            elseif isa(arg, 'replab.rep.TensorRep')
                factorM = cellfun(@(f) self.multiplicities(f), arg.factors, 'uniform', 0);
                mults = factorM{1};
                for i = 2:length(factorM)
                    mults = self.multiplicityProduct(mults, factorM{i});
                end
            elseif isa(arg, 'replab.Rep')
                for i = 1:n
                    mults(i) = self.character(i).dotRep(arg);
                end
            else
                error('Invalid argument');
            end
            mults = round(mults);
        end


        function res = imap(self, f)
        % Maps the character table under an isomorphism
        %
        % Example:
        %   >>> D6a = replab.PermutationGroup.of([3 2 1], [2 3 1]);
        %   >>> D6b = replab.PermutationGroup.of([1 4 3 2], [1 3 4 2]);
        %   >>> f = D6a.isomorphismByImages(D6b, 'preimages', D6a.generators, 'images', D6b.generators);
        %   >>> Ca = D6a.characterTable;
        %   >>> Cb = Ca.imap(f);
        %   >>> Cb.laws.checkSilent;
        %
        % Args:
        %   f (`.FiniteIsomorphism`): Isomorphism with ``self.group.isSubgroupOf(f.source)``
        %
        % Returns:
        %   `.CharacterTable`: The character table of the subgroup in the image of the isomorphism
            if self.group.order < f.source.order
                f = f.restrictedSource(self.group);
            end
            classes1 = self.classes.imap(f);
            group1 = f.target;
            characters1 = self.characters;
            irreps1 = cell(1, self.nIrreps);
            for i = 1:self.nIrreps
                if ~isempty(self.irreps{i})
                    irreps1{i} = self.irreps{i}.imap(f);
                end
            end
            res = replab.CharacterTable(group1, classes1, characters1, 'irreps', irreps1, 'classNames', self.classNames, 'irrepNames', self.irrepNames);
            if self.inCache('kronecker')
                res.cache('kronecker', self.kronecker, 'error');
            end
        end

    end

end
