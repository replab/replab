classdef CharacterTable < replab.Obj
% Describes the real or complex character table of a finite group
%
% This contains several pieces of information. Some of these pieces, such as the character values, are essentially unique
% for a given group (up to permutation of characters and conjugacy classes). Other pieces of information, such as the
% names of conjugacy classes/irreducible representations may depend on the context (the mathematical or physical are, chemistry, etc).
%
% The character table can optionally store the explicit images of the irreducible representations. Of course, these irreps are
% defined up to the choice of basis. For example, if the character table is of type `.ComplexCharacterTable`, it can contain
% irreducible representations of real-type, those images are not necessarily written using real coefficients, even if they
% could.
%
% The arguments passed to the constructor for irreducible representations (`.irreps`) or generalized Kronecker coefficients
% (`.kronecker`) can include function handles; it that case, the corresponding objects are only constructed when they are
% requested and the computed value is cached.

    properties (SetAccess = protected)
        group % (`+replab.FiniteGroup`): Group represented by character table
        field % ({'R', 'C'}): Whether the character table is real or complex
        classes % (`.ConjugacyClasses`): Conjugacy classes of `.group`
        classNames % (cell(1,nClasses) of charstring): Names of conjugacy classes
        irrepNames % (cell(1,nIrreps) of charstring): Names of the irreducible representations/characters
        values % (`.cyclotomic` (nIrreps,nClasses)): Character values
    end

    properties (Access = protected)
        irreps_ % (cell(1, nClasses) of ``[]`` or `.RepByImages` or function_handle): Explicit representations (can contain empty values)
        kronecker_ % (integer(\*,\*,\*) or function_handle): Kronecker coefficients
    end

    methods

        function self = CharacterTable(group, field, classes, values, varargin)
        % Constructs a character table
        %
        % Args:
        %   group (`.FiniteGroup`): Group represented by character table
        %   field ({'R', 'C'}): Field over which the characters are defined
        %   classes (`.ConjugacyClasses`): Conjugacy classes of `.group`
        %   values (`.cyclotomic` (nClasses, nClasses)): Character values
        %
        % Keyword Args:
        %   irreps (cell(1,\*) of ``[]`` or `+replab.RepByImages` or function_handle, optional): Explicit matrix representations (can contain empty values)
        %   classNames (cell(1,\*) of charstring, optional): Names of conjugacy classes
        %   irrepNames (cell(1,\*) of charstring, optional): Names of irreducible representations
        %   kronecker (integer(\*,\*,\*) or function_handle, optional): Kronecker coefficients
            assert(isa(group, 'replab.FiniteGroup'));
            assert(strcmp(field, 'R') || strcmp(field, 'C'));
            assert(isa(classes, 'replab.ConjugacyClasses'));
            assert(classes.group == group);
            values = replab.cyclotomic(values);
            assert(isa(values, 'replab.cyclotomic'));
            nIrreps = size(values, 1);
            nClasses = size(values, 2);
            assert(classes.nClasses == nClasses);
            args = struct('irrepNames', {replab.CharacterTable.defaultIrrepNames(nIrreps)}, ...
                          'classNames', {replab.CharacterTable.defaultClassNames(classes.classElementOrders)}, ...
                          'irreps', {cell(1, nIrreps)}, 'kronecker', []);
            args = replab.util.populateStruct(args, varargin);
            if isempty(args.kronecker)
                args.kronecker = @() replab.ct.computeKronecker(self);
            end
            assert(length(irrepNames) == nIrreps);
            assert(length(classNames) == nClasses);
            assert(all(cellfun(@(n) ischar(n), irrepNames)));
            assert(all(cellfun(@(n) ischar(n), classNames)));
            assert(all(cellfun(@(x) isempty(x) || isa(x, 'replab.Rep') || isa(x, 'function_handle'), args.irreps)));
            assert(isdouble(args.kronecker) || isa(args.kronecker, 'function_handle'));
            self.group = group;
            self.field = field;
            self.classes = classes;
            self.values = values;
            self.irrepNames = args.irrepNames;
            self.classNames = args.classNames;
            self.irreps_ = args.irreps;
        end

        function K = kronecker(self)
        % Returns the Kronecker coefficients corresponding to this character table
        %
        % This returns an integer matrix $K$ such that $K(i,j,k)$ is the multiplicity of the $i$-th irrep in
        % the product of the $j$-th and $k$-th irrep.
        %
        % Returns:
        %   integer(\*,\*,\*): Kronecker coefficients
            if isa(self.kronecker_, 'function_handle')
                f = self.kronecker_;
                self.kronecker_ = f();
            end
            K = self.kronecker_;
        end

        function b = overR(self)
        % Returns whether the character table is defined over the real field
        %
        % Returns:
        %   logical: True if the character table is defined over the real field
            b = strcmp(self.field, 'R');
        end

        function b = overC(self)
        % Returns whether the character table is defined over the complex field
        %
        % Returns:
        %   logical: True if the character table is defined over the complex field
            b = strcmp(self.field, 'C');
        end

    end

    methods % Conjugacy classes

        function ind = identityConjugacyClassIndex(self)
        % Returns the index of the conjugacy class containing the identity
        %
        % Returns:
        %   integer: Index of (the column containing) the identity conjugacy class
            ind = self.classes.classIndex(self.group.identity);
        end

        function n = nClasses(self)
        % Returns the number of conjugacy classes in this character table
        %
        % Returns:
        %   integer: Number of conjugacy classes
            n = self.classes.nClasses;
        end

    end

    methods % Characters

        function c = character(self, ind)
        % Returns an irreducible character in this character table
        %
        % Args:
        %   ind (integer): Index of the character
        %
        % Returns:
        %   `.Character`: Indexed character
            c = replab.Character(self.classes, self.values(ind, :));
        end

        function C = characters(self)
        % Returns the characters in this table
        %
        % Returns:
        %   cell(1,\*) of `.Character`: List of characters contained in this character table
            C = arrayfun(@(i) self.character(i), 1:self.nCharacters, 'uniform', 0);
        end

        function ind = linearCharacterIndices(self)
        % Returns the indices of the linear characters
        %
        % The linear characters correspond to representation of dimension ``1``.
        %
        % Returns:
        %   integer(1,\*): Indices of the linear characters
            ind = find(self.values(:, self.identityConjugacyClassIndex) == 1);
        end

        function n = nCharacters(self)
        % Returns the number of characters in this character table
        %
        % Returns:
        %   integer: Number of characters
            n = self.nClasses;
        end

    end

    methods % Representations

        function dec = decomposeRep(self, rep)
        % Computes the decomposition of a representation into irreducibles using the explicit irreps contained in this table
        %
        % Args:
        %   rep (`.Rep`): Representation to decompose, with ``rep.field == self.field``
        %
        % Returns:
        %   `.Irreducible`: Irreducible decomposition
            assert(rep.field == self.field, 'The representation must be defined on the same field as the character table');
            assert(self.group == rep.group);
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

        function n = nIrreps(self)
        % Returns the number of
            n = self.;
        end


    end

    methods % Transformations of the character table

% $$$         function ct = directProduct(self, ct2)
% $$$         % Returns the direct product of character tables
% $$$         %
% $$$         % Args:
% $$$         %   ct2 (`+replab.ComplexCharacterTable`): character table with which to take direct product
% $$$         %
% $$$         % Returns:
% $$$         %   ct (`+replab.ComplexCharacterTable`): The direct product of the character tables
% $$$             ct = replab.ct.directProduct(self, ct2);
% $$$         end

        function ct = forIsomorphicGroup(self, group)
        % Returns the character table updated to describe the representations of an isomorphic group
        %
        % Args:
        %   group (`.FiniteGroup`): Group isomorphic to `.group`
        %
        % Returns:
        %   `.CharacterTable`: Updated character table


        % TODO
        % TODO
        % TODO
        % TODO
        % TODO
        % TODO
        end

    end

    methods % Implementations

    end

    methods (Static)

    end

end
