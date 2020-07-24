classdef AbstractGroup < replab.NiceFiniteGroup
% Describes a group whose elements are written as product of generators
%
% The group is the quotient of a free group by the normal closure of a set of relators. The relators are
% written using words in the free group.

    properties (Access = protected)
        groupId % (integer): Unique group id
    end

    properties (SetAccess = protected)
        names % (cell(1,\*) of charstring): Generator names
        permutationGroup % (`.PermutationGroup`): Permutation group realization of this abstract group
    end

    methods (Static)

        function A = make(names, relators)
            gens = replab.fp.permutationGeneratorsForRelators(names, relators);
            pg = replab.PermutationGroup.of(gens{:});
            A = replab.AbstractGroup(names, pg, relators);
        end

        function [A varargout] = parsePresentation(str)
        % Creates an abstract group from a presentation string
        %
        % Returns the finite group generators as additional output arguments.
        %
        % Example:
        %   >>> [G, x] = replab.AbstractGroup.parsePresentation('< x | x^3 = 1 >');
        %
        % Args:
        %   str (charstring): Single-line description string
        %
        % Returns:
        %   `.AbstractGroup`: The parsed abstract group
            [ok, names, relatorLetters] = replab.fp.Parser.parsePresentation(str);
            assert(ok, 'Error in given presentation string');
            relatorLetters = cellfun(@(r) replab.fp.reduceLetters(r), relatorLetters, 'uniform', 0);
            mask = cellfun(@isempty, relatorLetters);
            relatorLetters = relatorLetters(~mask);
            relators = cell(1, length(relatorLetters));
            for i = 1:length(relatorLetters)
                relators{i} = replab.fp.printLetters(relatorLetters{i}, names, ' ');
            end
            A = replab.AbstractGroup.make(names, relators);
            if nargout > 1
                for i = 1:length(A.nGenerators)
                    varargout{i} = A.generator(i);
                end
            end
        end

    end

    methods

        function self = AbstractGroup(names, permutationGroup, relators)
        % Creates an abstract group from generator names, images and optional relators
        %
        % Args:
        %   names (cell(1,\*) of charstring): Generator names
        %   permutationGroup (`.PermutationGroup`): Permutation group realization of this abstract group
        %   relators (cell(1,\*) of charstring, optional): Relators
            self.type = self;
            self.groupId = replab.globals.nextUniqueId;
            self.identity = '1';
            self.names = names;
            self.generators = names;
            self.permutationGroup = permutationGroup;
            if nargin >= 3
                self.cache('relators', relators, '=');
            end
        end

        function r = relators(self)
            r = self.cached('relators', @() self.computeRelators);
        end

        function r = computeRelators(self)
            r = replab.fp.relatorsForPermutationGroup(self.permutationGroup, self.names);
        end

        function letters = toLetters(self, word)
        % Parses word letters from word as a string
        %
        % Example:
        %   >>> A = replab.AbstractGroup({'x'}, {[2 3 1]}, {'x^3'});
        %   >>> isequal(A.toLetters('x^2'), [1 1])
        %       1
        %
        % Args:
        %   word (charstring): Word as a string
        %
        % Returns:
        %   integer(1,\*): Word letters
        %
        % Raises:
        %   An error if the string is malformed
            [ok, tokens] = replab.fp.Parser.lex(word, self.names);
            assert(ok, 'Unknown tokens in string');
            [pos, letters] = replab.fp.Parser.word(tokens, 1);
            assert(pos > 0, 'Malformed word');
            assert(tokens(1, pos) == replab.fp.Parser.types.END, 'Badly terminated word');
        end

        function word = fromLetters(self, letters)
        % Prints a word formed of letters as a string
        %
        %   >>> A = replab.AbstractGroup({'x'}, {[2 3 1]}, {'x^3'});
        %   >>> A.fromLetters([1 1])
        %       'x^2'
        %
        % Args:
        %   letters (integer(1,\*)): Word letters
        %
        % Returns:
        %   charstring: Word as a string
            word = replab.fp.printLetters(letters, self.names, ' ');
        end

        function img = computeImage(self, word, target, targetGeneratorImages)
        % Computes the image of this word using the given generator images
        %
        % Does not verify the validity of the implied homomorphism.
        %
        % Args:
        %   word (charstring): Word
        %   target (`+replab.Group`): Target group
        %   targetGeneratorImages (cell(1,\*) of elements of ``target``): Images of the generators of this group
        %
        % Returns:
        %   permutation`: Computed image
            letters = self.toLetters(word);
            img = target.identity;
            for i = 1:length(letters)
                l = letters(i);
                if l > 0
                    img = target.compose(img, targetGeneratorImages{l});
                else
                    img = target.composeWithInverse(img, targetGeneratorImages{-l});
                end
            end
        end


        function l = imagesDefineMorphism(self, target, targetGeneratorImages)
        % Checks whether the given images satisfy the relations of the presentation of this group
        %
        % If it returns true, it means those images describe a valid homomorphism from this `.AbstractGroup`
        % to the given target group.
        %
        % Args:
        %   target (`+replab.Group`): Target group
        %   targetGeneratorImages (cell(1,\*) of elements of ``target``): Images of the generators of this group
        %
        % Returns:
        %   logical: True if the images verify the presentation
            nR = length(self.relators);
            for i = 1:nR
                r = self.relators{i};
                if ~target.isIdentity(self.computeImage(r, target, targetGeneratorImages))
                    l = false;
                    return
                end
            end
            l = true;
        end

    end

    methods % Implementations

        function res = eq(self, rhs)
            res = (self.groupId == rhs.groupId);
        end

        function res = ne(self, rhs)
            res = self.groupId ~= rhs.groupId;
        end

        % Domain

        function b = eqv(self, x, y)
            b = all(self.niceImage(x) == self.niceImage(y));
        end

        function l = laws(self)
            l = replab.AbstractGroupLaws(self);
        end

        % Monoid

        function z = compose(self, x, y)
            xl = self.toLetters(x);
            yl = self.toLetters(y);
            zl = replab.fp.composeLetters(xl, yl);
            z = self.fromLetters(zl);
        end

        % Group

        function z = inverse(self, x)
            xl = self.toLetters(x);
            zl = replab.fp.inverseLetters(xl);
            z = self.fromLetters(zl);
        end

        % NiceFiniteGroup

        function perm = niceImage(self, word)
        % Computes the image of this word using the given generator images
        %
        % Does not verify the validity of the implied homomorphism.
        %
        % Args:
        %   word (charstring): Word
        %
        % Returns:
        %   permutation`: Computed image
            letters = self.toLetters(word);
            pg = self.permutationGroup;
            perm = pg.identity;
            for i = 1:length(letters)
                l = letters(i);
                if l > 0
                    perm = pg.compose(perm, pg.generator(l));
                else
                    perm = pg.composeWithInverse(perm, pg.generator(-l));
                end
            end
        end

    end

end
