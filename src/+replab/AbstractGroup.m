classdef AbstractGroup < replab.NiceFiniteGroup
% Describes a group whose elements are written as products of generators
%
% The elements of that group obey an equivalence relation, which can be given in two different ways:
%
% - A list of relations form a group presentation (see `<https://en.wikipedia.org/wiki/Presentation_of_a_group>`_
%
% - The equivalence relation is given by an homomorphism into an explicit realization of the group.
%
% We now give two corresponding examples.
%
% Example:
%   >>> [G, x] = replab.AbstractGroup.parsePresentation('< x | x^3 = 1 >');
%   >>> x
%       x =
%       'x'
%   >>> G.compose(x, x)
%       'x^2'
%   >>> G.order
%       3
%
% Example:
%   >>> G = replab.AbstractGroup({'s', 't'}, replab.S(3));
%   >>> G.order
%       6
%   >>> G.elements.at(3)
%       't'
%
% Abstract groups can also be created by isomorphisms from an explicit group.
%
% Example:
%   >>> G = replab.S(3);
%   >>> f = G.abstractGroupIsomorphism({'s' 't'});
%   >>> f.imageElement([2 3 1])
%       's'

    properties (Access = protected)
        groupId % (integer): Unique group id
    end

    properties (SetAccess = protected)
        generatorNames % (cell(1,\*) of charstring): Generator names
        permutationGroup % (`.PermutationGroup`): Permutation group realization of this abstract group
        name % (charstring or ``[]``): Group name
    end

    methods (Static)

        function A = make(generatorNames, relatorWords)
            relators = cellfun(@(w) replab.fp.parseLetters(w, generatorNames), relatorWords, 'uniform', 0);
            gens = replab.fp.permutationGeneratorsForRelators(length(generatorNames), relators);
            pg = replab.PermutationGroup.of(gens{:});
            A = replab.AbstractGroup(generatorNames, pg, relatorWords);
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
            [ok, generatorNames, relatorLetters] = replab.fp.Parser.parsePresentation(str);
            assert(ok, 'Error in given presentation string');
            relatorLetters = cellfun(@(r) replab.fp.Letters.reduce(r), relatorLetters, 'uniform', 0);
            mask = cellfun(@isempty, relatorLetters);
            relatorLetters = relatorLetters(~mask);
            relators = cell(1, length(relatorLetters));
            for i = 1:length(relatorLetters)
                relators{i} = replab.fp.Letters.print(relatorLetters{i}, generatorNames, ' ');
            end
            A = replab.AbstractGroup.make(generatorNames, relators);
            if nargout > 1
                for i = 1:length(A.nGenerators)
                    varargout{i} = A.generator(i);
                end
            end
        end

    end

    methods (Access = protected)

        function r = computeRelators(self)
            r = replab.fp.relatorsForPermutationGroup(self.permutationGroup, self.generatorNames);
            r = cellfun(@(w) self.fromLetters(w), r, 'uniform', 0);
        end

        function m = computeNiceMorphism(self)
            if self.order <= 65536
                m = replab.nfg.AbstractGroupIsomorphismEnumeration.make(self);
            else
                m = replab.nfg.AbstractGroupIsomorphismChain(self);
            end
        end

    end

    methods

        function self = AbstractGroup(generatorNames, permutationGroup, relators, name)
        % Creates an abstract group from generator names, images and optional relators
        %
        % Args:
        %   generatorNames (cell(1,\*) of charstring): Generator names
        %   permutationGroup (`.PermutationGroup`): Permutation group realization of this abstract group
        %   relators (cell(1,\*) of charstring, optional): Relators
            self.type = self;
            self.groupId = replab.globals.nextUniqueId;
            self.identity = '1';
            self.generatorNames = generatorNames;
            self.generators = generatorNames;
            self.permutationGroup = permutationGroup;
            if nargin >= 3 && ~isempty(relators)
                self.cache('relators', relators, 'error');
            end
            if nargin >= 4 && ~isempty(name);
                self.name = name;
            else
                self.name = 'Abstract group';
            end
        end

        function A1 = withRenamedGenerators(self, generatorNames1)
        % Returns a modified copy of this abstract group with the generators renamed
        %
        % Args:
        %   generatorNames1 (cell(1,\*) of charstring): New generator names
        %
        % Returns:
        %   `.AbstractGroup`: Updated copy
            A1 = replab.AbstractGroup(generatorNames1, self.permutationGroup, self.cachedOrEmpty('relators'));
            if self.inCache('niceMorphism')
                A1.cache('niceMorphism', self.niceMorphism.withUpdatedSource(A1), 'error');
            end
        end

        function r = relators(self)
        % Returns the list of relators defining this abstract group
        %
        % Returns:
        %   cell(1,\*) of charstring: Words defining the relators
            r = self.cached('relators', @() self.computeRelators);
        end

        function s = presentationString(self)
            r = self.relators;
            gens = strjoin(self.generatorNames, ', ');
            rels = strjoin(self.relators, ' = ');
            s = sprintf('< %s | %s = 1 >', gens, rels);
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
            [ok, tokens] = replab.fp.Parser.lex(word, self.generatorNames);
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
            word = replab.fp.Letters.print(letters, self.generatorNames, ' ');
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
        %   permutation: Computed image
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

        % Str

        function s = shortStr(self, maxColumns)
            if self.inCache('relators')
                s = self.presentationString;
            else
                s = 'AbstractGroup';
            end
        end

        function h = headerStr(self)
            h = self.name;
        end

        % Domain

        function b = eqv(self, x, y)
            b = all(self.niceImage(x) == self.niceImage(y));
        end

        function l = laws(self)
            l = replab.laws.AbstractGroupLaws(self);
        end

        % Monoid

        function z = compose(self, x, y)
            xl = self.toLetters(x);
            yl = self.toLetters(y);
            zl = replab.fp.Letters.compose(xl, yl);
            z = self.fromLetters(zl);
        end

        % Group

        function z = inverse(self, x)
            xl = self.toLetters(x);
            zl = replab.fp.Letters.inverse(xl);
            z = self.fromLetters(zl);
        end

        % FiniteSet

        function b = hasSameTypeAs(self, rhs)
            b = self.type.groupId == rhs.type.groupId;
        end

        % NiceFiniteGroup

        function perm = niceImage(self, word)
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
