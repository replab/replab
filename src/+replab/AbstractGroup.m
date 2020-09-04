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
%   >>> f = G.abstractMorphism({'s' 't'});
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
        % Creates an abstract group from generator names and relators
        %
        % Args:
        %   generatorNames (cell(1,\*) of charstring): Names of the generators
        %   relatorWords (cell(1,\*) of charstring): Relators of the abstract group
        %
        % Returns:
        %   `.AbstractGroup`: The constructed abstract group
            if length(generatorNames) == 0
                % Trivial group case
                S1 = replab.S(1);
                A = replab.AbstractGroup(cell(1, 0), S1.trivialSubgroup, cell(1, 0));
                return
            end
            relators = cellfun(@(w) replab.fp.Letters.parse(w, generatorNames), relatorWords, 'uniform', 0);
            [gens order] = replab.fp.permutationRealizationForRelators(length(generatorNames), relators);

            ds = length(gens{1});
            isId = cellfun(@(g) isequal(g, 1:ds), gens);
            if any(isId)
                error('One of the generators in the presentation is the identity');
            end
            pg = replab.PermutationGroup.of(gens{:});
            pg.cache('order', order, '==');
            A = replab.AbstractGroup(generatorNames, pg, relatorWords);
            A.cache('order', order, '==');
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
                for i = 1:A.nGenerators
                    varargout{i} = A.generator(i);
                end
            end
        end

    end

    methods (Access = protected)

        function r = computeRelators(self)
            r = replab.fp.relatorsForPermutationGroup(self.permutationGroup);
            r = cellfun(@(w) self.imageLetters(fliplr(w)), r, 'uniform', 0);
        end

    end

    methods (Access = protected)

        function G = computeNiceGroup(self)
            G = self.permutationGroup;
        end

        function m = computeNiceMorphism(self)
            m = replab.mrp.AbstractGroupNiceIsomorphism(self);
        end

        function A = computeDefaultAbstractGroup(self)
            if isequal(self.generatorNames, self.defaultGeneratorNames)
                A = self;
            else
                A = self.withRenamedGenerators(self.defaultGeneratorNames);
            end
        end

        function m = computeDefaultAbstractMorphism(self)
            if isequal(self.generatorNames, self.defaultGeneratorNames)
                m = replab.FiniteIsomorphism.identity(self);
            else
                m = replab.mrp.AbstractGroupRenamingIsomorphism(self, self.abstractGroup);
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
        %   name (charstring): Name of the abstract group
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
            rels = self.cachedOrEmpty('relators');
            if ~isempty(rels)
                rels = cellfun(@(r) replab.fp.Letters.print(self.factorizeLetters(r), generatorNames1), rels, 'uniform', 0);
            end
            A1 = replab.AbstractGroup(generatorNames1, self.permutationGroup, rels);
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
            if isempty(self.relators)
                s = '< | >'; % empty presentation
            else
                gens = strjoin(self.generatorNames, ', ');
                rels = strjoin(self.relators, ' = ');
                s = sprintf('< %s | %s = 1 >', gens, rels);
            end
        end

        function res = simplify(self, word)
        % Attempts to simplify the given word
        %
        % Args:
        %   word (charstring): Word to simplify
        %
        % Returns:
        %   charstring: Simplified word
            res = self.niceMorphism.preimageElement(self.niceMorphism.imageElement(word));
            if length(res) > length(word)
                res = word;
            end
        end


        function letters = factorizeLetters(self, word)
        % Parses word letters from word as a string
        %
        % Example:
        %   >>> A = replab.AbstractGroup({'x'}, {[2 3 1]}, {'x^3'});
        %   >>> isequal(A.factorizeLetters('x^2'), [1 1])
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

        function word = imageLetters(self, letters)
        % Prints a word formed of letters as a string
        %
        %   >>> A = replab.AbstractGroup({'x'}, {[2 3 1]}, {'x^3'});
        %   >>> A.imageLetters([1 1])
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
            letters = self.factorizeLetters(word);
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
            xl = self.factorizeLetters(x);
            yl = self.factorizeLetters(y);
            zl = replab.fp.Letters.compose(xl, yl);
            z = self.imageLetters(zl);
        end

        % Group

        function z = inverse(self, x)
            xl = self.factorizeLetters(x);
            zl = replab.fp.Letters.inverse(xl);
            z = self.imageLetters(zl);
        end

        % FiniteSet

        function b = hasSameTypeAs(self, rhs)
            b = self.type.groupId == rhs.type.groupId;
        end

        % NiceFiniteGroup

        function perm = niceImage(self, word)
            letters = self.factorizeLetters(word);
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

    methods (Access = protected)

        function l = isMorphismByImages_(self, target, preimages, images)
            hasSameGenerators = length(preimages) == self.nGenerators && ...
                all(arrayfun(@(i) self.eqv(preimages{i}, self.generator(i)), 1:self.nGenerators));
            if hasSameGenerators
                nR = length(self.relators);
                for i = 1:nR
                    r = self.factorizeLetters(self.relators{i});
                    g = target.composeLetters(images, r);
                    if ~target.isIdentity(g)
                        l = false;
                        return
                    end
                end
                l = true;
            else
                l = isMorphismByImages_@replab.FiniteGroup(self, target, preimages, images);
            end
        end

    end

end
