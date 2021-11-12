classdef AbstractGroup < replab.gen.FiniteGroup
% Finite group defined using a presentation: a set of generators and a set of relations among generators
%
% The name of generator always start with a letter (a-z or A-Z), and then contain either letters (a-z or A-Z),
% digits (0-9) or underscores. Words in the generators are strings (charstring) containing generators separated
% by spaces, that sequence being understood as a product of generators.
%
% In a word, generators can be taken to an integer power, using the syntax ``x^n`` where ``x`` is the name of a generator
% and ``n`` is an integer. In the word syntax, we optionally accept explicit multiplication operators (``*``),
% division operators (then ``x/y`` is understood as ``x y^-1``), parenthesis, and commutators ``[x, y] = x^-1 y^-1 x y``.
%
% The abstract group is defined by a set of relations, that define equivalence relations between words. For example,
% a cyclic group is defined by the generator ``x`` and the relation ``x^n = 1`` where ``n`` is the order of the group.
% Then we deduce the relations ``x^(n+m) = x^m`` and that ``x^-m = x^(n-m)`` for any integer ``m``.
%
% Informally, the abstract group is the "largest" group with the sets of generators subject to these relations. Note
% that RepLAB rewrites all relations as ``lhs = 1``, and those left hand sides are named "relators".
%
% See `<https://en.wikipedia.org/wiki/Presentation_of_a_group>`_ for a detailed discussion.
%
% Internally, RepLAB computes a realization of the abstract group as a permutation group, and delegates the computations
% to that permutation group.
%
% Example:
%   >>> [G, x] = replab.AbstractGroup.fromPresentation('< x | x^3 = 1 >');
%   >>> x
%       'x'
%   >>> G.compose(x, x)
%       'x^2'
%   >>> G.order
%       3
%
% Example:
%   >>> G = replab.AbstractGroup({'x'}, {'x^3'});
%   >>> G.order
%       3
%
% Abstract groups can also be created by isomorphisms from an explicit group.
%
% Example:
%   >>> G = replab.S(3);
%   >>> f = G.abstractMorphism({'s' 't'});
%   >>> f.imageElement([2 3 1])
%       's'

    properties (SetAccess = protected)
       name % (charstring): Group name
       inAtlas % (logical): Whether this group is already part of the atlas
    end

    methods (Static)

        function [A, varargout] = fromPresentation(str, varargin)
        % Creates an abstract group from a presentation string
        %
        % Returns the finite group generators as additional output arguments.
        %
        % Additional arguments are passed to the `.AbstractGroup` constructor.
        %
        % Example:
        %   >>> [G, x] = replab.AbstractGroup.fromPresentation('< x | x^3 = 1 >');
        %   >>> [G, x, y] = replab.AbstractGroup.fromPresentation('< x, y | x^3 = y^2 = x y x y^-1 = 1 >');
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
            A = replab.AbstractGroup(generatorNames, relators, varargin{:});
            if nargout > 1
                for i = 1:A.nGenerators
                    varargout{i} = A.generator(i);
                end
            end
        end

    end

% $$$     methods (Access = protected)
% $$$
% $$$         function G = computeNiceGroup(self)
% $$$             G = self.permutationGroup;
% $$$         end
% $$$
% $$$         function m = computeNiceMorphism(self)
% $$$             m = replab.mrp.AbstractGroupNiceIsomorphism(self);
% $$$         end
% $$$
% $$$         function A = computeAbstractGroup(self)
% $$$             A = self;
% $$$         end
% $$$
% $$$         function m = computeAbstractMorphism(self)
% $$$             m = replab.FiniteIsomorphism.identity(self);
% $$$         end
% $$$
% $$$         function R = computeRecognize(self)
% $$$             R = [];
% $$$             if ~self.inAtlas
% $$$                 R = replab.Atlas.recognize(self);
% $$$             end
% $$$         end
% $$$
% $$$         function R = computeFastRecognize(self)
% $$$             R = [];
% $$$             if ~self.inAtlas
% $$$                 R = self.niceGroup.fastRecognize;
% $$$                 if ~isempty(R)
% $$$                     R = R.andThen(self.niceMorphism.inverse);
% $$$                 end
% $$$             end
% $$$         end
% $$$
% $$$     end

    methods

        function self = AbstractGroup(generatorNames, relators, varargin)
        % Creates an abstract group from generator names and relators
        %
        % Args:
        %   generatorNames (cell(1,\*) of charstring): Generator names
        %   relators (cell(1,\*) of charstring, optional): Relators
        %
        % Keyword Args:
        %   name (charstring, optional): Group name, optional
        %   order (integer, optional): Group order
        %   permutationGenerators (cell(1,\*) of permutation): Realization of the generators using permutations
        %   inAtlas (logical, optional): Whether this group is part of the atlas
            args = struct('order', 0, 'permutationGenerators', 'none', 'generatorNames', 'none', 'name', 'Abstract group', 'inAtlas', false);
            args = replab.util.populateStruct(args, varargin);
            assert(isequal(args.generatorNames, 'none'), 'Cannot provide generatorNames argument in addition');
            % parse relators into flat form
            relatorsLetter = relators;
            for i = 1:length(relators)
                if ischar(relators{i})
                    relatorsLetter{i} = replab.fp.Letters.parse(relators{i}, generatorNames);
                end
            end
            % Get the target group
            if isequal(args.permutationGenerators, 'none')
                [permutationGenerators, domainSize, order] = replab.fp.permutationRealizationForRelators(length(generatorNames), relatorsLetter);
                nice = replab.PermutationGroup(domainSize, permutationGenerators, 'order', order, 'relators', relators, 'generatorNames', generatorNames);
                if ~isequal(args.order, 0)
                    assert(order == args.order, 'The given order does not match the computed order');
                end
            else % Permutation generators are given
                domainSize = length(args.permutationGenerators{1});
                if isequal(args.order, 0)
                    nice = replab.PermutationGroup(domainSize, args.permutationGenerators, 'relators', relators, 'generatorNames', generatorNames);
                else
                    nice = replab.PermutationGroup(domainSize, args.permutationGenerators, 'relators', relators, 'generatorNames', generatorNames, 'order', args.order);
                end
            end

            type = replab.abstract.FiniteGroupType(generatorNames, nice);

            self@replab.gen.FiniteGroup(type, generatorNames, 'nice', nice, 'niceIsomorphism', type.isomorphism);
            self.type.isomorphism.setSource(self);
            self.name = args.name;
            self.inAtlas = args.inAtlas;
        end

% $$$         function res = simplify(self, word)
% $$$         % Attempts to simplify the given word
% $$$         %
% $$$         % Args:
% $$$         %   word (charstring): Word to simplify
% $$$         %
% $$$         % Returns:
% $$$         %   charstring: Simplified word
% $$$             res = self.niceMorphism.preimageElement(self.niceMorphism.imageElement(word));
% $$$             if length(res) > length(word)
% $$$                 res = word;
% $$$             end
% $$$         end
% $$$
% $$$         function letters = factorizeLetters(self, word)
% $$$         % Parses word letters from word as a string
% $$$         %
% $$$         % Example:
% $$$         %   >>> A = replab.AbstractGroup({'x'}, {'x^3'});
% $$$         %   >>> isequal(A.factorizeLetters('x^2'), [1 1])
% $$$         %       1
% $$$         %
% $$$         % Args:
% $$$         %   word (charstring): Word as a string
% $$$         %
% $$$         % Returns:
% $$$         %   integer(1,\*): Word letters
% $$$         %
% $$$         % Raises:
% $$$         %   An error if the string is malformed
% $$$             [ok, tokens] = replab.fp.Parser.lex(word, self.generatorNames);
% $$$             assert(ok, 'Unknown tokens in string');
% $$$             [pos, letters] = replab.fp.Parser.word(tokens, 1);
% $$$             assert(pos > 0, 'Malformed word');
% $$$             assert(tokens(1, pos) == replab.fp.Parser.types.END, 'Badly terminated word');
% $$$         end
% $$$
% $$$         function word = imageLetters(self, letters)
% $$$         % Prints a word formed of letters as a string
% $$$         %
% $$$         %   >>> A = replab.AbstractGroup({'x'}, {[2 3 1]}, {'x^3'});
% $$$         %   >>> A.imageLetters([1 1])
% $$$         %       'x^2'
% $$$         %
% $$$         % Args:
% $$$         %   letters (integer(1,\*)): Word letters
% $$$         %
% $$$         % Returns:
% $$$         %   charstring: Word as a string
% $$$             word = replab.fp.Letters.print(letters, self.generatorNames, ' ');
% $$$         end
% $$$
% $$$         function img = computeImage(self, word, target, targetGeneratorImages)
% $$$         % Computes the image of this word using the given generator images
% $$$         %
% $$$         % Does not verify the validity of the implied homomorphism.
% $$$         %
% $$$         % Args:
% $$$         %   word (charstring): Word
% $$$         %   target (`+replab.Group`): Target group
% $$$         %   targetGeneratorImages (cell(1,\*) of elements of ``target``): Images of the generators of this group
% $$$         %
% $$$         % Returns:
% $$$         %   permutation: Computed image
% $$$             letters = self.factorizeLetters(word);
% $$$             img = target.identity;
% $$$             for i = 1:length(letters)
% $$$                 l = letters(i);
% $$$                 if l > 0
% $$$                     img = target.compose(img, targetGeneratorImages{l});
% $$$                 else
% $$$                     img = target.composeWithInverse(img, targetGeneratorImages{-l});
% $$$                 end
% $$$             end
% $$$         end
% $$$
% $$$         function m = renamingMorphism(self, newNames)
% $$$         % Returns a morphism from this abstract group with the generators renamed
% $$$         %
% $$$         % Args:
% $$$         %   newNames (cell(1,\*) of charstring): New generator names
% $$$         %
% $$$         % Returns:
% $$$         %   `.AbstractGroup`: Updated copy
% $$$             m = replab.mrp.AbstractGroupRenamingIsomorphism(self, self.withGeneratorNames(newNames));
% $$$         end

    end

    methods % Implementations

        % Str

        function s = shortStr(self, maxColumns)
            s = [self.name ' ' self.presentation];
        end

        function h = headerStr(self)
            h = self.name;
        end

        function names = hiddenFields(self)
            names = hiddenFields@replab.gen.FiniteGroup(self);
            names{1,end+1} = 'type';
            names{1,end+1} = 'inAtlas';
        end
% $$$
% $$$         function [names, values] = additionalFields(self)
% $$$             [names, values] = additionalFields@replab.NiceFiniteGroup(self);
% $$$             names{1,end+1} = 'relators';
% $$$             values{1,end+1} = self.relators;
% $$$         end
% $$$
% $$$         % Domain
% $$$
% $$$         function b = eqv(self, x, y)
% $$$             b = all(self.niceImage(x) == self.niceImage(y));
% $$$         end
% $$$
% $$$         function l = laws(self)
% $$$             l = replab.laws.AbstractGroupLaws(self);
% $$$         end
% $$$
% $$$         % Monoid
% $$$
% $$$         function z = compose(self, x, y)
% $$$             xl = self.factorizeLetters(x);
% $$$             yl = self.factorizeLetters(y);
% $$$             zl = replab.fp.Letters.compose(xl, yl);
% $$$             z = self.imageLetters(zl);
% $$$         end
% $$$
% $$$         % Group
% $$$
% $$$         function z = inverse(self, x)
% $$$             xl = self.factorizeLetters(x);
% $$$             zl = replab.fp.Letters.inverse(xl);
% $$$             z = self.imageLetters(zl);
% $$$         end
% $$$
% $$$         % FiniteSet
% $$$
% $$$         function b = hasSameTypeAs(self, rhs)
% $$$             b = self.type.groupId == rhs.type.groupId;
% $$$         end
% $$$
% $$$         % FiniteGroup
% $$$
% $$$         function A1 = withGeneratorNames(self, newNames)
% $$$         % Returns a modified copy of this abstract group with the generators renamed
% $$$         %
% $$$         % Note: the abstract group returned by this method is not equal to the original abstract group. This differs
% $$$         % from the behavior of `.withGeneratorNames` called on any `.FiniteGroup` which is not an abstract group.
% $$$         %
% $$$         % Args:
% $$$         %   newNames (cell(1,\*) of charstring): New generator names
% $$$         %
% $$$         % Returns:
% $$$         %   `.AbstractGroup`: Updated copy
% $$$             if isequal(self.generatorNames, newNames)
% $$$                 A1 = self;
% $$$                 return
% $$$             end
% $$$             rels = cellfun(@(r) replab.fp.Letters.print(self.factorizeLetters(r), newNames), self.relators, 'uniform', 0);
% $$$             args = {};
% $$$             if self.inCache('order')
% $$$                 args{1,end+1} = 'order';
% $$$                 args{1,end+1}=  self.order;
% $$$             end
% $$$             if self.inCache('permutationGroup')
% $$$                 args{1,end+1} = 'permutationGenerators';
% $$$                 args{1,end+1} = self.permutationGroup.generators;
% $$$             end
% $$$             A1 = replab.AbstractGroup(newNames, rels, args{:});
% $$$         end
% $$$
% $$$         function A = abstractGroup(self, generatorNames)
% $$$             if nargin < 2 || isempty(generatorNames)
% $$$                 generatorNames = self.generatorNames;
% $$$             end
% $$$             A = self.withGeneratorNames(generatorNames);
% $$$         end
% $$$
% $$$         function m = abstractMorphism(self, generatorNames)
% $$$             if nargin < 2 || isempty(generatorNames)
% $$$                 generatorNames = self.generatorNames;
% $$$             end
% $$$             m = self.renamingMorphism(generatorNames);
% $$$         end
% $$$
% $$$         % NiceFiniteGroup
% $$$
% $$$         function perm = niceImage(self, word)
% $$$             letters = self.factorizeLetters(word);
% $$$             pg = self.permutationGroup;
% $$$             perm = pg.identity;
% $$$             for i = 1:length(letters)
% $$$                 l = letters(i);
% $$$                 if l > 0
% $$$                     perm = pg.compose(perm, pg.generator(l));
% $$$                 else
% $$$                     perm = pg.composeWithInverse(perm, pg.generator(-l));
% $$$                 end
% $$$             end
% $$$         end

    end

end
