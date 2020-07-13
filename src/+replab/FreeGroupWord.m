classdef FreeGroupWord < replab.Str
% Describes a word in a free group
%
% In a free group, this word describes a product of the group generators and their inverses.
% By convention, the word is composed of a sequence of letters, stored as an integer row vector.
% The generators have indices $1, 2, ..., n$, while the indices $-1, -2, ..., -n$ describe the
% inverses of those generators.
%
% Also, empty words must have ``reducedLetters = zeros(1, 0)`` and not ``reducedLetters = []``
%
% In a free group, the words are always written in their reduced form. Words can always be compared
% by their letters.

    properties
        group % (`+replab.FreeGroup`): Group this word is part of
        reducedLetters % (integer(1,\*)): Contents of the word, see documentation in `.Word`
    end

    methods (Access = protected) % Protected constructor

        function self = FreeGroupWord(group, reducedLetters)
            assert(isequal(reducedLetters, replab.fp.reduceLetters(reducedLetters)));
            self.group = group;
            self.reducedLetters = reducedLetters;
        end

    end

    methods (Static) % Word construction

        function word = parse(group, string)
        % Constructs a word from a string
        %
        % Example:
        %   >>> [F x y] = replab.FreeGroup.of('x', 'y');
        %   >>> replab.FreeGroupWord.parse(F, 'x (x y)^2 / x')
        %       x^2 y x y x^-1
        %
        % Args:
        %   group (`.FreeGroup`): Free group in which to parse this word
        %   str (charstring): String describing a word
        %
        % Returns:
        %   `.FreeGroupWord`: The parsed word
        %
        % Raises:
        %   An error if the string is malformed
            P = replab.fp.Parser;
            [res, tokens] = P.lex(string, group.generatorNames);
            assert(res, 'Unknown tokens in string');
            pos = 1;
            [pos, letters] = P.word(tokens, pos);
            assert(pos > 0, 'Malformed word');
            assert(tokens(1, pos) == P.types.END, 'Badly terminated word');
            word = replab.FreeGroupWord.make(group, letters);
        end

        function word = make(group, letters)
        % Constructs a word from a sequence of letters
        %
        % The constructed word is reduced.
        %
        % Example:
        %   >>> [F x y] = replab.FreeGroup.of('x', 'y');
        %   >>> replab.FreeGroupWord.make(F, [1 2 -2 1])
        %       x^2
        %
        % Args:
        %   group (`.FreeGroup`): Free group in which to parse this word
        %   letters (integer(1,\*)): Sequence of letters with indices in ``{-r,...,-1,1,...,r}`` where ``r`` is `.rank`
        %
        % Returns:
        %   `.FreeGroupWord`: The reduced word
            rl = replab.fp.reduceLetters(letters);
            word = replab.FreeGroupWord(group, rl);
        end

        function word = empty(group)
        % Constructs the empty word, i.e. the identity
        %
        % Example:
        %   >>> [F x y] = replab.FreeGroup.of('x', 'y');
        %   >>> replab.FreeGroupWord.empty(F).length
        %       0
        %
        % Args:
        %   group (`.FreeGroup`): Free group in which to parse this word
        %
        % Returns:
        %   `.FreeGroupWord`: The empty word
            assert(isa(group, 'replab.FreeGroup'));
            word = replab.FreeGroupWord(group, zeros(1, 0));
        end

    end

    methods % Implementations

        function l = ne(self, rhs)
        % Word non-equality comparison
        %
        % See `.eq`
        %
        % Args:
        %   rhs (`+replab.FreeGroupWord`): Other word
        %
        % Returns:
        %   logical: True if this is not equal to the given word
            l = ~(self == eq);
        end

        function l = eq(self, rhs)
        % Word equality comparison
        %
        % Will return true if
        %
        % - this word is part of the same free group as the other word
        % - their reduced sequences of letters match
        %
        % Args:
        %   rhs (`+replab.FreeGroupWord`): Other word
        %
        % Returns:
        %   logical: True if this is equal to the given word
            if isempty(self.reducedLetters) % TODO: remove
                assert(isequal(self.reducedLetters, zeros(1, 0)));
            end
            if isempty(rhs.reducedLetters) % TODO: remove
                assert(isequal(rhs.reducedLetters, zeros(1, 0)));
            end
            l = (self.group.id == rhs.group.id) && isequal(self.reducedLetters, rhs.reducedLetters);
        end

        function z = mtimes(self, rhs)
        % Word multiplication
        %
        % Args:
        %   rhs (`.FreeGroupWord`): Other word such that ``rhs.group == self.group``
        %
        % Returns:
        %   `.FreeGroupWord`: The product
            assert(self.group == rhs.group, 'Must multiply words of the same group');
            letters = replab.fp.composeLetters(self.reducedLetters, rhs.reducedLetters);
            z = replab.FreeGroupWord(self.group, letters);
        end

        function z = inv(self)
        % Word inverse
        %
        % Returns a word ``iw`` such that ``iw * self = self * iw = identity``.
        %
        % Returns:
        %   `.FreeGroupWord`: The inverse
            z = replab.FreeGroupWord(self.group, replab.fp.inverseLetters(self.reducedLetters));
        end

        function z = mrdivide(self, rhs)
        % Word multiplication by inverse
            letters = replab.fp.composeLetters(self.reducedLetters, replab.fp.inverseLetters(rhs.reducedLetters));
            z = replab.FreeGroupWord(self.group, letters);
        end

        function z = mpower(self, n)
        % Word power
            z = self.group.composeN(self, n);
        end

    end


    methods % FreeGroupWord methods

        function l = isEmpty(self)
        % Returns whether this word is the empty word
        %
        % Returns:
        %   logical: True if this word is empty
            l = isempty(self.reducedLetters);
        end

        function l = length(self)
        % Number of letters composing this word
        %
        % Returns:
        %   integer: Number of letters in this word
            l = length(self.reducedLetters);
        end

        function s = toString(self)
        % Returns a string representation of this word
        %
        % That string can be parsed back by `.FPGroup.parse`.
        %
        % Example:
        %   >>> [F x y] = replab.FreeGroup.of('x', 'y');
        %   >>> w = x*y/x;
        %   >>> s = w.toString;
        %   >>> w1 = F.parse(s);
        %   >>> w == w1
        %       1
        %
        % Returns:
        %   charstring: String representing the word
            if self.length == 0
                s = '1';
                return
            end
            s = '';
            sep = '';
            i = 1;
            word = self.reducedLetters;
            names = self.group.generatorNames;
            while i <= length(word)
                e = sign(word(i)); % exponent
                l = abs(word(i));
                assert(e ~= 0);
                i = i + 1;
                while i <= length(word) && abs(word(i)) == l
                    e = e + sign(word(i));
                    i = i + 1;
                end
                if e ~= 0
                    s = [s sep names{l}];
                    sep = ' ';
                    if e ~= 1
                        s = [s '^' num2str(e)];
                    end
                end
            end
        end

        function s = headerStr(self)
            s = self.toString;
        end

        function s = shortStr(self, maxColumns)
            s = self.toString;
        end

        function s = longStr(self, maxRows, maxColumns)
            s = {self.toString};
        end

        function g = computeImage(self, target, generatorImages)
        % Computes the image of this word using the given generator images
        %
        % Does not verify the validity of the implied homomorphism.
        %
        % Args:
        %   target (`.Group`): Target group
        %   generatorImages (cell(1,\*) of ``target`` elements): Images of the free group generators in the target group
        %
        % Returns:
        %   element of ``target``: Computed image
            letters = self.reducedLetters;
            g = target.identity;
            for i = 1:length(letters)
                l = letters(i);
                if l > 0
                    g = target.compose(g, generatorImages{l});
                else
                    g = target.composeWithInverse(g, generatorImages{-l});
                end
            end
        end

    end

end
